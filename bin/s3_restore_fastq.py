#!/usr/bin/env python

"""
Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
This file is part of DRAGoN.

This source code is licensed under the MIT License found in the
LICENSE file in the root directory of this source tree.
"""

import argparse
import asyncio
import collections.abc
import functools
import logging
import os
import sys
import typing

import aioboto3
from urllib3.util import parse_url

if typing.TYPE_CHECKING:
    from types_aiobotocore_s3 import S3Client
    from types_aiobotocore_s3.literals import StorageClassType
    from types_aiobotocore_s3.type_defs import (
        HeadObjectOutputTypeDef,
        RestoreRequestTypeDef,
    )

logger = logging.getLogger("s3_restore_fastq")

# Typing generics for function annotation later on
_A = typing.TypeVar("_A")
_K = typing.TypeVar("_K")
_R = typing.TypeVar("_R")
_T = typing.TypeVar("_T")


class S3RestoreException(Exception):
    pass


class CLI(argparse.Namespace):
    dry_run: bool = False
    _threads: int = 1
    debug: int = logging.INFO
    fastq_filenames: list[str]

    def __init__(self, args=None):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "-n", "--dry-run", dest="dry_run", action="store_true", default=CLI.dry_run
        )
        parser.add_argument(
            "-t",
            "--threads",
            dest="_threads",
            type=int,
            default=CLI._threads,
        )
        parser.add_argument(
            "-d",
            "--debug",
            dest="debug",
            action="store_const",
            const=logging.DEBUG,
            default=CLI.debug,
        )
        parser.add_argument("fastq_filenames", nargs="+")
        parser.parse_args(args, self)
        logger.setLevel(self.debug)
        self.to_write = []

    @functools.cached_property
    def threads(self):
        return asyncio.Semaphore(self._threads)

    def decompose_s3_path(self, path: str) -> tuple[str, str]:
        scheme, _, bucket, _, key, *_ = parse_url(path)
        assert scheme == "s3"
        return bucket, key.lstrip("/")

    def is_storage_class_archive(self, storage: "StorageClassType | None") -> bool:
        return storage in {"DEEP_ARCHIVE", "GLACIER"}

    def is_intelligent_tiering_archive(self, head: "HeadObjectOutputTypeDef") -> bool:
        return head.get("StorageClass") == "INTELLIGENT_TIERING" and head.get(
            "ArchiveStatus"
        ) in {"DEEP_ARCHIVE_ACCESS", "ARCHIVE_ACCESS"}

    def classify_storage_class(
        self, head: "HeadObjectOutputTypeDef"
    ) -> tuple[bool, bool]:
        return self.is_storage_class_archive(
            head.get("StorageClass")
        ), self.is_intelligent_tiering_archive(head)

    async def is_restore_in_progress(
        self,
        client: "S3Client",
        bucket: str,
        key: str,
        *,
        head: "HeadObjectOutputTypeDef | None" = None,
    ) -> tuple[str, str] | None:
        if head is None:
            logger.debug("await client.head_object(Bucket=%r, Key=%r)", bucket, key)
            head = await client.head_object(Bucket=bucket, Key=key)
        if head.get("Restore", "").startswith('ongoing-request="true"'):
            return bucket, key

    async def should_restore_file(
        self, client: "S3Client", filename: str
    ) -> tuple[str, str] | None:
        try:
            bucket, key = self.decompose_s3_path(filename)
            logger.debug("await client.head_object(Bucket=%r, Key=%r)", bucket, key)
            head: "HeadObjectOutputTypeDef" = await client.head_object(
                Bucket=bucket, Key=key
            )
            is_archive, is_intelligent = self.classify_storage_class(head)
            logger.info(
                "is_archive: %s, is_intelligent: %s", is_archive, is_intelligent
            )
            if is_archive or is_intelligent:
                if not await self.is_restore_in_progress(
                    client, bucket, key, head=head
                ):
                    restore_request: "RestoreRequestTypeDef" = (
                        {}
                        if is_intelligent
                        else {"Days": 3, "GlacierJobParameters": {"Tier": "Standard"}}
                    )
                    if not self.dry_run:
                        logger.debug(
                            "await client.restore_object(Bucket=%r, Key=%r, RestoreRequest=%r)",
                            bucket,
                            key,
                            restore_request,
                        )
                        await client.restore_object(
                            Bucket=bucket, Key=key, RestoreRequest=restore_request
                        )
                return bucket, key
        except Exception as e:
            raise S3RestoreException(
                f"Error checking status of {filename}: {e}"
            ) from None

    async def map_with_concurrency(
        self,
        coro_func: collections.abc.Callable[
            [_T, _A, _K], collections.abc.Coroutine[typing.Any, typing.Any, _R]
        ],
        iterable: collections.abc.Iterable[_T],
        *args: _A,
        return_exceptions=False,
        **kwargs: _K,
    ) -> list[_R]:
        """Applies the coroutine function to each element of `iterable` in turn and runs them concurrently, respecting the user-specified --threads.

        Args:
            coro_func (collections.abc.Callable[ [_T, _A, _K], collections.abc.Coroutine[typing.Any, typing.Any, _R] ]): The coroutine to apply
            iterable (collections.abc.Iterable[_T]): The iterable to apply to
            args: Additional args to pass to each instance of coro_func
            return_exceptions (bool, optional): If True, will return exceptions in the gathered list. Otherwise, will gather them at the end and raise them as an ExceptionGroup. Defaults to False.
            kwargs: Additional kwargs to pass to each instance of coro_func

        Returns:
            list[_R]: _description_
        """

        async def with_concurrency(coro):
            async with self.threads:
                return await coro

        result: list[_R | BaseException] = await asyncio.gather(
            *(with_concurrency(coro_func(x, *args, **kwargs)) for x in iterable),
            return_exceptions=True,
        )
        if not return_exceptions and (
            exceptions := [x for x in result if isinstance(x, BaseException)]
        ):
            logging.critical("One or more jobs failed")
            for e in exceptions:
                logging.critical(None, exc_info=e)
            sys.exit(1)
        return result

    async def main(self):
        if not (
            to_test := [
                file for file in self.fastq_filenames if file.startswith("s3://")
            ]
        ):
            logger.info("Nothing to do")
            return 0

        # Troubleshoot credentials acquisition
        session = aioboto3.Session()
        if not (credentials := await session.get_credentials()):
            logger.error(
                'Unable to locate credentials. You can configure credentials by running "aws configure".'
            )
            return 1

        frz = await credentials.get_frozen_credentials()
        if (
            env_access_key := os.getenv("AWS_ACCESS_KEY_ID")
        ) is not None and env_access_key != frz.access_key:
            logger.error(
                "Found credentials but they do not match the configured environment. Exiting now.",
            )
            return 1

        async def inner1(fn: str):
            return await self.should_restore_file(client, fn)

        async def inner2(tup: tuple[str, str]):
            return await self.is_restore_in_progress(client, *tup)

        async with session.client("s3") as client:
            client: "S3Client"

            logger.info("Checking whether fastq files need to be restored ...")
            restoring = set(
                tup
                for tup in await self.map_with_concurrency(
                    inner1,
                    to_test,
                    return_exceptions=False,
                )
                if tup
            )
            if not restoring:
                logger.info(
                    "No fastq files found on S3 archive tier and all rows valid"
                )
                return 0

            if self.dry_run:
                logger.info("Dry run, exiting")
                return 0

            logger.info(
                "Restoring %d objects. This can take up to 12 hours. Better grab a coffee.",
                len(restoring),
            )

        while restoring:
            logger.info("Waiting to restore %d objects ...", len(restoring))
            await asyncio.sleep(600)  # check every 10 minutes

            async with session.client("s3") as client:
                client: "S3Client"

                restoring = set(
                    tup
                    for tup in await self.map_with_concurrency(
                        inner2, restoring, return_exceptions=False
                    )
                    if tup
                )

        return 0


if __name__ == "__main__":
    logging.basicConfig(format="[%(asctime)s] %(name)s: %(levelname)s: %(message)s")
    cli = CLI()
    try:
        sys.exit(asyncio.run(cli.main()))
    except Exception as e:
        sys.tracebacklimit = 0
        logger.critical("One or more errors occurred", exc_info=e)
        logger.critical("Check the session credentials maybe?")
        sys.exit(1)
