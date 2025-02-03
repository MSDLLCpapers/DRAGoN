#!/usr/bin/env python

import argparse
import pathlib
import typing
import re
import shutil
import tempfile
import contextlib
import os
import io
import sys
import fnmatch
from collections import defaultdict

import boto3
if typing.TYPE_CHECKING:
    from boto3_type_annotations.s3 import Object

sys.path.append(str(pathlib.Path(__file__).parent.resolve()))
import bin.dragon_report as dragrep


s3 = boto3.resource('s3')
s3pat = re.compile('s3://(?P<Bucket>[^/]+)/(?P<Key>.+$)')


@contextlib.contextmanager
def working_in(workdir: os.PathLike):
    prev_dir = os.getcwd()
    os.chdir(workdir)
    yield
    os.chdir(prev_dir)


def maybe_s3(path: os.PathLike):
    if m := s3pat.match(path):
        # Must download to BinaryIO but want to return TextIO
        fileobj = io.BytesIO()
        s3.Object(*m.groups()).download_fileobj(fileobj)
        return io.StringIO(fileobj.getvalue().decode())
    return open(path)


class CLI(argparse.Namespace):
    nextflow_log: typing.TextIO
    outdir: pathlib.Path

    _parser = argparse.ArgumentParser()
    _parser.add_argument('nextflow_log', type=maybe_s3)
    _parser.add_argument('outdir', type=pathlib.Path)

    def __init__(self, args=None):
        self.__class__._parser.parse_args(args, self)
        self.to_stage: defaultdict[str, list[tuple[pathlib.Path | 'Object'], str]] = defaultdict(list)
        self.ambiguous = ['Unique']
        self.client = boto3.client('s3')
        self._workdir: tempfile.TemporaryDirectory = None
    
    @property
    def workdir(self):
        return self._workdir and pathlib.Path(self._workdir.name)

    def stage_one(self, workdir: 'pathlib.Path | Object', relsource: str):
        dest = self.workdir / relsource
        dest.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(workdir, pathlib.Path):
            dest.symlink_to(workdir / relsource)
        else:
            s3.Object(workdir.bucket_name, workdir.key + '/' + relsource).download_file(dest)
    
    def parse_log(self):
        patterns: dict[str, tuple[tuple[str, str]]] = {
            'Preprocessing:Demultiplex': (('dmux', '*.demux.json'),),
            'DRAGoN:FeatureCounts': (('features', '*.summary'),),
            'DRAGoN:Deduplicate': (('matrices', 'DRAGoN.out.*'),),
        }
        patterns['DRAGoN:pDRAGoN'] = sum((patterns[f'DRAGoN:{step}'] for step in ('FeatureCounts', 'Deduplicate')), ())
        for line in self.nextflow_log:
            if 'nextflow run' in line:
                if m := re.search('--Dedup.ambiguous ([^-]+)', line):
                    self.ambiguous = m[1].split()
                break
        for match in re.finditer(
            '\[Task monitor\] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler\[(jobId: \d+; )?id: \d+; name: (?P<name>.+?); status: COMPLETED; exit: 0; error: -; workDir: (?P<workdir>.+?)( started: \d+; exited: .+?; )?\]',
            self.nextflow_log.read()
        ):
            jobname = match['name'].split()[0]
            if jobname not in patterns:
                continue
            if match['workdir'].startswith('s3://'):
                bucketname, key = s3pat.match(match['workdir']).groups()
                workobj = s3.Object(bucketname, key)
                resp = self.client.list_objects_v2(Bucket=bucketname, Prefix=key + '/')
                for result in resp['Contents']:
                    obj = result['Key']
                    for arg, pattern in patterns[jobname]:
                        if fnmatch.fnmatchcase(obj, '**/' + pattern):
                            self.to_stage[arg].append((workobj, obj.replace(key + '/', '')))
                            break
            else:
                workdir = pathlib.Path(match['workdir'])
                for arg, pattern in patterns[jobname]:
                    self.to_stage[arg].extend((workdir, str(x.relative_to(workdir))) for x in workdir.glob(pattern))
        if missing_keys := {'dmux', 'matrices', 'features'} - set(self.to_stage):
            raise KeyError(', '.join(missing_keys))

    def stage(self):
        self._workdir = tempfile.TemporaryDirectory()
        self.workdir.mkdir(exist_ok=True)
        for objects in self.to_stage.values():
            for workdir, obj in objects:
                self.stage_one(workdir, obj)
    
    def run(self):
        args = []
        for key, value in self.to_stage.items():
            args.append(f'--{key}')
            args.extend(set(x.split('/')[0] for _, x in value))
        with working_in(self.workdir):
            dragrep.CLI(args + ['--ambiguous'] + self.ambiguous).main()
    
    def unstage(self):
        self._workdir.cleanup()
        self._workdir = None
    
    def publish(self):
        self.outdir.mkdir(parents=True, exist_ok=True)
        output_files = ['DRAGoN.out', 'DRAGoNreport.txt', 'featureCounts.summary', 'dedup_stats.tsv', 'demux.json']
        for file in output_files:
            actual_file = self.workdir / file
            (shutil.copy if actual_file.is_file() else shutil.copytree)(actual_file, self.outdir / file)
    
    def __enter__(self):
        self.stage()
        return self
    
    def __exit__(self, exc_type, exc, exc_traceback):
        self.unstage()
    
    def main(self):
        self.parse_log()
        with self:
            self.run()
            self.publish()


if __name__ == '__main__':
    CLI().main()
