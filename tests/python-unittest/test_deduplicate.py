import os
import pathlib
import shutil
import sys
import unittest

import numpy as np
import scipy.io as scio
import scipy.sparse as sp

project_root = pathlib.Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root / "bin"))
import collapse_umis


def all_close(A, B, atol=1e-6):
    if np.array_equal(A.shape, B.shape) == 0:
        return False
    diff = A - B
    _, _, values = sp.find(diff)
    return values.size == 0 or np.allclose(values, 0, atol=atol)


class DeduplicateTest(unittest.TestCase):
    def test_deduplicate_2000(self):
        test_root = project_root / "src/test-data"
        self.assertTrue(os.path.exists(test_root / "tmp1.bam"))
        self.assertTrue(os.path.exists(test_root / "tmp1_MarkDups_ref.bam"))
        self.assertTrue(os.path.exists(test_root / "genes.gtf"))
        self.assertTrue(os.path.exists(test_root / "Simulated_1_barcodes.txt"))
        self.assertTrue(os.path.exists(test_root / "DRAGoN.out.ref"))
        try:
            os.remove(test_root / "tmp1_MarkDups.bam")
        except FileNotFoundError:
            pass
        try:
            shutil.rmtree(test_root / "DRAGoN.out")
        except FileNotFoundError:
            pass
        with collapse_umis.Deduplicator(
            test_root / "tmp1.bam",
            1,
            test_root / "Simulated_1_barcodes.txt",
            test_root / "genes.gtf",
            umilen=12,
            mmap_strategies=frozenset(collapse_umis.MultiMapStrategy),
        ) as worker:
            worker.get_umis()
        worker.dump_matrix(test_root / "DRAGoN.out")

        # Compare unique counts matrices
        cmpmat = scio.mmread(test_root / "DRAGoN.out/matrix.mtx")
        refmat = scio.mmread(test_root / "DRAGoN.out.ref/matrix.mtx")
        self.assertTrue(all_close(cmpmat, refmat))
        for strategy in worker.mmap_strategies:
            if strategy is not collapse_umis.MultiMapStrategy.UNIQUE:
                cmpmat = scio.mmread(
                    test_root / f"DRAGoN.out/UniqueAndMult-{strategy.value}.mtx"
                )
                refmat = scio.mmread(
                    test_root / f"DRAGoN.out.ref/UniqueAndMult-{strategy.value}.mtx"
                )
                self.assertTrue(all_close(cmpmat, refmat))


if __name__ == "__main__":
    unittest.main()
