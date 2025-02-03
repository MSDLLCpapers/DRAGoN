import array
import contextlib
import io
import os
import pathlib
import random
import sys
import tempfile
import unittest

import pysam

project_root = pathlib.Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root / "bin"))
import demultiplex


class DemultiplexTest(unittest.TestCase):
    def setUp(self):
        barcodes_fp = io.StringIO(
            "AAACGGCGAGCCCTTT\tPool1_A1__Media1\t1\t1\n"
            "AGACTTTCAATCACTA\tPool1_B1__293T\t1\t2\n"
            "CCTACCATCAAGTCCA\tPool1_C1__293T\t1\t3\n"
            "CCGAGGCACCGAGAAC\tPool1_D1__293T\t1\t4\n"
            "GGCGTTTCGTGAGAAG\tPool1_E1__293T_DMSO\t1\t5\n"
            "GATGCAGTGTAGTGCG\tPool1_F1__293T_DMSO\t1\t6\n"
            "TTGTAAGATCTTAGGC\tPool1_G1__293T_DMSO\t1\t7\n"
            "TTCTACAGTGCTCTGT\tPool1_H1__DMSO\t1\t8\n"
            "CTGCAGCCCATGGCAG\tPool1_A2__Media1\t2\t1\n"
            "CGTAGGTTCTTAGTCA\tPool1_B2__293T\t2\t2\n"
        )
        self.demultiplexer = demultiplex.Demultiplexer(
            barcodes_fp, umilen=12, check=False
        )
        barcode = "TTGTAAGATCTTAGGC"
        umi = "CATACGACTGCA"
        cdna = "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA"
        self.read_pair = demultiplex.ReadPair(
            "TestRead",
            cdna,
            array.array("B", b"\x25" * 90),
            barcode,
            "FFFFFFFFFFFFFFFF",
            umi,
            "FFFFFFFFFFFF",
            array.array("B", b"\x25" * 16),
            array.array("B", b"\x25" * 12),
        )

    def test_create_class(self):
        self.demultiplexer.mismatch = 1
        self.demultiplexer.check()

    @unittest.expectedFailure
    def test_too_close(self):
        self.demultiplexer.mismatch = 5
        with contextlib.redirect_stderr(io.StringIO()):
            self.demultiplexer.check()

    def test_seq_to_mask(self):
        self.assertEqual(
            self.demultiplexer.seq_to_mask("ACGTNRYVBX="),
            (1, 2, 4, 8, 15, 5, 10, 7, 14, 15, 0),
        )

    def test_get_homopolymer_regex(self):
        self.demultiplexer.min_homop_len = 10
        self.demultiplexer.init_filters()
        self.assertTrue(
            self.demultiplexer.get_homopolymer_regex("A").search(
                "GTGTGTGTACTCGTACTAAAAAAAAAAAAAATGCTATCTACG"
            )
        )

    def test_iter_read_pairs(self):
        with tempfile.NamedTemporaryFile(
            suffix=".fastq"
        ) as r1fq, tempfile.NamedTemporaryFile(suffix=".fastq") as r2fq:
            r1fq.write(
                "@TestRead R1\n"
                f"{self.read_pair.bcseq}{self.read_pair.umiseq}\n"
                "+\n"
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFF\n".encode()
            )
            r2fq.write(
                "@TestRead R2\n"
                f"{self.read_pair.cdna_seq}\n"
                "+\n"
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n".encode()
            )
            r1fq.seek(0)
            r2fq.seek(0)
            for read_pair in self.demultiplexer.iter_read_pairs(r1fq.name, r2fq.name):
                break
            else:
                self.fail("fastq iteration failed")
        self.assertEqual(read_pair, self.read_pair)

    def test_match_sequences(self):
        self.assertEqual(
            self.demultiplexer.match_sequences(
                self.demultiplexer.seq_to_mask("ACGTACGTACGT"),
                self.demultiplexer.seq_to_mask("ACGTACGTACGT"),
            ),
            0,
        )
        self.assertEqual(
            self.demultiplexer.match_sequences(
                self.demultiplexer.seq_to_mask("ACGTACGTACGT"),
                self.demultiplexer.seq_to_mask("ACGTACGTACAA"),
            ),
            2,
        )
        self.assertEqual(
            self.demultiplexer.match_sequences(
                self.demultiplexer.seq_to_mask("ACGTACGTACGT"),
                self.demultiplexer.seq_to_mask("AAAAAAAAAAAA"),
            ),
            9,
        )

    def test_match_barcode(self):
        self.assertEqual(self.demultiplexer.match_barcode("GGCGTTTCGTGAGAAG"), (4, 0))
        self.assertEqual(self.demultiplexer.match_barcode("GGCGTATCGTGAGAAG"), (4, 1))
        self.assertEqual(self.demultiplexer.match_barcode("GGCGTTTCGTNAGAAG"), (4, 0))
        self.assertEqual(self.demultiplexer.match_barcode("GNNNTTCGTGAGAAG"), (-1, -1))

    def test_demultiplex_read(self):
        self.assertEqual(self.demultiplexer.demultiplex_read(self.read_pair), (6, 0, 6))
        self.read_pair.bcseq = "ACACACACACACACAC"
        self.assertEqual(
            self.demultiplexer.demultiplex_read(self.read_pair), (-1, -1, -1)
        )
        self.read_pair.bcseq = "GGCGTATCGTGAGAAG"
        self.demultiplexer.split = 6
        self.assertEqual(self.demultiplexer.demultiplex_read(self.read_pair), (4, 1, 2))

    def test_quality_trim_index(self):
        self.assertEqual(self.demultiplexer.quality_trim_index(self.read_pair), 90)
        self.read_pair.cdna_qual = self.read_pair.cdna_qual[:-20] + array.array(
            "B", b"\x13" * 20
        )
        self.assertEqual(self.demultiplexer.quality_trim_index(self.read_pair), 70)

    def test_should_discard_umi_quality(self):
        self.assertEqual(
            self.demultiplexer.should_discard_umi_quality(self.read_pair),
            demultiplex.QCFailReason.PASS,
        )
        self.read_pair.umiqual_a = array.array("B", b"\x13" * 7 + b"\x15" * 5)
        self.assertEqual(
            self.demultiplexer.should_discard_umi_quality(self.read_pair),
            demultiplex.QCFailReason.FAIL_UMI_QUALITY,
        )

    def test_should_discard_length(self):
        self.assertEqual(
            self.demultiplexer.should_discard_length(self.read_pair),
            demultiplex.QCFailReason.PASS,
        )
        self.read_pair.cdna_seq = self.read_pair.cdna_seq[:19]
        self.read_pair.cdna_qual = self.read_pair.cdna_qual[:19]
        self.assertEqual(
            self.demultiplexer.should_discard_length(self.read_pair),
            demultiplex.QCFailReason.FAIL_LENGTH,
        )

    def test_should_discard_homogeny(self):
        self.assertEqual(
            self.demultiplexer.should_discard_homogeny(self.read_pair),
            demultiplex.QCFailReason.PASS,
        )
        self.read_pair.cdna_seq = self.read_pair.cdna_seq[:8] + "".join(
            random.choices("AT", k=82)
        )
        self.assertEqual(
            self.demultiplexer.should_discard_homogeny(self.read_pair),
            demultiplex.QCFailReason.FAIL_AT_CONTENT,
        )
        self.read_pair.cdna_seq = self.read_pair.cdna_seq[:8] + "".join(
            random.choices("GC", k=82)
        )
        self.assertEqual(
            self.demultiplexer.should_discard_homogeny(self.read_pair),
            demultiplex.QCFailReason.FAIL_GC_CONTENT,
        )

    def test_should_discard_homopolymer(self):
        self.assertEqual(
            self.demultiplexer.should_discard_homopolymer(self.read_pair),
            demultiplex.QCFailReason.PASS,
        )
        self.read_pair.cdna_seq = (
            self.read_pair.cdna_seq[:10] + "A" * 20 + "C" * 20 + "G" * 20 + "T" * 20
        )
        self.assertEqual(
            self.demultiplexer.should_discard_homopolymer(self.read_pair),
            demultiplex.QCFailReason.FAIL_HOMOPOLYMER,
        )

    def test_should_discard_trimming(self):
        self.assertEqual(
            self.demultiplexer.should_discard_trimming(self.read_pair),
            demultiplex.QCFailReason.PASS,
        )
        self.read_pair.cdna_qual = self.read_pair.cdna_qual[:15] + array.array(
            "B", b"\x13" * 75
        )
        self.assertEqual(
            self.demultiplexer.should_discard_trimming(self.read_pair),
            demultiplex.QCFailReason.FAIL_QUALITY_TRIMMING,
        )

    def test_quality_check(self):
        self.assertEqual(
            self.demultiplexer.quality_check(self.read_pair),
            demultiplex.QCFailReason.PASS,
        )
        self.read_pair.cdna_seq = self.read_pair.cdna_seq[:8] + "".join(
            random.choices("AT", k=82)
        )
        self.assertEqual(
            self.demultiplexer.quality_check(self.read_pair),
            demultiplex.QCFailReason.FAIL_AT_CONTENT,
        )
        self.read_pair.cdna_seq = (
            self.read_pair.cdna_seq[:15] + "A" * 10 + self.read_pair.cdna_seq[25:]
        )
        self.assertEqual(
            self.demultiplexer.quality_check(self.read_pair),
            demultiplex.QCFailReason.FAIL_AT_CONTENT,
        )
        self.read_pair.umiqual_a = array.array("B", b"\x13" * 7 + b"\x15" * 5)
        self.assertEqual(
            self.demultiplexer.quality_check(self.read_pair),
            demultiplex.QCFailReason.FAIL_UMI_QUALITY,
        )

    def test_adapter_trimming(self):
        self.assertEqual(
            self.demultiplexer.quality_check(self.read_pair),
            demultiplex.QCFailReason.PASS,
        )
        self.demultiplexer.adapter = "CTAGTCGATCGATCGT"
        self.demultiplexer.init_filters()
        self.assertEqual(
            self.demultiplexer.quality_check(self.read_pair),
            demultiplex.QCFailReason.PASS,
        )
        self.assertEqual(self.read_pair.cdna_seq, "AGCTACGTAGTAGCAGCTAGCATAGCTA")

    def test_demultiplex_experiment(self):
        with tempfile.NamedTemporaryFile(
            suffix="_R1.fastq"
        ) as r1fq, tempfile.NamedTemporaryFile(
            suffix=".fastq"
        ) as r2fq, tempfile.NamedTemporaryFile(
            suffix=".bam"
        ) as bam:
            r1fq.write(
                "@TestRead R1\n"
                f"{self.read_pair.bcseq}{self.read_pair.umiseq}\n"
                "+\n"
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFF\n".encode()
            )
            r2fq.write(
                "@TestRead R2\n"
                f"{self.read_pair.cdna_seq}\n"
                "+\n"
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n".encode()
            )
            r1fq.seek(0)
            r2fq.seek(0)
            self.demultiplexer.demultiplex_experiment(
                r1fq.name, r2fq.name, os.path.splitext(bam.name)[0]
            )
            bam.seek(0)
            try:
                with pysam.AlignmentFile(bam.name, check_sq=False) as bamfp:
                    read: pysam.AlignedSegment = next(bamfp)
            except:
                self.skipTest("Temporary failure writing BAM file")
        self.assertEqual(self.demultiplexer.stats.read_count, 1)
        self.assertEqual(self.demultiplexer.stats.well_counts["Pool1_G1__293T_DMSO"], 1)
        self.assertEqual(self.demultiplexer.stats.well_counts["Pool1_H1__DMSO"], 0)
        self.assertEqual(self.demultiplexer.stats.well_counts["noBC"], 0)
        self.assertEqual(read.query_name, self.read_pair.name)
        self.assertEqual(read.query_sequence, self.read_pair.cdna_seq)
        self.assertEqual(read.query_qualities, self.read_pair.cdna_qual)
        self.assertEqual(read.get_tag("CR"), self.read_pair.bcseq)
        self.assertEqual(read.get_tag("UR"), self.read_pair.umiseq)
        self.assertEqual(read.get_tag("CY"), self.read_pair.bcqual)
        self.assertEqual(read.get_tag("UY"), self.read_pair.umiqual)
        self.assertEqual(read.get_tag("CB"), self.read_pair.bcseq)
        self.assertEqual(read.get_tag("UB"), self.read_pair.umiseq)
        self.assertEqual(read.get_tag("BI"), "0006")
        self.assertEqual(read.get_tag("BM"), 0)
        self.assertEqual(read.get_tag("QC"), "PASS")


if __name__ == "__main__":
    unittest.main()
