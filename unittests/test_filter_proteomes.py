import unittest
import os
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import tempfile

from pipeline_util.filter_proteomes import extract_gene_name, filter_proteome, write_filtered_proteome


class TestProteomeFilter(unittest.TestCase):
    def setUp(self):
        """Create a temporary FASTA file for testing."""
        self.temp_fasta = tempfile.NamedTemporaryFile(mode="w+", delete=False)
        records = [
            SeqRecord(
                Seq("MKT"),
                id="sp|PROT1|GENE1",
                description="sp|PROT1|GENE1 GN=GENE1 PE=1 SV=1",
            ),
            SeqRecord(
                Seq("MKTA"),
                id="tr|PROT2|GENE1",
                description="tr|PROT2|GENE1 GN=GENE1 PE=2 SV=1",
            ),
            SeqRecord(
                Seq("MK"),
                id="sp|PROT3|GENE2",
                description="sp|PROT3|GENE2 GN=GENE2 PE=1 SV=1",
            ),
            SeqRecord(
                Seq("MKAA"),
                id="tr|PROT4|GENE3",
                description="tr|PROT4|GENE3 GN=GENE3 PE=2 SV=1",
            ),
            SeqRecord(
                Seq("MKAAA"),
                id="tr|PROT5|GENE3",
                description="tr|PROT5|GENE3 GN=GENE3 PE=2 SV=1",
            ),
        ]
        SeqIO.write(records, self.temp_fasta, "fasta")
        self.temp_fasta.close()

        self.output_fasta = tempfile.NamedTemporaryFile(delete=False).name

    def tearDown(self):
        """Clean up temporary files."""
        os.remove(self.temp_fasta.name)
        if os.path.exists(self.output_fasta):
            os.remove(self.output_fasta)

    def test_extract_gene_name(self):
        record = SeqRecord(
            Seq("MKT"),
            id="sp|PROT1|GENE1",
            description="sp|PROT1|GENE1 GN=GENE1 PE=1 SV=1",
        )
        self.assertEqual(extract_gene_name(record), "GENE1")

        record_no_gn = SeqRecord(
            Seq("MKT"), id="sp|PROT6|UNKNOWN", description="sp|PROT6|UNKNOWN PE=1 SV=1"
        )
        self.assertEqual(extract_gene_name(record_no_gn), "sp|PROT6|UNKNOWN")

    def test_filter_proteome(self):
        best_records = filter_proteome(self.temp_fasta.name)
        # GENE1 should keep the reviewed entry
        self.assertIn("GENE1", best_records)
        self.assertTrue(best_records["GENE1"][2])  # reviewed

        # GENE2 has only one reviewed entry
        self.assertIn("GENE2", best_records)
        self.assertEqual(str(best_records["GENE2"][0].seq), "MK")

        # GENE3 should keep the longest unreviewed entry
        self.assertIn("GENE3", best_records)
        self.assertEqual(str(best_records["GENE3"][0].seq), "MKAAA")

    def test_write_filtered_proteome(self):
        best_records = filter_proteome(self.temp_fasta.name)
        write_filtered_proteome(best_records, self.output_fasta)
        # Verify output file exists and contains correct number of sequences
        sequences = list(SeqIO.parse(self.output_fasta, "fasta"))
        self.assertEqual(len(sequences), 3)
        gene_names = [extract_gene_name(r) for r in sequences]
        self.assertIn("GENE1", gene_names)
        self.assertIn("GENE2", gene_names)
        self.assertIn("GENE3", gene_names)


if __name__ == "__main__":
    unittest.main()
