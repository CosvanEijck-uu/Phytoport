import unittest
import pandas as pd
from io import StringIO
from pipeline_util.fetch_mnemonic_genes import retrieve_genes


class TestRetrieveGenes(unittest.TestCase):
    def setUp(self):
        # Sample TSV content
        self.tsv_data = """gene1|GEN1_HUMAN
gene2|GEN2_MOUSE
gene3|GEN3_HUMAN
gene4|GEN4_YEAST
gene5|GEN5_MOUSE
"""

        self.tsv_file = StringIO(self.tsv_data)

    def test_retrieve_genes(self):
        # Pandas can read from StringIO object as if it were a file
        df = pd.read_csv(self.tsv_file, sep="\t", header=None)
        df.to_csv("temp.tsv", sep="\t", header=False, index=False)

        genes1, genes2 = retrieve_genes("temp.tsv", "HUMAN", "MOUSE")

        self.assertEqual(sorted(genes1), ["GEN1_HUMAN", "GEN3_HUMAN"])
        self.assertEqual(sorted(genes2), ["GEN2_MOUSE", "GEN5_MOUSE"])

    def tearDown(self):
        import os
        if os.path.exists("temp.tsv"):
            os.remove("temp.tsv")


if __name__ == "__main__":
    unittest.main()
