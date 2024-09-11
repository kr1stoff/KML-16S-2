import unittest
from pathlib import Path
from kml_16s_2 import prepare_fastq_by_sample_table, get_names_by_sample_table, generata_metadata_table

workdir = '/data/mengxf/Project/KML240906_16S_pipeline/result/24090601'
sample_table = '/data/mengxf/GitHub/KML-16S-2/template/sample_table.xlsx'


class MyTestCase(unittest.TestCase):
    # def test_prepare(self):
    #     prepare_fastq_by_sample_table(workdir, sample_table)
    #     self.assertEqual(True, True)  # add assertion here

    def test_get(self):
        names = get_names_by_sample_table(sample_table)
        self.assertEqual(names != [], True)

    def test_metadata(self):
        generata_metadata_table(workdir, sample_table)
        self.assertEqual(True, True)


if __name__ == '__main__':
    unittest.main()
