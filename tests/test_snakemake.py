import unittest
from kml_16s_2 import get_names_by_sample_table
from kml_16s_2 import create_snakemake_configfile


class MyTestCase(unittest.TestCase):
    def test_create(self):
        workdir = '/data/mengxf/Project/KML240906_16S_pipeline/result/24090601'
        sample_table = '/data/mengxf/GitHub/KML-16S-2/template/sample_table.xlsx'
        sample_names = get_names_by_sample_table(sample_table)
        create_snakemake_configfile(sample_names, workdir)
        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
