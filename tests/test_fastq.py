import unittest
from pathlib import Path
from kml_16s_2 import prepare_fastq_by_sample_table, get_names_by_sample_table


class MyTestCase(unittest.TestCase):
    def test_prepare(self):
        workdir = '/data/mengxf/Project/KML240906_16S_pipeline/result/24090601'
        sample_table = '/data/mengxf/GitHub/KML-16S-2/template/sample_table.xlsx'
        prepare_fastq_by_sample_table(workdir, sample_table)
        self.assertEqual(
            Path('/data/mengxf/Project/KML240906_16S_pipeline/result/24090601/.rawdata/KTND240173_1.fastq.gz').exists(),
            True)  # add assertion here

    def test_get(self):
        sample_table = '/data/mengxf/GitHub/KML-16S-2/template/sample_table.xlsx'
        names = get_names_by_sample_table(sample_table)
        self.assertEqual(names != [], True)


if __name__ == '__main__':
    unittest.main()
