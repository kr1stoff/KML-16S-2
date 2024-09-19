#!/usr/bin/env python
"""
@Time ： 2024/9/6
@Auth ： kristoff
@IDE ：PyCharm
@Motto：Continuous learning
@LastModified : 2024/9/19
"""

import click
import logging
from pathlib import Path
from kml_16s_2 import prepare_fastq_by_sample_table
from kml_16s_2 import generata_metadata_table
from kml_16s_2 import get_names_by_sample_table
from kml_16s_2 import create_snakemake_configfile
from kml_16s_2 import run_snakemake

logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


@click.command()
@click.option('-i', '--sample_table', required=True, type=click.Path(exists=True),
              help='输入样本表格. 四列 "name group fq1 fq2".')
@click.option('-o', '--workdir', default='kml_16s_2_out', help='结果输出目录.')
@click.help_option('-h', '--help')
def main(sample_table, workdir):
    """KML DJH 16S 分析流程"""
    logging.info('开始分析')

    prepare_fastq_by_sample_table(workdir, sample_table)
    generata_metadata_table(workdir, sample_table)
    sample_names = get_names_by_sample_table(sample_table)
    create_snakemake_configfile(sample_names, workdir)
    run_snakemake(workdir)

    logging.info('分析结束')


if __name__ == '__main__':
    main()
