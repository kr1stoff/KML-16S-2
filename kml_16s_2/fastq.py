from pathlib import Path
import pandas as pd
from subprocess import run
import re
import logging


def prepare_fastq_by_sample_table(workdir, sample_table: str) -> None:
    """
    在项目目录下面准备 fastq 文件
    - 如果未压缩 link, 如果压缩 zcat
    - 支持 .tsv 和 .xlsx 格式
    :param workdir:         工作目录
    :param sample_table:     样本信息表 SampleSheet
    :return:
    """
    logging.info('在项目目录下面准备 fastq 文件')
    # 创建 {workdir}/.rawdata
    Path(workdir).joinpath('.rawdata').mkdir(exist_ok=True, parents=True)

    df = sample_table_to_dataframe(sample_table)

    # 软链接或解压
    for row in df.iterrows():
        name, fastq1, fastq2 = row[1]
        copy_fastq(workdir, name, fastq1, fastq2)


def get_names_by_sample_table(sample_table: str) -> list:
    """
    获取样本名列表
    :param sample_table:
    :return sample_names:   样本名列表
    """
    logging.info('获取样本名列表')
    df = sample_table_to_dataframe(sample_table)
    return df.iloc[:, 0].to_list()


def sample_table_to_dataframe(sample_table: str) -> pd.DataFrame:
    """
    输入 SampleSheet 转成 DataFrame 格式

    :param sample_table:
    :return df: SampleSheet 转的 DataFrame
    """
    if sample_table.endswith('.xlsx'):
        df = pd.read_excel(sample_table, header=None)
    elif sample_table.endswith('.tsv'):
        df = pd.read_table(sample_table, sep='\t', header=None)
    else:
        raise ValueError(f'sample_table 扩展名必须是 .xlsx or .tsv : {sample_table}')

    # 检查 SampleSheet
    check_sample_table(df)

    return df


def copy_fastq(workdir, name, fq1: str, fq2) -> None:
    """
    复制或压缩 fastq 到 .rawdata 目, 按照指定格式明明
    :param workdir:     分析解雇目录
    :param name:        样本名
    :param fq1:         fastq1
    :param fq2:         fastq2
    """
    if fq1.endswith('.gz'):
        cml = f"""
        cp {fq1} {workdir}/.rawdata/{name}_1.fastq.gz
        cp {fq2} {workdir}/.rawdata/{name}_2.fastq.gz
        """
    else:
        cml = f"""
        gzip -c {fq1} > {workdir}/.rawdata/{name}_1.fastq.gz
        gzip -c {fq2} > {workdir}/.rawdata/{name}_2.fastq.gz
        """
    logging.debug(cml)
    run(cml, shell=True, executable='/bin/bash', capture_output=True)


def check_sample_table(df) -> None:
    """
    检查 SampleSheet 文件, 输入 SampleSheet 转的 DataFrame
    :param df:
    """
    for row in df.iterrows():
        name, fastq1, fastq2 = row[1]

        # 检查名称
        pattern = r'[\\/:*?"<>| ]'
        assert not re.search(pattern, name), f'样本名称含有非法字符 (\\/:*?"<>| ) : {name}'

        # 检查 fastq 是否存在
        assert Path(fastq1).exists(), f'fastq1 不存在 : {fastq1}'
        assert Path(fastq2).exists(), f'fastq2 不存在 : {fastq2}'

        # 检查 fastq1 和 fastq2 是否相同
        assert fastq1 != fastq2, f'fastq1 和 fastq2 相同 : {fastq1} - {fastq2}'
