import yaml
from pathlib import Path
from subprocess import run
import logging
from kml_16s_2 import get_conda_env_dict
from kml_16s_2 import get_threads_dict


def create_snakemake_configfile(sample_names, workdir):
    """
    创建 snakemake 配置文件
    :param sample_names:    样本名列表
    :param workdir:         分析结果目录
    :return:
    """
    logging.info('创建 snakemake 配置文件')
    dir_temp = Path(workdir).joinpath('.temp')
    dir_temp.mkdir(exist_ok=True, parents=True)

    dict_smk = {
        'workdir': workdir,
        'samples': sample_names,
        'threads': get_threads_dict(),
        'conda': get_conda_env_dict(),
    }

    with open(f'{workdir}/.temp/snakemake.yaml', 'w') as f:
        yaml.dump(dict_smk, f)
