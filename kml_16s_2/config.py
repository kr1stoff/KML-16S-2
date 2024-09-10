from pathlib import Path
import yaml
import logging
import os
import math


def get_conda_env_dict() -> dict:
    """
    获取 Conda 环境字典
    :return:    环境字典
    """
    logging.info('获取环境字典')
    yaml_conda_env = Path(__file__).resolve().parents[1].joinpath('config/conda_env.yaml')

    with open(yaml_conda_env) as f:
        dict_conda_env = yaml.safe_load(f)

    return dict_conda_env


def get_threads_dict() -> dict:
    """
    获取最大线程数, 高线程分配为 max * 2 / 3, 低线程为 high / 2
    :return dict_thr:   high_threads, low_threads 高/低线程分配数
    """
    logging.info('获取线程数字典')
    max_threads = os.cpu_count()
    high_threads = math.ceil(max_threads * 2 / 3)
    low_threads = math.floor(high_threads / 4)
    dict_thr = {'high': high_threads, 'low': low_threads}

    return dict_thr
