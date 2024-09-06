from pathlib import Path
import yaml
import logging


def get_conda_env_dict() -> dict:
    """获取环境字典"""
    logging.info('获取环境字典')
    yaml_conda_env = Path(__file__).resolve().parent.joinpath('config/conda_env.yaml')

    with open(yaml_conda_env) as f:
        dict_conda_env = yaml.safe_load(f)

    return dict_conda_env
