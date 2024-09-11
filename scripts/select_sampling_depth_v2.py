"""
#athour: dengjunhao
#datetime: 2024-03-29
简介：
该脚本用于16S分析的ASV/OTU的序列计数特征表(丰度表需要暂不支持, 需要微调代码)中计算达到特定比例的ASV/OTU时所需要的抽样深度。

目的：
用于分析计算多个样本选取合适的抽平深度

原理:
1. 根据OTU数目表, 构建一个和总reads数目一样长度的数组表
2. 用二分法选取抽样深度点, 从1中构建的数组中进行无放回抽样, 统计OTU种类:
    2.1 随机抽样取多次(默认10次)的平均值
    2.2 当OTU计数低于所需要的百分比深度, 用二分法增加抽样深度, 并且把最低抽样范围定为此值;
    2.3 当OTU计数高于所需要的百分比深度, 用二分法减少抽样深度, 并且把最大抽样范围定为此值;
    2.4 不断重复上述3个步骤, OTU计数接近所需要的百分比深度(默认的差值小于0.0001), 或者重复次数超过限制次数(默认是100次), 结束抽样
    2.5 最后的抽样深度为最接近 ASV/OTU 达到一定百分比时的抽样深度。
3. 剔除较小的异常值的数据量后(mean ± 3std), 取最小的样本的数据量与上述中获取的抽样深度比较, 两者取最大值
"""

import os, sys
import numpy as np
import pandas as pd
from itertools import chain
from collections import Counter
import math

#设置抽样的随机数种子, 能够保证每次分析结果完全一致
#如果不设置此值, 结果可能会有细微的误差
np.random.seed(123)

def sampling_depth(otu_reads_data, depth=1, iters=10):    
    #构建抽样完整的数组
    new_otu_reads_data = np.array(list(chain(*[
        otu_reads*[otu_idx] for otu_idx, otu_reads in enumerate(otu_reads_data) if otu_reads!=0]
        )))

    #降序操作, 能够减少内存消耗(numpy.sort无降序参数, 取反排序在取绝对值实现)
    new_otu_reads_data.sort()
    new_otu_reads_data = abs(np.sort(-new_otu_reads_data))
    
    #不放回随机抽样
    #从new_otu_reads_data数组中抽取特定数目depth的元素, 重复iters=10次, 取平均值为抽样结果
    observerd_otu_counts = []
    for i in range(iters):
        observerd_otu = len(Counter(np.random.choice(new_otu_reads_data, size=depth)))
        #observerd_otu = len(Counter(np.random.choice(newlist, size=depth, replace=True))) #replace=True 表示有放回抽样
        observerd_otu_counts.append(observerd_otu)

    #返还平均值
    return np.array(observerd_otu_counts).mean()

import click
@click.command()
@click.option("--feature_table", "-f", required=True, type=click.Path(exists=True),
    help="输出的ASV/OTu特征表结果")
@click.option("--percentage", "-p", default=0.9, type=click.FLOAT,
    help="抽样预期达到的ASV/OTU占总ASV/OTU的百分比,取值范围(0, 1), default, 0.9")
@click.option("--diff_percentage_cutoff", "-d", default=0.0001, type=click.FLOAT,
    help="与抽样预期ASV/OTU百分比误差值, 在此误差值范围内结束抽样, 推荐高选取比percentage高一个级别的精度收敛速度快(但是误差偏大), default:0.0001")
@click.option("--iters", "-i", default=10, type=click.INT,
    help="每次抽样的迭代次数, default: 10")
@click.option("--max_times", "-m", default=100, type=click.INT,
    help="二分法挑选抽样深度最大迭代次数，default: 100")
@click.option("--output", "-o", required=True, type=str,
    help="输出低于预期抽平/采样深度的样本信息到次文件")
def sampling_depth_table(
    feature_table,
    percentage,
    diff_percentage_cutoff,
    iters,
    max_times,
    output,
    ):
    #读取特征表数据
    df = pd.read_table(feature_table, keep_default_na=None, header=1)
    df = df[list(df.columns)[1:]]

    sample_reads = df.sum(axis=0)

    #最大最小数据量(样本总reads)异常值
    max_outlier = sample_reads.mean() + 3*sample_reads.std()
    min_outlier = sample_reads.mean() - 3*sample_reads.std()

    #样本特征值ASVs计数
    feature_sum = (df>0).sum(axis=0)

    #从最小值和最小异常值中取大值作为抽平筛选的取样深度
    valid_min_depth = int(max(sample_reads.min(), min_outlier))

    breakpoint_depths = {}
    for sampleid in df.columns:
        if sampleid.startswith("#"):continue #"OTU ID"列是OTU名称, 不少数值, 不进行计数
        #if sampleid != "L1S105":continue

        #获取每个样本的OTU的reads计数
        #如果是OTU丰度表, 需要把所有丰度值乘以一个系数(例如:1000)取整, 才能进行下列随机抽样
        otu_reads_data = df[sampleid].astype(int)
    
        #样本总reads数目
        total_reads = otu_reads_data.sum()
        #样本总ASV/OTU数目
        total_otu   = len(otu_reads_data[otu_reads_data>0])

        #随机抽样深度最大最小值范围, 后续抽样回逐步缩小次范围, 使得抽样后的otu数目达到特定百分比
        max_target_depth = total_reads
        min_target_depth = 0
    
        #开始的随机抽样深度. 用二分法逼近, math.ceil向上取整
        target_depth = math.ceil(0.5 * max_target_depth)
    
        #进行无放回随机抽样的OTU计数
        target_observed = sampling_depth(otu_reads_data, depth=target_depth, iters=iters)
    
        #抽样的OTU百分比与目标百分比的差值
        diff_percentage = target_observed/total_otu - percentage

        try_times = 0
        while abs(diff_percentage) > diff_percentage_cutoff: 
            if diff_percentage > 0:
                #当OTU计数高于所需要的百分比深度, 用二分法减少抽样深度, 并且把最大抽样范围定为此值
                max_target_depth = target_depth
                target_depth = target_depth - math.ceil(0.5 * (target_depth - min_target_depth))
            else:
                #当OTU计数低于所需要的百分比深度, 用二分法增加抽样深度, 并且把最低抽样范围定为此值
                min_target_depth = target_depth
                target_depth = target_depth + math.ceil(0.5 * (max_target_depth - target_depth))        
            target_observed = sampling_depth(otu_reads_data, depth=target_depth, iters=iters)        
            diff_percentage = target_observed/total_otu - percentage
            try_times += 1
            if try_times > max_times:
                break

        #print(sampleid, target_depth, total_reads, total_otu, try_times, diff_percentage)
        #print(sampleid, target_depth)
        breakpoint_depths[sampleid] = target_depth
    
    #推荐抽平深度
    p_sampling_depth = max(max(breakpoint_depths.values()), valid_min_depth)    

    #输出低于抽平深度样本
    deldata = sample_reads[sample_reads<p_sampling_depth]

    outdata = {
        "min_depth":    int(sample_reads.min()), #必须是整数, 后续可能要用到, 浮点数会报错
        "max_depth":    int(sample_reads.max()), #须是整数, 后续可能要用到
        "mean_depth":   sample_reads.mean(),
        "median_depth": int(sample_reads.median()),  #中位数
        "std_depth":    sample_reads.std(),     #方差
        "p_sampling_depth": p_sampling_depth,
        "filtered_samples_depth(auto)": dict(map(lambda x,y:(x, int(y)), deldata.index, deldata.values))
        }
    import json
    with open(output, 'w') as w:
        json.dump(outdata, w, indent=4)

if __name__ == "__main__":
    sampling_depth_table()
