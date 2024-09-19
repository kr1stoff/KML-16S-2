# -*- coding: utf-8 -*
import sys, os, re
import pandas as pd

"""
用于qiime2频数特征表中生成krona的输入格式文件
    all表示输出所有样本
"""

def table_to_krona(inputfile, target_samples=["all"]): #L1S105,L1S140
    df = pd.read_csv(inputfile, keep_default_na=False, header=1, sep="\t")    
    headers = list(df.columns)
    samplelist = [h for h in headers if h not in ["#OTU ID"]]
    if target_samples == ["all"]:
        target_samples = samplelist
    else:
        target_samples = list(set(target_samples) & set(samplelist))

    sps_df = df[["#OTU ID"] + target_samples]
    for row in sps_df.iterrows():
        index = row[0]
        data = row[1]
        taxs = [re.sub("^[kdpcofgs]__", "", tag) for tag in data["#OTU ID"].split(";")]
        sp_counts = list(data[target_samples])
        #筛除所有结果都是0的物种列
        if sum(sp_counts) == 0:
            continue            

        #过滤部分注释为___, Unassigned等的行
        levels = [("unknown" if tax=="__" else tax) for tax in taxs]
        if levels.count("unknown") >= 2: continue
        if levels.count("Unassigned") >= 1: continue

        out = str(sum(sp_counts)) + "\t" + "\t".join(levels)
        print(out)

if __name__ == "__main__":
    table_to_krona(sys.argv[1], target_samples=sys.argv[2:])
