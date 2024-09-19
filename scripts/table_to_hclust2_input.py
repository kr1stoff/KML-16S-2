# -*- coding: utf-8 -*
import sys, os, re
import pandas as pd
import math
"""
用于从qiime2丰度表中生成hclust2的输入格式文件
要求:
    去除无意义的菌名(__, unclutrured, metagenome, unidentifie, unknown)
    菌名列去重
    并且更加样本数目, 物种数目, hclust2下一步参数ftop, 自动生成下一步比较合适的hclust2参数值 flabel_size, slabel_size, cell_aspect_ratio
"""

def table_to_hclust2(inputfile, ftop=100):
    df = pd.read_csv(inputfile, keep_default_na=False, header=1, sep="\t")    
    headers = list(df.columns)
    "表头可能有taxnomy(lefse需要)非样信息, 需要去除"
    samplelist = [h for h in headers if h not in ["#OTU ID", "taxonomy"]]
    trim_headers =["Name"] + samplelist
    print("\t".join(trim_headers))

    valid_row_counts = 0
    uniq_taxid = []
    for row in df.iterrows():
        index = row[0]
        data = row[1]
        #得到物种各层级名称
        taxs = [re.sub("^[kdpcofgs]__", "", tag) for tag in data["#OTU ID"].split(";")]
        freqs = [str(data[h]) for h in headers if h not in ["#OTU ID", "taxonomy"]]
        #最后一个层级
        end_level_name = (taxs[-1] if taxs[-1] != "__" else "unknown")
        if end_level_name == "unknown":continue

        #去除无意义的菌名
        flt_tags = ["uncultured", "metagenome", "unidentified", "human_gut"]
        if any(re.search(flt_tag, end_level_name, re.I) for flt_tag in flt_tags):
            continue
        
        #菌名列去重
        if end_level_name in uniq_taxid:
            continue
        uniq_taxid.append(end_level_name)

        valid_row_counts += 1
        print("%s\t%s"%(end_level_name, "\t".join(freqs)))
    
    #自动生成下一步比较合适的hclust2参数值
    valid_row_counts = min(ftop, valid_row_counts)
    auto_slalel_size =  int(min(10, math.ceil(300.0/len(samplelist))))
    auto_flabel_size =  int(min(10, math.ceil(0.5*300.0/valid_row_counts)))
    auto_cell_ratio =  0.5*len(samplelist)/valid_row_counts
    #sys.stderr.write("%s\t%s\n"%(len(samplelist), valid_row_counts))
    sys.stderr.write("%s\t%s\t%s\n"%(auto_slalel_size, auto_flabel_size, auto_cell_ratio))

if __name__ == "__main__":
    table_to_hclust2(sys.argv[1], ftop=int(sys.argv[2]))
