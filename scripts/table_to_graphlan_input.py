# -*- coding: utf-8 -*
import sys, os, re
import pandas as pd
import math

"""
用于从qiime2丰度表中生成graphlan的输入格式文件
参考:https://github.com/biobakery/biobakery/wiki/metaphlan3
"""

def table_to_graphlan(inputfile):
    df = pd.read_csv(inputfile, keep_default_na=False, header=1, sep="\t")    
    headers = list(df.columns)
    samplelist = [h for h in headers if h not in ["#OTU ID", "taxonomy"]]
    trim_headers =["Name"] + samplelist
    print("\t".join(trim_headers))

    valid_row_counts = 0
    uniq_taxs = []
    for row in df.iterrows():
        index = row[0]
        data = row[1]
        taxs = [tag for tag in data["#OTU ID"].split(";")]
        taxs = [re.sub("\.", "_", tax) for tax in taxs]

        freqs = [str(data[h]) for h in headers if h not in ["#OTU ID", "taxonomy"]]
        end_level_name = (taxs[-1] if taxs[-1] != "__" else "unknown")
        if end_level_name == "unknown":continue

        flt_tags = ["uncultured", "metagenome", "unidentified", "human_gut"]
        if any(re.search(flt_tag, end_level_name, re.I) for flt_tag in flt_tags):
            continue
        
        new_taxs = ".".join(taxs[:-1]) + "." + end_level_name
        if new_taxs in uniq_taxs:
            continue
        uniq_taxs.append(end_level_name)

        valid_row_counts += 1
        print("%s\t%s"%(new_taxs, "\t".join(freqs)))
    
if __name__ == "__main__":
    table_to_graphlan(sys.argv[1])
