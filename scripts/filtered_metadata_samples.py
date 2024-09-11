# -*- coding: utf-8 -*
import os, sys, re

"""
抽平后可能删掉某些样本, 因此需要生成删掉这些样本信息metadata表, 方便后续分析
"""
def metadata_samples_filter(metafile, statsfile, depth):
    #筛选低于depth的样本名称
    filter_samples = []
    with open(statsfile, 'r') as f:
        for ln, line in enumerate(f):
            if line.startswith("#"):continue
            cols = [l.strip() for l in line.split("\t")]
            if ln == 0:
                header = cols
                sample_idx = header.index("sample-id")
                final_depth_idx = header.index("non-chimeric")
            else:
                sample = cols[sample_idx]
                final_depth = cols[final_depth_idx]
                if int(final_depth) < depth:
                    filter_samples.append(sample)
    
    with open(metafile, 'r') as f:
        for ln, line in enumerate(f):
            cols = [l.strip() for l in line.split("\t")]
            #原样输出#开头行, 不做任何处理
            if line.startswith("#"):
                print(line.strip())
                continue

            if ln == 0:
                #从表头获取样本所在列index
                print(line.strip())
                header = cols
                sample_idx = 0
                for sp_header in ["sample-id", "sample id"]:
                    try:
                        sample_idx = cols.index(sp_header)
                        continue
                    except:
                        pass
            else:
                if cols[sample_idx] in filter_samples:
                    sys.stderr.write("###filtered line %s: %s\n"%(ln, line.strip()))
                else:
                    print(line.strip())


if __name__ == "__main__":
    metadata_samples_filter(
        sys.argv[1], 
        sys.argv[2],
        int(float(sys.argv[3])))
