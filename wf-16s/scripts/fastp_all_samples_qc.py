import sys
from pathlib import Path
import json
import pandas as pd
import numpy as np


def fastp_all_samples_qc(files_fastp_json, outfile=None):
    title = ["Sample", "Clean_Reads", "Total_Base",
             "Q20", "Q30", "Q20_Rate", "Q30_Rate",
             "Average_Length", "GC"]
    df = pd.DataFrame(columns=title)
    for js_path in files_fastp_json:
        js_data = json.loads(open(js_path, "r").read())
        sample = Path(js_path).stem
        mean_lengths = np.array(
            [v for k, v in js_data["summary"]["after_filtering"].items() if k.endswith("mean_length")])
        out = [
            sample,
            js_data["summary"]["after_filtering"]["total_reads"],
            js_data["summary"]["after_filtering"]["total_bases"],
            js_data["summary"]["after_filtering"]["q20_bases"],
            js_data["summary"]["after_filtering"]["q30_bases"],
            js_data["summary"]["after_filtering"]["q20_rate"],
            js_data["summary"]["after_filtering"]["q30_rate"],
            mean_lengths.mean(),
            js_data["summary"]["after_filtering"]["gc_content"],
        ]
        df.loc[len(df)] = out

    if outfile:
        df.to_csv(outfile, index=False, sep="\t")
    else:
        print("\t".join(df.columns))
        for i in range(df.shape[0]):
            print("\t".join([str(v) for v in df.iloc[i].values]))
    return df


if __name__ == "__main__":
    fastp_all_samples_qc(snakemake.input, outfile=snakemake.output[0])
