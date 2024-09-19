import sys, os, re

def combine_alpha_diversity_stat(inputfiles):
    stats = {}
    headers = []
    samples = []
    for inputfile in inputfiles:
        with open(inputfile, 'r') as f:
            for ln, line in enumerate(f):
                cols = [l.strip() for l in line.split("\t")]
                if ln == 0:
                    header = cols[1]
                    headers.append(header)
                    stats[header] = {}
                else:
                    sample, value = cols[:2]
                    stats[header][sample] = value
                    if sample not in samples:
                        samples.append(sample)
    print("SampleID\t" + "\t".join(headers))
    for sample in samples:
        values = [stats[header].get(sample, "") for header in headers]
        print(sample + "\t" + "\t".join(values))

if __name__ == "__main__" :
    combine_alpha_diversity_stat(sys.argv[1:])

