import dill as pickle
import operator
import os

ensGeneMap_f = snakemake.input.map_file

with open(ensGeneMap_f, 'rb') as fp:
    eventToInfo = pickle.load(fp)

dirs = snakemake.input.sample_list
outdir = 'samples/miso/summary_vs'

for dir in dirs:
    print(dir)
    fname = os.path.basename(dir)
    sample1, sample2 = fname.split("_vs_")
    bf_f = os.path.join(os.getcwd(), dir, 'bayes-factors', fname + ".miso_bf")
    data = []
    for line in open(bf_f):
        if not line.startswith("event_name"):
            vals = line.strip().split("\t")
            bf = vals[8]
            if ',' in bf:
                maxbf = max(map(float, bf.split(",")))
            else:
                maxbf = float(bf)
            vals.append(maxbf)
            data.append(vals)
    data.sort(key = operator.itemgetter(-1), reverse = True)
    out = open(os.path.join(outdir, fname + ".miso_bf"), 'w')
    out.write("\t".join(["#event_name", sample1 + "_posterior_mean", sample1 + "_ci_low", sample1 + "_ci_high", \
                         sample2 + "_posterior_mean", sample2 + "_ci_low", sample2 + "_ci_high", "diff",
                         "bayes_factor", "isoforms", sample1 + "_counts", sample1 + "_assigned_counts", \
                         sample2 + "_counts", sample2 + "_assigned_counts", "chrom", "strand", "mRNA_starts",
                         "mRNA_ends", "max_bf"]) + "\n")
    for item in data:
        event = item[0]
        try:
            item.extend(eventToInfo[item[0]])
        except:
            item.extend(["n/a", "n/a", "n/a"])
        out.write("\t".join(map(str, item)) + "\n")
    out.close()
