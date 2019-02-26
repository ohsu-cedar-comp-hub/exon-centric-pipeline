import os


def getSummarySingle(summary_f):
    """Generates a dictionary containing summary information about a non-comparison summary MISO file
    Args:
        summary_f (str/path): path to MISO file generated from MISO -summarize flag
    Returns:
       Dictionary where the MISO event is the key and a single dictionary as the value containing the event's mean, low, high, isoforms, counts, assigned counts, chrom, strand, starts, ends values
    """
    eventToInfo = {}
    for line in open(summary_f):
        if not line.startswith("event_name"):
            try:
                vals = line.strip().split("\t")
                event = vals[0]
                s1m = vals[1]
                s1l = vals[2]
                s1h = vals[3]
                isoforms = vals[4]
                s1c = vals[5]
                s1ac = vals[6]
                s1ac = sum([int(x.split(":")[1]) for x in s1ac.split(",")])
                chrom = vals[7]
                strand = vals[8]
                starts = vals[9]
                ends = vals[9]

                eventToInfo[event] = {}
                eventToInfo[event]['mean'] = list(map(float, s1m.split(",")))
                eventToInfo[event]['low'] = list(map(float, s1l.split(",")))
                eventToInfo[event]['high'] = list(map(float, s1h.split(",")))
                eventToInfo[event]['isoforms'] = isoforms
                eventToInfo[event]['counts'] = s1c
                eventToInfo[event]['assigned counts'] = s1ac
                eventToInfo[event]['chrom'] = chrom
                eventToInfo[event]['strand'] = strand
                eventToInfo[event]['starts'] = starts
                eventToInfo[event]['ends'] = ends
            except:
                pass
        else:
            header = line[1:].strip().split("\t")

    return eventToInfo, header


def getSummary(summary_f):
    """Generates a dictionary containing summary information about a comparison summary MISO file
    Args:
        summary_f (str/path): path to MISO comparison summary file generated from MISO -summarize flag
    Returns:
       Dictionary where the MISO event is the key and a single dictionary as the value containing the event's mean, low, high,delta psi, bayes factor, isoforms, counts, assigned counts, chrom, strand, starts, ends
       values for both sample one and two
    """
    eventToInfo = {}
    for line in open(summary_f):
        if not line.startswith("#"):
            try:
                vals = line.strip().split("\t")
                event = vals[0]
                s1m = vals[1]
                s1l = vals[2]
                s1h = vals[3]
                s2m = vals[4]
                s2l = vals[5]
                s2h = vals[6]
                dpsi = vals[7]
                bf = vals[8]
                isoforms = vals[9]
                s1c = vals[10]
                s1ac = vals[11]
                s1ac = sum([int(x.split(":")[1]) for x in s1ac.split(",")])
                s2c = vals[12]
                s2ac = vals[13]
                s2ac = sum([int(x.split(":")[1]) for x in s2ac.split(",")])
                maxbf = vals[18]
                gene = vals[19]
                symb = vals[20]
                desc = vals[21]

                eventToInfo[event] = {}
                eventToInfo[event]['sample 1 mean'] = list(map(float, s1m.split(",")))
                eventToInfo[event]['sample 1 low'] = list(map(float, s1l.split(",")))
                eventToInfo[event]['sample 1 high'] = list(map(float, s1h.split(",")))
                eventToInfo[event]['sample 2 mean'] = list(map(float, s2m.split(",")))
                eventToInfo[event]['sample 2 low'] = list(map(float, s2l.split(",")))
                eventToInfo[event]['sample 2 high'] = list(map(float, s2h.split(",")))
                eventToInfo[event]['delta psi'] = list(map(float, dpsi.split(",")))
                eventToInfo[event]['bayes factor'] = list(map(float, bf.split(",")))
                eventToInfo[event]['isoforms'] = isoforms
                eventToInfo[event]['sample 1 counts'] = s1c
                eventToInfo[event]['sample 1 assigned counts'] = s1ac
                eventToInfo[event]['sample 2 counts'] = s2c
                eventToInfo[event]['sample 2 assigned counts'] = s2ac
                eventToInfo[event]['max bf'] = float(maxbf)
                eventToInfo[event]['gene'] = gene
                eventToInfo[event]['symb'] = symb
                eventToInfo[event]['desc'] = desc
            except:
                print(event)
        else:
            header = line[1:].strip().split("\t")

    return eventToInfo, header


groupToSamples = {}
samples = []

summarydir = 'samples/miso/summary_vs'
out_f = snakemake.output[0]

samples = snakemake.params.samples
print(samples)
print(len(samples), 'samples')
# Iterate through all samples and save psi values
eventMaster = {}  # event -> sample -> psi values and bfs

for s in samples:
    print(s)
    f = os.path.join(summarydir, s + ".miso_summary")
    eventToInfo, header = getSummarySingle(f)
    for e in eventToInfo:
        if e not in eventMaster:
            eventMaster[e] = {}
        eventMaster[e][s] = [eventToInfo[e]['low'][0], eventToInfo[e]['mean'][0], eventToInfo[e]['high'][0]]
# Iterate through all comparisons and save bfs
for i in range(len(samples)):
    for j in range(i + 1, len(samples)):
        f1 = os.path.join(summarydir, samples[i] + "_vs_" + samples[j] + ".miso_bf")
        f2 = os.path.join(summarydir, samples[j] + "_vs_" + samples[i] + ".miso_bf")
        if os.path.exists(f1):
            f = f1
        else:
            f = f2
        print(f)
        eventToInfo, header = getSummary(f)
        for e in eventToInfo:
            eventMaster[e][samples[i] + "_vs_" + samples[j]] = eventToInfo[e]['max bf']
            if 'gene' not in eventMaster[e]:
                eventMaster[e]['gene'] = eventToInfo[e]['gene']
                eventMaster[e]['symb'] = eventToInfo[e]['symb']
                eventMaster[e]['desc'] = eventToInfo[e]['desc']

out = open(out_f, 'w')
out.write("#Event\t")

for i in range(len(samples)):
    out.write("%s\t%s\t%s\t" % (samples[i] + "_low", samples[i] + "_mean", samples[i] + "_high"))

for i in range(len(samples)):
    for j in range(i + 1, len(samples)):
        out.write(samples[i] + "_vs_" + samples[j] + "\t")

out.write("gene\tsymb\tdesc\n")

for e in eventMaster:
    out.write(e + "\t")
    for i in range(len(samples)):
        if samples[i] in eventMaster[e]:
            out.write("\t".join(list(map(str, eventMaster[e][samples[i]]))) + "\t")
        else:
            out.write("n/a\tn/a\tn/a\t")
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            compname = samples[i] + "_vs_" + samples[j]
            if compname in eventMaster[e]:
                out.write(str(eventMaster[e][compname]) + "\t")
            else:
                out.write("1\t")
    try:
        out.write(eventMaster[e]['gene'] + "\t" + eventMaster[e]['symb'] + "\t" + eventMaster[e]['desc'] + "\n")
    except:
        out.write("n/a\tn/a\tn/a\n")
out.close()


