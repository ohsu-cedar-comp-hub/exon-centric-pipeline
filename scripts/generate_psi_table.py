import os
import re


def parseCounts(countsfield):
    """Helper function to parse counts field.
    Args:
        countsfield (slice): column of MISO summary file containing countsfield
    Returns:
       Total number of counts for a given event
    """
    cttot = 0
    iterator = re.finditer("\(.+?\)\:[0-9]+", countsfield)
    for item in iterator:
        isoform, ct = countsfield[item.start(): item.end()].split(":")
        cttot += int(ct)
    return cttot


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


def psiTable(summarydir, out_f, includelist_f=False, minct=1):
    """ Create a table of psi values, where rows are events and columns are samples. You can specify the samples to include in the includelist file.
    Args:
        summarydir (str/path): MISO summary directory
        out_f (str/path): Label for the resulting text file
        includelist_f (bool): An optional file where samples are listed one per line and will dictate which samples are included in the final table
        minct (int): The minimum number of reads mapped to an exon junction to be included in the table
    Returns:
       Nothing. Generates a table of order psi values, where raw read counts for a given event are greater than minct for all samples
    """
    minct = int(minct)
    comps = [f for f in os.listdir(summarydir) if f.endswith("miso_bf")]
    eventToSample = {}
    eventToGene = {}
    samples = {}
    for comp in comps:
        print(comp)
        eventToInfo, header = getSummary(os.path.join(summarydir, comp))
        for event in eventToInfo:
            if event not in eventToSample:
                eventToSample[event] = {}

            sample1 = header[1]
            sample2 = header[4]
            samples[sample1] = 0
            samples[sample2] = 0

            bf = eventToInfo[event]['bayes factor']
            s1c = eventToInfo[event]['sample 1 counts']
            s2c = eventToInfo[event]['sample 2 counts']

            ct1 = parseCounts(s1c)
            ct2 = parseCounts(s2c)

            if min([ct1, ct2]) >= minct:
                if event not in eventToGene:
                    eventToGene[event] = [eventToInfo[event]['gene'],
                                          eventToInfo[event]['symb'], \
                                          eventToInfo[event]['desc']]
                if sample1 not in eventToSample[event]:
                    eventToSample[event][sample1] = \
                        eventToInfo[event]['sample 1 mean'][0]
                if sample2 not in eventToSample[event]:
                    eventToSample[event][sample2] = \
                        eventToInfo[event]['sample 2 mean'][0]

    if includelist_f is not False:
        samples = []
        for line in open(includelist_f):
            samples.append(line.strip())
    else:
        samples = list(samples.keys())
        samples.sort()

    out = open(out_f, 'w')
    print(len(samples), 'samples')
    out.write("#Event\t" + "\t".join(samples) + "\n")
    for event in eventToSample:
        useme = True
        for s in samples:
            if s not in eventToSample[event]:
                useme = False
                break
        if useme:
            out.write(event + "\t")
            out.write("\t".join(map(str, [eventToSample[event][s] \
                                          for s in samples])) + "\t")
            out.write("\t".join(eventToGene[event]) + "\n")
    out.close()


summary_dir = 'samples/miso/summary_vs'
out_f = snakemake.output[0]
psiTable(summary_dir, out_f)


