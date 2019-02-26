from pylab import *
import os
import operator
import sys


def monotonic(consolidated_f, groups_f, minbf, nshuffles, out_f):
    """ Find events that change monotonically (significantly).
    Args:
        consolidated_f (str/path): Consolidated summary file
        groups_f (str/path): Directory containing <groups_f>.txt. File denotes sample and specified group, respectfully, one sample per line. Tab delimited.
        minbf (int): Denotes the minimum bayes factor metric to filter events by.
        nshuffles (int): Denotes the number of randomized shuffles used to generate a Z-score for a given significant event.
        out_f (str): Name of Consolidated summary file to be generated
    Returns:
        Nothing. Generates a text file containing the consolidated Monotonic summary information
    """
    minbf = float(minbf)
    nshuffles = int(nshuffles)
    groupToSamples = {}
    grouporder = []
    samples = []
    for line in open(groups_f):
        sample, g = line.strip().split()
        samples.append(sample)
        if g not in grouporder:
            grouporder.append(g)
        try:
            groupToSamples[g].append(sample)
        except:
            groupToSamples[g] = [sample]
    groups = groupToSamples.keys()
    ctllabel = grouporder[0]
    explabel = grouporder[-1]
    # Identify between-group blocks and get their indices for reference
    groupidx = []
    lowidx = 0
    for g in grouporder:
        n = len(groupToSamples[g])
        groupidx.append([lowidx, lowidx + n])
        lowidx += n
    compidx = []
    for i in range(len(groups)):
        idx1 = groupidx[i]
        for j in range(i + 1, len(groups)):
            # Analyze each block of groups
            idx2 = groupidx[j]
            for k in range(idx1[0], idx1[1]):
                for l in range(idx2[0], idx2[1]):
                    compidx.append([k, l])
    compidx = list(zip(*compidx))
    data = []
    print("Reading data...")
    counter = 0
    for line in open(consolidated_f):
        vals = line.strip().split("\t")
        if line.startswith("#"):
            header = vals
        else:
            sampleToBF = {}
            event = vals[0]
            symb = vals[-2]
            sampleToPsi = {}
            skipme = False
            for i in range(len(header)):
                if "_mean" in header[i]:
                    sample = header[i][:-5]
                    if sample in samples:
                        if vals[i] == 'n/a':
                            skipme = True
                            break
                        else:
                            sampleToPsi[sample] = float(vals[i])
                elif "_vs_" in header[i]:
                    s1, s2 = header[i].split("_vs_")
                    bf = float(vals[i])
                    if s1 not in sampleToBF:
                        sampleToBF[s1] = {}
                    if s2 not in sampleToBF:
                        sampleToBF[s2] = {}
                    sampleToBF[s1][s2] = bf
                    sampleToBF[s2][s1] = bf
            if not skipme:
                mat = zeros((len(samples), len(samples)), dtype = 'int')
                for i in range(len(samples)):
                    for j in range(i + 1, len(samples)):
                        try:
                            psi1 = sampleToPsi[samples[i]]
                            psi2 = sampleToPsi[samples[j]]
                            bf = sampleToBF[samples[i]][samples[j]]
                            dpsi = psi2 - psi1
                            if bf >= minbf:
                                mat[i, j] = sign(dpsi)
                        except:
                            pass
                shuffledsamples = list(samples)
                mat = zeros((len(samples), len(samples), nshuffles), dtype = 'int')
                for n in range(nshuffles + 1):
                    for i in range(len(samples)):
                        for j in range(i + 1, len(samples)):
                            try:
                                psi1 = sampleToPsi[shuffledsamples[i]]
                                psi2 = sampleToPsi[shuffledsamples[j]]
                                bf = sampleToBF[shuffledsamples[i]][shuffledsamples[j]]
                                dpsi = psi2 - psi1
                                if bf >= minbf:
                                    mat[i, j, n] = sign(dpsi)
                            except:
                                pass
                    shuffle(shuffledsamples)
                # Get the mean value in exp vs. ctl
                ctlpsi = array([sampleToPsi[x] for x in sampleToPsi if x in groupToSamples[ctllabel]])
                exppsi = array([sampleToPsi[x] for x in sampleToPsi if x in groupToSamples[explabel]])
                signs = mat[compidx[0], compidx[1], :]
                signs = signs.sum(axis = 0)
                m = signs[1:].mean()
                stdev = signs[1:].std()
                z = (signs[0] - m) / stdev
                if signs[0] == 0:
                    z = 0
                data.append([event, symb, signs[0], m, stdev, round(ctlpsi.mean(), 2), round(exppsi.mean(), 2),round(exppsi.mean() - ctlpsi.mean(), 2), z])
                counter += 1
    data.sort(key = operator.itemgetter(-1), reverse = True)
    out = open(out_f, 'w')
    out.write("#Event\tSymb\tTrueval\tMean\tStd\t%s_psi\t%s_psi\tdelta_psi\tZ-score\n" % (ctllabel, explabel))
    for item in data:
        out.write("\t".join(map(str, item)) + "\n")
    out.close()

consolidated_f = snakemake.input.summary
groups_f = snakemake.input.group_f 
minbf = snakemake.params.minbayes
nshuffles = snakemake.params.shuffles
out_f = snakemake.output[0]
monotonic(consolidated_f, groups_f, minbf, nshuffles, out_f)
