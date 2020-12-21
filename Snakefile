__author__ = "Joey Estabrook"
__email__ = "estabroj@ohsu.edu"
__license__ = "MIT"

"""Computation Hub exocentric splice site processing pipeline"""


import datetime
import sys
import os
import pandas as pd
import json
from pylab import *
from itertools import combinations,product

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"omic_config.yaml"
project_id = config["project_id"]
with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())
rule_dirs.pop(rule_dirs.index('__default__'))

SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fastq.gz")
controls = ['20-00077', '20-00309', '20-00452', '20-00488', '20-00515', '20-00519', '20-00516', '20-00518', '20-00517', '20-00530', '20-00521', '20-00522', '20-00524', '20-00534', '20-00535', '20-00526', '20-00531', '20-00532', '20-00520', '20-00527', '20-00533', '20-00525', '20-00529', '20-00523', '20-00528', '20-00061']
asxl_mutants = ['13-00281', '13-00522', '13-00659', '14-00021', '14-00240', '14-00504', '14-00608', '14-00670', '14-00832', '15-00057', '15-00073', '15-00123', '15-00171', '15-00395', '15-00383', '15-00479', '15-00572', '15-00578', '15-00777', '15-00870', '15-00976', '15-00990', '16-00145']

#comb1,comb2 = zip(*list(combinations(SAMPLES[:10],2)))
comb1,comb2 = zip(*list(product(controls,asxl_mutants)))

for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)

result_dirs = ['diffexp','tables']
for rule in result_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'results',rule)):
        log_out = os.path.join(os.getcwd(), 'results', rule)
        os.makedirs(log_out)
        print(log_out)


def message(mes):
    sys.stderr.write("|--- " + mes + "\n")


def filter_monotonic(m_fh):
    monotonic_file = pd.read_csv(m_fh,sep='\t',index_col=0)
    mono_regulated = monotonic_file[(abs(monotonic_file['delta_psi']) >= .20 )].index.tolist()

    return mono_regulated


def get_contrast(wildcards):
    """Return each contrast provided in the configuration file"""
    return config["diffexp"]["contrasts"][wildcards.contrast]

def get_mean(wildcards):
    """Return mean value of insert length distribution"""
    fh = "samples/miso/analysis/{}/insert-dist/Aligned.sortedByCoord.out.bam.insert_len".format(wildcards.sample)
    try:
        with open(fh) as f:
            mean,sdev,dispersion,num_pairs = f.readline().strip()[1:].split(',')
        mean = float(mean.split('=')[1])
    except:
        mean = 158.7

    return mean

def get_sdev(wildcards):
    """Return sdev value of insert length distribution"""
    fh = "samples/miso/analysis/{}/insert-dist/Aligned.sortedByCoord.out.bam.insert_len".format(wildcards.sample)
    try:
        with open(fh) as f:
            mean,sdev,dispersion,num_pairs = f.readline().strip()[1:].split(',')
        sdev = float(sdev.split('=')[1])
    except:
        sdev = 43.6
    return sdev


for sample in SAMPLES:
    message("Sample " + sample + " will be processed")


rule all:
    input:
        expand("results/tables/{project_id}_STAR_mapping_statistics.txt", project_id = config['project_id']),
        "data/{project_id}_counts.txt".format(project_id=config['project_id']),
        expand("samples/star_notrim/{sample}_bam/Aligned.sortedByCoord.out.bam.bai",sample=SAMPLES),
        "samples/miso/settings/{type}_map.pkl".format(type=config['event_type']),
        expand("samples/miso/analysis/{sample}/summary/{sample}.miso_summary",sample=SAMPLES),
        expand("samples/miso/comparisons/{sample1}_vs_{sample2}/bayes-factors/{sample1}_vs_{sample2}.miso_bf",sample1=comb1,sample2=comb2),
        expand("samples/miso/sashimi_plots/{contrast}/all_plots_all_launched.txt",contrast = config["diffexp"]["contrasts"]),
        expand("samples/miso/bams/{sample}.bai",sample=SAMPLES),
        "samples/miso/results/{project_id}_psi_table.txt".format(project_id=config['project_id'])

include: "rules/align_rmdp.smk"
include: "rules/miso.smk"
