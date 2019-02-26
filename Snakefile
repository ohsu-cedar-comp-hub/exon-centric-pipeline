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


timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"omic_config.yaml"
project_id = config["project_id"]
with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())
rule_dirs.pop(rule_dirs.index('__default__'))

SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fastq.gz")

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


for sample in SAMPLES:
    message("Sample " + sample + " will be processed")


rule all:
    input:
        expand("results/tables/{project_id}_STAR_mapping_statistics.txt", project_id = config['project_id']),
        "data/{project_id}_counts.txt".format(project_id=config['project_id']),
        expand("samples/star_notrim/{sample}_bam/Aligned.sortedByCoord.out.bam.bai",sample=SAMPLES),
        "samples/miso/settings/{type}_map.pkl".format(type=config['event_type']),
        expand("samples/miso/sashimi_plots/{contrast}/all_plots_all_launched.txt",contrast = config["diffexp"]["contrasts"]),
        expand("samples/miso/bams/{sample}.bai",sample=SAMPLES),
        "samples/miso/results/{project_id}_psi_table.txt".format(project_id=config['project_id'])

include: "rules/align_rmdp.smk"
include: "rules/miso.smk"
