#!/usr/bin/env python2.7
import os
import subprocess

bam = snakemake.input[0]
events = snakemake.params.events
read_len = snakemake.params.read_len

final_out = snakemake.output[0]
out = os.path.join('/'.join(final_out.split('/')[:-2]))




bash_cmd = "miso --run {events} {bam} --output-dir {out} --read-len {read_len}".format(events=events, bam=bam, out=out, read_len=read_len)
print('bass command:',bash_cmd)

print('Calling process')
process = subprocess.call(bash_cmd.split())

print('exiting process')
out_final = open(final_out,'w')
out_final.close()
