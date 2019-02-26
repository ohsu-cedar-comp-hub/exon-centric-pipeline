import os
import subprocess

bam = snakemake.input[0]
insert_dist_f = snakemake.input[1]
events = snakemake.params.events
read_len = snakemake.params.read_len

final_out = snakemake.output[0]
out = os.path.join('/'.join(final_out.split('/')[:-2]))


with open(insert_dist_f) as f:
    mean,sdev,dispersion,num_pairs = f.readline().strip()[1:].split(',')

mean = float(mean.split('=')[1])
sdev = float(sdev.split('=')[1])

bash_cmd = "miso --run {events} {bam} --output-dir {out} --read-len {read_len} --paired-end {mean} {sdev}".format(events=events, bam=bam, out=out, read_len=read_len, mean=mean, sdev=sdev)
print('bass command:',bash_cmd)

print('Calling process')
process = subprocess.call(bash_cmd.split())

print('exiting process')
out_final = open(final_out,'w')
out_final.close()
