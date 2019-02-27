import os
import subprocess

input = snakemake.input[0]
path = snakemake.output[0]

input = os.path.join('/'.join(input.split('/')[:-2]))
out = os.path.join('/'.join(path.split('/')[:-2]))

bash_cmd ="summarize_miso --summarize-samples {input} {out}".format(input=input,out=out)

process = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
out_shell, error = process.communicate()
