import os
import pandas as pd
import subprocess


def ensure_dir(relnm):
    """ Accept relative filepath string, create it if it doesnt already exist
        return filepath string

    Args:
        relnm (str) : Relative name/path

    Returns:
        relnm (str)

    """

    d = os.path.join(os.getcwd(), relnm)
    if not os.path.exists(d):
        print('--- path does not exist : {} ---'.format(d))
        print('--- constructing path : {} ---'.format(d))
        os.makedirs(d)

    return relnm


def filter_monotonic(m_fh):
    monotonic_file = pd.read_csv(m_fh,sep='\t',index_col=0)
    mono_regulated = monotonic_file[(abs(monotonic_file['delta_psi']) >= .20 )].index.tolist()

    return mono_regulated


monotonic_file = snakemake.input.mono_fh
settings_file = snakemake.input.ps_fh
events = snakemake.params.events

mono_regulated = filter_monotonic(monotonic_file)
output = snakemake.output[0]
output_dir = os.path.join(*output.split('/')[:-1])
ensure_dir(output_dir)

print(output)


if len(mono_regulated) > 0:
    for event in mono_regulated[-2:-1]:
        shell_cmd = "sashimi_plot  --plot-event {event} {events} {settings_file} --output-dir {output_dir}".format(event=event, events=events, settings_file=settings_file, output_dir=output_dir)
        print(shell_cmd)
        process = subprocess.Popen(shell_cmd.split(), stdout=subprocess.PIPE)
        out_shell, error = process.communicate()
else:
    pass

print(output)
out = os.path.join(os.getcwd(),output)
print(out)
out_f = open(out,'w')
out_f.close()
