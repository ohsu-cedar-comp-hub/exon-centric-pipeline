import os
import textwrap
import pandas as pd

map_stats = pd.read_csv(snakemake.input[0],sep='\t',index_col=0)
map_stats.columns = [c.split('_')[0] for c in map_stats.columns]
meta = pd.read_csv(snakemake.params[0],sep='\t',index_col=0)
meta = meta.reindex(map_stats.columns)
contrast = snakemake.params.contrast
condition = snakemake.params.condition[0]
out_file = snakemake.output[0]

print(contrast)
print(condition)
print(out_file)

map_colors = ["#CC0011","#FF8800"]

sub_meta = meta.loc[meta[condition].isin(contrast)].sort_values(condition)
cmap = dict(zip(contrast,map_colors))
sub_meta['color'] = sub_meta.Condition.map(cmap)

bam_directory = os.path.join(os.getcwd(),'samples','miso','bams')
miso_prefix = os.path.join(os.getcwd(),'samples','miso','analysis')

sashimi_settings = """
[data]
# directory where BAM files are
bam_prefix = {bam_directory}
# directory where MISO output is
miso_prefix = {miso_prefix}

bam_files = {bam_files}

miso_files = {miso_files}

[plotting]
# Dimensions of figure to be plotted (in inches)
fig_width = 9
fig_height = 7
# Factor to scale down introns and exons by
intron_scale = 30
exon_scale = 4
# Whether to use a log scale or not when plotting
logged = False
font_size = 6

bar_posteriors = False

# Axis tick marks
nyticks = 3
nxticks = 4

# Whether to show axis labels
show_ylabel = False
show_xlabel = True

# Whether to plot posterior distributions inferred by MISO
show_posteriors = True

# Whether to plot the number of reads in each junction
number_junctions = True

resolution = .5
posterior_bins = 40
gene_posterior_ratio = 5

# List of colors for read densities of each sample
colors = {colors}

# Number of mapped reads in each sample
# (Used to normalize the read density for RPKM calculation)
coverages = {coverages}

# Bar color for Bayes factor distribution
# plots (--plot-bf-dist)
# Paint them blue
bar_color = "b"

# Bayes factors thresholds to use for --plot-bf-dist
bf_thresholds = [0, 1, 2, 5, 10, 20]
"""


out_f = open(out_file,'w')

reformatted_cmd = textwrap.dedent(sashimi_settings).strip()
context = {'bam_directory':bam_directory,'miso_prefix':miso_prefix, 'bam_files':sub_meta.index.tolist(),'miso_files':sub_meta.index.tolist(),'colors':sub_meta['color'].values.tolist(),'coverages':list(map(int,map_stats.loc['                   Uniquely mapped reads number |', sub_meta.index].values.tolist()))}
out_f.write(reformatted_cmd.format(**context))
out_f.close()

