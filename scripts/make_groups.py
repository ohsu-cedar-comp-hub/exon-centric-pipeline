import os
import pandas as pd


contrast = snakemake.params.contrast
meta_f = snakemake.params.meta
condition = snakemake.params.condition[0]
out_f = snakemake.output[0]

meta = pd.read_csv(meta_f,sep='\t',index_col=0)

sub_meta = meta[(meta[condition].isin(contrast))].sort_values(condition)
filt_meta = sub_meta.loc[:,condition].reset_index()
print(filt_meta)
filt_meta.to_csv(out_f,header=False,index=False,sep='\t')
