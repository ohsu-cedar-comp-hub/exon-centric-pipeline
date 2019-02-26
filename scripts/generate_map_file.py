import pandas as pd
import dill as pickle



events = snakemake.params[0]
relnm = snakemake.output[0]

gff = pd.read_csv(events, sep='\t', skiprows=1, header=None)
gff_filt = gff[(gff[2].isin(['gene']))]
event_type = gff_filt[1]
gff_expand = gff_filt[8].str.split(';',expand=True)
event = gff_expand[0].str.split('=', expand=True)[1]
ens = gff_expand[3].str.split('=', expand=True)[1]
gene_symb = gff_expand[4].str.split('=', expand=True)[1]

joined = pd.concat([event, gene_symb, ens, event_type],axis=1)
joined.columns=['Event', 'Gene', 'Ensembl', 'Type']
joined_dict=joined.set_index('Event').T.to_dict('list')


with open(relnm, 'wb') as f:
    pickle.dump(joined_dict, f, protocol = -1)
