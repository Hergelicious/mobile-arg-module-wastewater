"""Bacterial reference for the pre-specified ALR normalisation: raw bacterial mOTU fragment counts per sample
(kingdom_motus_pad_agg.csv, Zenodo 14652833), summed over lanes and keyed by genepid. The length-adjusted value is inf
for four records, so the raw count is used."""
import re, pandas as pd
U='/mnt/user-data/uploads/'
lng=pd.read_csv(U+'kingdom_motus_pad_agg.csv'); b=lng[lng.name=='Bacteria'].copy()
meta=pd.read_excel(U+'41467_2025_66070_MOESM4_ESM.xlsx').drop_duplicates('complete_name')
run2gp=dict(zip(meta.ena_run_acc.astype(str),meta.genepid))
def gp(x):
    m=re.search(r'_(\d{7})_',str(x)); return int(m.group(1)) if m else run2gp.get(str(x))
b['genepid']=b.run_accession.map(gp)
b.dropna(subset=['genepid']).groupby('genepid').fragmentCountAln.sum().to_csv('bacterial_reference.csv')
