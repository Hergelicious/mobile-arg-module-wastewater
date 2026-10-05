"""Robustness of the module score to sequencing depth, normalisation and annotation/database choice.
Uses only data already analysed; reports agreement with the score as computed in the manuscript.
Nothing here recomputes or replaces any previously reported result."""
import warnings; warnings.filterwarnings("ignore")
import json, numpy as np, pandas as pd
from scipy.stats import spearmanr
rng=np.random.default_rng(606)
U='/mnt/user-data/uploads/'
S3=json.load(open('frozen_spec3.json')); mapped=S3['gene_mapping']['module_mapped']
g=pd.read_csv('gene_table_final.csv',index_col=0); core=list(g.index)
mm=pd.read_csv('module_membership.csv').set_index('gene').loc[core]
MOD=[x for x in core if bool(mm.loc[x,'in_module'])]
OUT={}
# ---------- sewage cohort: depth, normalisation ----------
m=pd.read_pickle('martiny_matrix.pkl')
meta=pd.read_excel(U+'41467_2025_66070_MOESM4_ESM.xlsx').drop_duplicates('complete_name').set_index('complete_name')
ref=pd.read_csv('bacterial_reference.csv',index_col=0).iloc[:,0]; gp=meta.genepid.to_dict()
keep=[i for i in m.index if i in meta.index and gp.get(i) in ref.index]; m=m.loc[keep]; meta=meta.loc[keep]
bac=np.array([ref[gp[i]] for i in keep],dtype=float)
prev=(m>0).mean(); MODs=sorted({p for x in MOD for p in mapped.get(x,[]) if p in m.columns and prev[p]>=0.5})
def alr(mat):
    X=mat.values.copy().astype(float)
    for i in range(X.shape[0]):
        r=X[i]; p=r[r>0]
        if len(p): r[r==0]=0.65*p.min()
    return pd.DataFrame(np.log(X)-np.log(bac)[:,None],index=mat.index,columns=mat.columns)
full=alr(m)[MODs].mean(1)
# depth: multinomial subsampling of each sample's ARG fragments
tot=m.sum(1).values
dep={}
for frac in [0.10,0.25,0.50,0.75]:
    sub=np.empty_like(m.values)
    for i in range(m.shape[0]):
        p=m.values[i]; s=p.sum()
        n=max(int(s*frac),1)
        sub[i]=rng.multinomial(n,p/s)*(s/n)      # rescale so units stay comparable
    sd=alr(pd.DataFrame(sub,index=m.index,columns=m.columns))[MODs].mean(1)
    genes=(pd.DataFrame(sub,index=m.index,columns=m.columns)>0).sum(1)
    dep[f'{int(frac*100)}%']={'rho_with_full_depth':round(float(spearmanr(sd,full)[0]),3),
                              'median_module_genes_detected':int(np.median((pd.DataFrame(sub,columns=m.columns)[MODs]>0).sum(1)))}
OUT['depth_subsampling_sewage']=dep
# normalisation: ALR vs CLR vs simple fraction of total ARGs
X=m.values.copy().astype(float)
for i in range(X.shape[0]):
    r=X[i]; p=r[r>0]
    if len(p): r[r==0]=0.65*p.min()
clr=pd.DataFrame(np.log(X)-np.log(X).mean(1,keepdims=True),index=m.index,columns=m.columns)[MODs].mean(1)
frac_sc=np.log(m[MODs].sum(1)/m.sum(1))
OUT['normalisation_sewage']={'ALR_vs_CLR_rho':round(float(spearmanr(full,clr)[0]),3),
                             'ALR_vs_fraction_rho':round(float(spearmanr(full,frac_sc)[0]),3),
                             'CLR_vs_fraction_rho':round(float(spearmanr(clr,frac_sc)[0]),3)}
# ---------- discovery cohort: annotation/database choice ----------
z=pd.read_excel(U+'41467_2025_59019_MOESM8_ESM.xlsx','Fig. S2b and Fig. S4',index_col=0)
G=z[core].astype(float); L=np.log10(G+G.replace(0,np.nan).min()/2)
card=L[MOD].mean(1)
cross=[x for x in MOD if mapped.get(x)]                       # module genes that exist in ResFinder/PanRes nomenclature
OUT['annotation_discovery']={'module_genes_CARD':len(MOD),'module_genes_with_cross_database_name':len(cross),
    'rho_CARD_vs_cross_database_subset':round(float(spearmanr(card,L[cross].mean(1))[0]),3)}
# depth proxy in discovery: rarefaction-effort covariate already available
OUT['note']='Depth subsampling is applied to the sewage counts, the only cohort with raw gene counts; the discovery cohort provides per-cell abundances only.'
json.dump(OUT,open('robustness_metric.json','w'),indent=1); print(json.dumps(OUT,indent=1))
