"""Minimal panel: how few genes reproduce the module score?
RULE FIXED BEFORE RUNNING. Selection uses the DISCOVERY cohort only:
  candidates = module genes detected in >=90% of discovery samples;
  rank by the within-country correlation of the gene with the full module score (discovery);
  take the top k for k = 5, 10, 15, 20.
Evaluation: (i) agreement with the full 77-gene module score (Spearman) in the discovery cohort and, out of sample,
in the sewage cohort; (ii) within-plant reliability (ICC) in discovery; (iii) whether the two properties the metric is
for survive - correlation with intI1 (discovery) and with national antibiotic consumption (sewage).
Reported whatever the outcome."""
import warnings; warnings.filterwarnings("ignore")
import json, numpy as np, pandas as pd
from scipy.stats import spearmanr
U='/mnt/user-data/uploads/'
g=pd.read_csv('gene_table_final.csv',index_col=0); core=list(g.index)
mm=pd.read_csv('module_membership.csv').set_index('gene').loc[core]
MOD=[x for x in core if bool(mm.loc[x,'in_module'])]
z=pd.read_excel(U+'41467_2025_59019_MOESM8_ESM.xlsx','Fig. S2b and Fig. S4',index_col=0)
G=z[core].astype(float); L=np.log10(G+G.replace(0,np.nan).min()/2)
full=L[MOD].mean(1)
ctry=z.CountryRegion_full
dm=lambda v: v-v.groupby(ctry).transform('mean')
cand=[x for x in MOD if (G[x]>0).mean()>=0.90]
rank=sorted(cand,key=lambda x: -spearmanr(dm(L[x]),dm(full))[0])
def icc1(v,grp):
    d=pd.DataFrame({'v':v,'g':grp}); d=d[d.g.map(d.g.value_counts())>=2]
    k=d.g.value_counts(); n=len(k); N=len(d)
    msb=(k*(d.groupby('g').v.mean()-d.v.mean())**2).sum()/(n-1); msw=((d.v-d.groupby('g').v.transform('mean'))**2).sum()/(N-n)
    kbar=(N-(k**2).sum()/N)/(n-1); return float((msb-msw)/(msb+(kbar-1)*msw))
rho_int=pd.read_csv('intI1_rho.csv',index_col=0).iloc[:,0]
# sewage cohort
m=pd.read_pickle('martiny_matrix.pkl')
meta=pd.read_excel(U+'41467_2025_66070_MOESM4_ESM.xlsx').drop_duplicates('complete_name').set_index('complete_name')
ref=pd.read_csv('bacterial_reference.csv',index_col=0).iloc[:,0]; gp=meta.genepid.to_dict()
keep=[i for i in m.index if i in meta.index and gp.get(i) in ref.index]; m=m.loc[keep]; meta=meta.loc[keep]
bac=np.array([ref[gp[i]] for i in keep],dtype=float)
S3=json.load(open('frozen_spec3.json')); mapped=S3['gene_mapping']['module_mapped']
X=m.values.copy()
for i in range(X.shape[0]):
    r=X[i]; p=r[r>0]
    if len(p): r[r==0]=0.65*p.min()
LD=pd.DataFrame(np.log(X)-np.log(bac)[:,None],index=m.index,columns=m.columns)
prev=(m>0).mean()
MODs=[p for x in MOD for p in mapped.get(x,[]) if p in m.columns and prev[p]>=0.5]
full_s=LD[sorted(set(MODs))].mean(1)
F=pd.read_csv('F_country_table.csv')
def sewage_sub(genes):
    ids=sorted({p for x in genes for p in mapped.get(x,[]) if p in m.columns})
    return (LD[ids].mean(1) if ids else None), len(ids)
OUT={'candidates':len(cand),'ranked_top20':rank[:20],'full_module_genes':len(MOD),'results':{}}
for k in [5,10,15,20]:
    sub=rank[:k]; sc=L[sub].mean(1)
    ss,nids=sewage_sub(sub)
    r_sew=float(spearmanr(ss,full_s)[0]) if ss is not None else None
    cons=None
    if ss is not None:
        cs=ss.groupby(meta.country.astype(str).values).mean(); v=F.country.map(cs); cons=list(map(float,spearmanr(v,F.ddd)[:2]))
    OUT['results'][k]={'genes':sub,'rho_with_full_discovery':float(spearmanr(sc,full)[0]),
        'ICC_discovery':round(icc1(sc,z.WWTPID),3),'rho_with_full_sewage':r_sew,'sewage_ids_matched':nids,
        'mean_rho_intI1':round(float(rho_int[sub].mean()),3),'consumption_rho_P':cons}
OUT['full_reference']={'ICC_discovery':round(icc1(full,z.WWTPID),3),'mean_rho_intI1':round(float(rho_int[MOD].mean()),3),
    'consumption_rho_P':list(map(float,spearmanr(F.country.map(full_s.groupby(meta.country.astype(str).values).mean()),F.ddd)[:2]))}
json.dump(OUT,open('minimal_panel.json','w'),indent=1)
print(json.dumps({k:v for k,v in OUT.items() if k!='ranked_top20'},indent=1)[:2600])
