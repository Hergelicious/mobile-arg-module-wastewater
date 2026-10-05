"""(B) Benchmark of candidate wastewater-resistome surveillance metrics.
Reads only existing outputs; writes benchmark.json. No previously reported result is recomputed or changed.

Metrics (all computed per sample in both cohorts where possible):
  total_ARG            total ARG abundance (per cell in sludge; ALR vs bacteria in sewage)   - current practice
  richness             number of core/detected ARG genes                                      - current practice
  acquired_sum         summed abundance of acquired (mobile-associated / ResFinder) genes     - intermediate
  module_score         mean abundance of the label-free module genes                          - proposed
  module_share         module genes as a fraction of the resistome                            - proposed
Criteria:
  reliability  within-plant ICC(1) across replicate samples (discovery cohort)
  stability    rank correlation of city means between earliest and latest sampling year (sewage cohort)
  response     rank correlation with national antibiotic consumption across countries (sewage cohort)
  country_R2   marginal R2 of country at city level (sewage cohort) - descriptive, not a claim of geography
"""
import warnings; warnings.filterwarnings("ignore")
import json, numpy as np, pandas as pd
from scipy.stats import spearmanr
U='/mnt/user-data/uploads/'
g=pd.read_csv('gene_table_final.csv',index_col=0); core=list(g.index)
mm=pd.read_csv('module_membership.csv').set_index('gene').loc[core]
MODg=[x for x in core if bool(mm.loc[x,'in_module'])]; MOBg=list(g.index[g.M_prior])
z=pd.read_excel(U+'41467_2025_59019_MOESM8_ESM.xlsx','Fig. S2b and Fig. S4',index_col=0)
G=z[core].astype(float); pc=G.replace(0,np.nan).min()/2
L=np.log10(G+pc)
sl=pd.DataFrame({'total':np.log10(z.ARG_per_cell),'richness':(G>0).sum(1),
                 'acquired_sum':np.log10(G[MOBg].sum(1)),'module_score':L[MODg].mean(1),
                 'module_share':np.log10(G[MODg].sum(1)/G[core].sum(1))},index=z.index)
def icc1(v,grp):
    d=pd.DataFrame({'v':v,'g':grp}); d=d[d.g.map(d.g.value_counts())>=2]
    k=d.g.value_counts(); n=len(k); N=len(d)
    msb=(k*(d.groupby('g').v.mean()-d.v.mean())**2).sum()/(n-1); msw=((d.v-d.groupby('g').v.transform('mean'))**2).sum()/(N-n)
    kbar=(N-(k**2).sum()/N)/(n-1); return float((msb-msw)/(msb+(kbar-1)*msw))
OUT={'reliability_discovery_ICC':{c:round(icc1(sl[c],z.WWTPID),3) for c in sl.columns}}
# sewage cohort
m=pd.read_pickle('martiny_matrix.pkl')
meta=pd.read_excel(U+'41467_2025_66070_MOESM4_ESM.xlsx').drop_duplicates('complete_name').set_index('complete_name')
ref=pd.read_csv('bacterial_reference.csv',index_col=0).iloc[:,0]; gp=meta.genepid.to_dict()
keep=[i for i in m.index if i in meta.index and gp.get(i) in ref.index]; m=m.loc[keep]; meta=meta.loc[keep]
bac=np.array([ref[gp[i]] for i in keep],dtype=float)
S=json.load(open('frozen_spec3.json')); sd2=pd.read_excel(U+'41467_2025_66070_MOESM5_ESM.xlsx')
acq=set(sd2[sd2.database=='resfinder'].gene)
MOD=[p for v in S['gene_mapping']['module_mapped'].values() for p in v if p in m.columns]
prev=(m>0).mean(); coreS=set(prev[prev>=0.5].index)
MODs=sorted(set(MOD)&coreS); ACQs=sorted(acq&coreS)
X=m.values.copy()
for i in range(X.shape[0]):
    r=X[i]; p=r[r>0]
    if len(p): r[r==0]=0.65*p.min()
LD=pd.DataFrame(np.log(X)-np.log(bac)[:,None],index=m.index,columns=m.columns)
tot=np.log(m.sum(1).values)-np.log(bac)
sw=pd.DataFrame({'total':tot,'richness':(m>0).sum(1).values,
                 'acquired_sum':np.log(m[ACQs].sum(1).values)-np.log(bac),
                 'module_score':LD[MODs].mean(1).values,
                 'module_share':np.log(m[MODs].sum(1).values/m.sum(1).values)},index=m.index)
city=sw.groupby(meta.city.astype(str).values).mean(); ctry=meta.groupby(meta.city.astype(str)).country.first().loc[city.index]
yr=meta.year
cy={c:sw[c].groupby([meta.city.astype(str).values,yr.values]).mean() for c in sw.columns}
multi=[c for c,n in meta.groupby(meta.city.astype(str)).year.nunique().items() if n>=2]
stab={}
for c,s_ in cy.items():
    d=s_.reset_index(); d.columns=['city','year','v']; d=d[d.city.isin(multi)].sort_values('year')
    f=d.groupby('city').first().v; l=d.groupby('city').last().v
    r,p=spearmanr(f,l); stab[c]=[round(float(r),3),float(p)]
OUT['stability_sewage_rho']=stab
F=pd.read_csv('F_country_table.csv')
resp={}
for c in sw.columns:
    cs=sw[c].groupby(meta.country.astype(str).values).mean(); v=F.country.map(cs)
    r,p=spearmanr(v,F.ddd); resp[c]=[round(float(r),3),float(p)]
OUT['response_to_consumption_rho']=resp
def r2c(v):
    v=pd.Series(v,index=city.index); m_=v.groupby(ctry.values).transform('mean')
    return float(1-((v-m_)**2).sum()/((v-v.mean())**2).sum())
OUT['country_R2_sewage']={c:round(r2c(city[c]),3) for c in sw.columns}
OUT['n']={'discovery_samples':len(z),'sewage_samples':len(m),'module_genes_discovery':len(MODg),'module_genes_sewage':len(MODs),'cities_multi_year':len(multi),'countries_consumption':len(F)}
json.dump(OUT,open('benchmark.json','w'),indent=1); print(json.dumps(OUT,indent=1))
