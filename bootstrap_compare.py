"""Statistical comparison of metrics, and a confidence interval for the USA+China exclusion result.
(1) Bootstrap differences in within-plant ICC between the module score and each conventional metric (resampling
    plants, 2,000 replicates).
(2) Bootstrap differences in the correlation with antibiotic consumption and in temporal stability (resampling
    countries and cities respectively).
(3) Bootstrap 95% CI for D with the USA and China excluded, resampling city clusters, so that the non-significant
    result is reported as an interval rather than as an absence of effect.
Reads existing outputs only."""
import warnings; warnings.filterwarnings("ignore")
import json, numpy as np, pandas as pd
from scipy.stats import spearmanr
rng=np.random.default_rng(99)
U='/mnt/user-data/uploads/'
g=pd.read_csv('gene_table_final.csv',index_col=0); core=list(g.index); M=g.M_prior.values.astype(bool)
mm=pd.read_csv('module_membership.csv').set_index('gene').loc[core]; MODg=[x for x in core if bool(mm.loc[x,'in_module'])]
z=pd.read_excel(U+'41467_2025_59019_MOESM8_ESM.xlsx','Fig. S2b and Fig. S4',index_col=0)
G=z[core].astype(float); L=np.log10(G+G.replace(0,np.nan).min()/2)
met=pd.DataFrame({'total':np.log10(z.ARG_per_cell),'richness':(G>0).sum(1),
                  'acquired_sum':np.log10(G[[x for x in core if bool(g.M_prior[x])]].sum(1)),
                  'module_score':L[MODg].mean(1),'module_share':np.log10(G[MODg].sum(1)/G[core].sum(1))},index=z.index)
plants=z.WWTPID
def icc1(v,grp):
    d=pd.DataFrame({'v':v,'g':grp}); d=d[d.g.map(d.g.value_counts())>=2]
    k=d.g.value_counts(); n=len(k); N=len(d)
    if n<3: return np.nan
    msb=(k*(d.groupby('g').v.mean()-d.v.mean())**2).sum()/(n-1); msw=((d.v-d.groupby('g').v.transform('mean'))**2).sum()/(N-n)
    kbar=(N-(k**2).sum()/N)/(n-1); return float((msb-msw)/(msb+(kbar-1)*msw))
OUT={'ICC_differences':{}}
uniq=plants.unique()
boot={m:[] for m in met.columns}
for _ in range(2000):
    pick=rng.choice(uniq,len(uniq),replace=True)
    idx=np.concatenate([np.where(plants.values==p)[0] for p in pick])
    lab=np.concatenate([[f'{p}_{i}']*int((plants.values==p).sum()) for i,p in enumerate(pick)])
    for m in met.columns: boot[m].append(icc1(met[m].values[idx],pd.Series(lab)))
B=pd.DataFrame(boot).dropna()
for m in ['total','richness','acquired_sum','module_share']:
    d=B.module_score-B[m]
    OUT['ICC_differences'][f'module_score_minus_{m}']={'point':round(float(icc1(met.module_score,plants)-icc1(met[m],plants)),3),
        'CI95':[round(float(np.percentile(d,2.5)),3),round(float(np.percentile(d,97.5)),3)],
        'P_two_sided':float(2*min((d<=0).mean(),(d>=0).mean()))}
# (2) sewage: consumption and stability differences
m_=pd.read_pickle('martiny_matrix.pkl')
meta=pd.read_excel(U+'41467_2025_66070_MOESM4_ESM.xlsx').drop_duplicates('complete_name').set_index('complete_name')
ref=pd.read_csv('bacterial_reference.csv',index_col=0).iloc[:,0]; gp=meta.genepid.to_dict()
keep=[i for i in m_.index if i in meta.index and gp.get(i) in ref.index]; m_=m_.loc[keep]; meta=meta.loc[keep]
bac=np.array([ref[gp[i]] for i in keep],dtype=float)
S3=json.load(open('frozen_spec3.json')); sd2=pd.read_excel(U+'41467_2025_66070_MOESM5_ESM.xlsx')
acq=set(sd2[sd2.database=='resfinder'].gene)
MOD=[p for v in S3['gene_mapping']['module_mapped'].values() for p in v if p in m_.columns]
prev=(m_>0).mean(); coreS=set(prev[prev>=0.5].index); MODs=sorted(set(MOD)&coreS); ACQs=sorted((acq&coreS)-set(MODs))
X=m_.values.copy().astype(float)
for i in range(X.shape[0]):
    r=X[i]; p=r[r>0]
    if len(p): r[r==0]=0.65*p.min()
LD=pd.DataFrame(np.log(X)-np.log(bac)[:,None],index=m_.index,columns=m_.columns)
cn=meta.country.astype(str).values
sw=pd.DataFrame({'total':np.log(m_.sum(1).values)-np.log(bac),'richness':(m_>0).sum(1).values,
                 'acquired_sum':np.log(m_[ACQs].sum(1).values)-np.log(bac),'module_score':LD[MODs].mean(1).values,
                 'module_share':np.log(m_[MODs].sum(1).values/m_.sum(1).values)},index=m_.index)
F=pd.read_csv('F_country_table.csv')
cs={c:sw[c].groupby(cn).mean() for c in sw.columns}
obs={c:spearmanr(F.country.map(cs[c]),F.ddd)[0] for c in sw.columns}
bd={c:[] for c in sw.columns}
for _ in range(2000):
    pick=rng.choice(F.index.values,len(F),replace=True); f=F.loc[pick]
    for c in sw.columns: bd[c].append(spearmanr(f.country.map(cs[c]),f.ddd)[0])
BD=pd.DataFrame(bd)
OUT['consumption_rho_differences']={}
for m in ['total','richness','module_share']:
    d=BD.module_score-BD[m]
    OUT['consumption_rho_differences'][f'module_score_minus_{m}']={'point':round(float(obs['module_score']-obs[m]),3),
        'CI95':[round(float(np.percentile(d,2.5)),3),round(float(np.percentile(d,97.5)),3)],'P_two_sided':float(2*min((d<=0).mean(),(d>=0).mean()))}
# (3) CI for D excluding USA and China
K=json.load(open('canonical.json'))
LPd=np.log10(G+G.replace(0,np.nan).min()/2).groupby(z.WWTPID).mean()
pm=z.groupby('WWTPID').agg(country=('CountryRegion_full','first'),city=('City','first')).loc[LPd.index]
sub=~pm.country.isin(['USA','China']).values
Y=LPd.values[sub]; ctry=pm.country.values[sub]; city=pm.city.values[sub]
def r2(Yv,lab):
    d=pd.DataFrame(Yv); mu=d.groupby(lab).transform('mean').values
    Yc=Yv-Yv.mean(0); return 1-((Yv-mu)**2).sum(0)/(Yc**2).sum(0)
obsD=float((lambda o:o[M].mean()-o[~M].mean())(r2(Y,ctry)))
cities=pd.unique(city); bs=[]
for _ in range(2000):
    pick=rng.choice(cities,len(cities),replace=True)
    idx=np.concatenate([np.where(city==c)[0] for c in pick])
    lab=np.concatenate([[f'{ctry[np.where(city==c)[0][0]]}_{i}' if False else ctry[np.where(city==c)[0][0]]]*int((city==c).sum()) for i,c in enumerate(pick)])
    o=r2(Y[idx],pd.Series(lab)); bs.append(float(o[M].mean()-o[~M].mean()))
OUT['D_excluding_USA_China']={'point':round(obsD,4),'CI95':[round(float(np.percentile(bs,2.5)),4),round(float(np.percentile(bs,97.5)),4)],
    'P_reported_in_primary':K['D_noBoth'][1],'n_plants':int(sub.sum()),'n_countries':int(pd.unique(ctry).size)}
json.dump(OUT,open('bootstrap_compare.json','w'),indent=1); print(json.dumps(OUT,indent=1))
