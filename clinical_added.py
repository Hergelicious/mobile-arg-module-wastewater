"""ADDED after the frozen test (labelled as such): the consumption covariate limits the adjusted model to ~24
countries, so here the same models are fitted with GDP per capita and WHO region only, which are available for most
countries. Question: does the association between sewage ARG abundance and clinical resistance survive adjustment for
national income, and is any of it specific to the module?"""
import warnings; warnings.filterwarnings("ignore")
import json, re, numpy as np, pandas as pd, statsmodels.api as sm
from scipy.stats import spearmanr
U='/mnt/user-data/uploads/'; nm=lambda s: re.sub(r'[^a-z]','',str(s).lower())
d=pd.read_csv('/tmp/glass/GLASS2022-master/compiled_WHO_GLASS_2022.csv'); b=d[d.Specimen=='BLOOD']
def outcome(path,drugs):
    s=b[(b.PathogenName==path)&(b.AbTargets.isin(drugs))].groupby('CountryTerritoryArea')[['Resistant','InterpretableAST']].sum()
    s=s[s.InterpretableAST>=20]; return (s.Resistant/s.InterpretableAST)
OUTC={'ecoli_3GC':outcome('Escherichia coli',['Cefotaxime','Ceftriaxone','Ceftazidime']),
      'ecoli_cipro':outcome('Escherichia coli',['Ciprofloxacin']),
      'kpneu_3GC':outcome('Klebsiella pneumoniae',['Cefotaxime','Ceftriaxone','Ceftazidime']),
      'kpneu_carb':outcome('Klebsiella pneumoniae',['Meropenem','Imipenem'])}
m=pd.read_pickle('martiny_matrix.pkl')
meta=pd.read_excel(U+'41467_2025_66070_MOESM4_ESM.xlsx').drop_duplicates('complete_name').set_index('complete_name')
ref=pd.read_csv('bacterial_reference.csv',index_col=0).iloc[:,0]; gp=meta.genepid.to_dict()
keep=[i for i in m.index if i in meta.index and gp.get(i) in ref.index]; m=m.loc[keep]; meta=meta.loc[keep]
bac=np.array([ref[gp[i]] for i in keep],dtype=float)
S3=json.load(open('frozen_spec3.json')); sd2=pd.read_excel(U+'41467_2025_66070_MOESM5_ESM.xlsx')
acq=set(sd2[sd2.database=='resfinder'].gene)
MOD=[p for v in S3['gene_mapping']['module_mapped'].values() for p in v if p in m.columns]
prev=(m>0).mean(); core=set(prev[prev>=0.5].index); MODs=sorted(set(MOD)&core); ACQs=sorted((acq&core)-set(MODs))
X=m.values.copy().astype(float)
for i in range(X.shape[0]):
    r=X[i]; p=r[r>0]
    if len(p): r[r==0]=0.65*p.min()
LD=pd.DataFrame(np.log(X)-np.log(bac)[:,None],index=m.index,columns=m.columns)
cn=meta.country.astype(str).values
W=pd.DataFrame({'module':LD[MODs].mean(1).groupby(cn).mean(),'acquired_other':LD[ACQs].mean(1).groupby(cn).mean(),
                'total':pd.Series(np.log(m.sum(1).values)-np.log(bac),index=m.index).groupby(cn).mean(),
                'region':meta.groupby(cn).Region.first()})
W['key']=[nm(i) for i in W.index]
w=pd.read_csv('/tmp/wdi/wdi_subset.csv'); gdp=w[w['Indicator Code']=='NY.GDP.PCAP.PP.KD'].set_index('Country Name')
yrs=[c for c in gdp.columns if c.isdigit() and int(c)<=2021]; gdp=gdp[yrs].ffill(axis=1)[yrs[-1]]
gdp.index=[nm(i) for i in gdp.index]
OUT={}
for name,res in OUTC.items():
    r=res.copy(); r.index=[nm(i) for i in r.index]
    df=W.copy(); df['res']=df.key.map(r); df['gdp']=df.key.map(gdp)
    df=df.dropna(subset=['res','gdp'])
    y=np.log(df.res.clip(0.01,0.99)/(1-df.res.clip(0.01,0.99)))
    base=pd.get_dummies(df.region,drop_first=True).astype(float); base['loggdp']=np.log(df.gdp)
    m0=sm.OLS(y,sm.add_constant(base)).fit()
    ent={'n':int(len(df)),'adjR2_covariates_only':float(m0.rsquared_adj),
         'gdp_only_spearman':list(map(float,spearmanr(np.log(df.gdp),df.res)))}
    for pred in ['module','acquired_other','total']:
        Xd=base.copy(); Xd[pred]=df[pred]; m1=sm.OLS(y,sm.add_constant(Xd)).fit()
        ent[pred]={'coef':float(m1.params[pred]),'P':float(m1.pvalues[pred]),'delta_adjR2':float(m1.rsquared_adj-m0.rsquared_adj)}
    OUT[name]=ent
json.dump(OUT,open('clinical_added.json','w'),indent=1)
for k,v in OUT.items():
    print(f"{k}: n={v['n']} | GDP-resistance rho={v['gdp_only_spearman'][0]:+.2f} (P={v['gdp_only_spearman'][1]:.3f})")
    for p_ in ['module','acquired_other','total']:
        t=v[p_]; print(f"    {p_:15s} coef={t['coef']:+.3f} P={t['P']:.4f}  dAdjR2={t['delta_adjR2']:+.3f}")
