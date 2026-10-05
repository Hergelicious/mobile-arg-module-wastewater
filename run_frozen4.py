import warnings; warnings.filterwarnings("ignore")
import json, hashlib, re, numpy as np, pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr
txt=open('frozen_spec4.json').read()
assert hashlib.sha256(txt.encode()).hexdigest()=='efc0f84a3cb5b40dec070c8dfab14f0938f5afc29076137d84039753b394b591','spec changed!'
U='/mnt/user-data/uploads/'; nm=lambda s: re.sub(r'[^a-z]','',str(s).lower())
# --- clinical outcomes
d=pd.read_csv('/tmp/glass/GLASS2022-master/compiled_WHO_GLASS_2022.csv')
b=d[d.Specimen=='BLOOD']
def outcome(path,drugs):
    s=b[(b.PathogenName==path)&(b.AbTargets.isin(drugs))].groupby('CountryTerritoryArea')[['Resistant','InterpretableAST']].sum()
    s=s[s.InterpretableAST>=20]; return (s.Resistant/s.InterpretableAST).rename('res'), s.InterpretableAST
OUTC={'ecoli_3GC':outcome('Escherichia coli',['Cefotaxime','Ceftriaxone','Ceftazidime']),
      'ecoli_cipro':outcome('Escherichia coli',['Ciprofloxacin']),
      'kpneu_3GC':outcome('Klebsiella pneumoniae',['Cefotaxime','Ceftriaxone','Ceftazidime']),
      'kpneu_carb':outcome('Klebsiella pneumoniae',['Meropenem','Imipenem'])}
# --- wastewater predictors (sewage cohort, ALR as pre-specified)
m=pd.read_pickle('martiny_matrix.pkl')
meta=pd.read_excel(U+'41467_2025_66070_MOESM4_ESM.xlsx').drop_duplicates('complete_name').set_index('complete_name')
ref=pd.read_csv('bacterial_reference.csv',index_col=0).iloc[:,0]; gp=meta.genepid.to_dict()
keep=[i for i in m.index if i in meta.index and gp.get(i) in ref.index]; m=m.loc[keep]; meta=meta.loc[keep]
bac=np.array([ref[gp[i]] for i in keep],dtype=float)
S3=json.load(open('frozen_spec3.json')); sd2=pd.read_excel(U+'41467_2025_66070_MOESM5_ESM.xlsx')
acq=set(sd2[sd2.database=='resfinder'].gene)
MOD=[p for v in S3['gene_mapping']['module_mapped'].values() for p in v if p in m.columns]
prev=(m>0).mean(); core=set(prev[prev>=0.5].index)
MODs=sorted(set(MOD)&core); ACQs=sorted((acq&core)-set(MODs))
X=m.values.copy().astype(float)
for i in range(X.shape[0]):
    r=X[i]; p=r[r>0]
    if len(p): r[r==0]=0.65*p.min()
LD=pd.DataFrame(np.log(X)-np.log(bac)[:,None],index=m.index,columns=m.columns)
cn=meta.country.astype(str).values
W=pd.DataFrame({'module':LD[MODs].mean(1).groupby(cn).mean(),'acquired_other':LD[ACQs].mean(1).groupby(cn).mean(),
                'total':(np.log(m.sum(1).values)-np.log(bac)).groupby(cn).mean() if False else pd.Series(np.log(m.sum(1).values)-np.log(bac),index=m.index).groupby(cn).mean(),
                'richness':pd.Series((m>0).sum(1).values,index=m.index).groupby(cn).mean(),
                'region':meta.groupby(cn).Region.first()})
W.index.name='country'; W['key']=[nm(i) for i in W.index]
# --- covariates
cons=pd.read_csv(U+'antibiotic-consumption-rate.csv'); vcol=cons.columns[-1]
cons=cons.sort_values('Year').groupby('Entity')[vcol].last(); cons.index=[nm(i) for i in cons.index]
w=pd.read_csv('/tmp/wdi/wdi_subset.csv')
gdp=w[w['Indicator Code']=='NY.GDP.PCAP.PP.KD'].set_index('Country Name')
yrs=[c for c in gdp.columns if c.isdigit() and int(c)<=2021]
gdp=gdp[yrs].ffill(axis=1)[yrs[-1]]; gdp.index=[nm(i) for i in gdp.index]
OUT={'spec_sha256':hashlib.sha256(txt.encode()).hexdigest(),'n_countries_wastewater':len(W),'results':{}}
for name,(res,nIso) in OUTC.items():
    r=res.copy(); r.index=[nm(i) for i in r.index]
    df=W.copy(); df['res']=df.key.map(r); df['cons']=df.key.map(cons); df['gdp']=df.key.map(gdp)
    df=df.dropna(subset=['res'])
    ok=df.dropna(subset=['cons','gdp'])
    rho,p=spearmanr(df.module,df.res)
    ent={'n_countries':int(len(df)),'n_with_covariates':int(len(ok)),'T1_spearman_module':[float(rho),float(p)]}
    for pred in ['module','acquired_other','total','richness']:
        ent[f'T1_spearman_{pred}']=list(map(float,spearmanr(df[pred],df.res)))
    if len(ok)>=25:
        y=np.log(ok.res.clip(0.01,0.99)/(1-ok.res.clip(0.01,0.99)))
        base=pd.get_dummies(ok.region,drop_first=True).astype(float)
        base['logcons']=np.log(ok.cons); base['loggdp']=np.log(ok.gdp)
        m0=sm.OLS(y,sm.add_constant(base)).fit()
        for pred in ['module','acquired_other','total']:
            Xd=base.copy(); Xd[pred]=ok[pred]
            m1=sm.OLS(y,sm.add_constant(Xd)).fit()
            ent[f'T2_{pred}']={'coef':float(m1.params[pred]),'P':float(m1.pvalues[pred]),
                               'adjR2_with':float(m1.rsquared_adj),'adjR2_without':float(m0.rsquared_adj),
                               'delta_adjR2':float(m1.rsquared_adj-m0.rsquared_adj)}
    OUT['results'][name]=ent
json.dump(OUT,open('frozen_results4.json','w'),indent=1)
r=OUT['results']['ecoli_3GC']
print('PRIMARY (E. coli 3GC): n =',r['n_countries'],'| Spearman module rho = %.3f, P = %.4f'%tuple(r['T1_spearman_module']))
for k,v in OUT['results'].items():
    print(f"\n{k}: n={v['n_countries']}, with covariates {v['n_with_covariates']}")
    for pred in ['module','acquired_other','total','richness']:
        a,b_=v[f'T1_spearman_{pred}']; print(f"   T1 {pred:15s} rho={a:+.3f} P={b_:.4f}")
    for pred in ['module','acquired_other','total']:
        if f'T2_{pred}' in v:
            t=v[f'T2_{pred}']; print(f"   T2 {pred:15s} coef={t['coef']:+.3f} P={t['P']:.4f} dAdjR2={t['delta_adjR2']:+.3f}")
