import warnings; warnings.filterwarnings("ignore")
import json, hashlib, re, numpy as np, pandas as pd, statsmodels.api as sm
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
txt=open('frozen_spec5.json').read()
assert hashlib.sha256(txt.encode()).hexdigest()=='8827f39982ab77ecdec194d1893345fd1f1f6d96f666490b0121019278d22948','spec changed!'
U='/mnt/user-data/uploads/'; nm=lambda s: re.sub(r'[^a-z]','',str(s).lower())
d=pd.read_csv('/tmp/glass/GLASS2022-master/compiled_WHO_GLASS_2022.csv'); b=d[d.Specimen=='BLOOD']
PH={'ecoli_3GC':('Escherichia coli',['Cefotaxime','Ceftriaxone','Ceftazidime']),'ecoli_cipro':('Escherichia coli',['Ciprofloxacin']),
    'kpneu_3GC':('Klebsiella pneumoniae',['Cefotaxime','Ceftriaxone','Ceftazidime']),'kpneu_carb':('Klebsiella pneumoniae',['Meropenem','Imipenem'])}
def cy_table(path,drugs):
    s=b[(b.PathogenName==path)&(b.AbTargets.isin(drugs))]
    t=s.groupby([s.CountryTerritoryArea.map(nm),'Year'])[['Resistant','InterpretableAST']].sum()
    t=t[t.InterpretableAST>=20]; t['res']=t.Resistant/t.InterpretableAST; return t.reset_index().rename(columns={'CountryTerritoryArea':'key'})
# sewage measures
m=pd.read_pickle('martiny_matrix.pkl')
meta=pd.read_excel(U+'41467_2025_66070_MOESM4_ESM.xlsx').drop_duplicates('complete_name').set_index('complete_name')
ref=pd.read_csv('bacterial_reference.csv',index_col=0).iloc[:,0]; gp=meta.genepid.to_dict()
keep=[i for i in m.index if i in meta.index and gp.get(i) in ref.index]; m=m.loc[keep]; meta=meta.loc[keep]
bac=np.array([ref[gp[i]] for i in keep],dtype=float)
S3=json.load(open('frozen_spec3.json')); sd2=pd.read_excel(U+'41467_2025_66070_MOESM5_ESM.xlsx')
acq=sd2[sd2.database=='resfinder']
MODids=[p for v in S3['gene_mapping']['module_mapped'].values() for p in v if p in m.columns]
prev=(m>0).mean(); core=set(prev[prev>=0.5].index)
MODs=sorted(set(MODids)&core); ACQs=sorted(set(acq.gene)&core-set(MODs))
X=m.values.copy().astype(float)
for i in range(X.shape[0]):
    r=X[i]; p=r[r>0]
    if len(p): r[r==0]=0.65*p.min()
LD=pd.DataFrame(np.log(X)-np.log(bac)[:,None],index=m.index,columns=m.columns)
key=pd.Series([nm(c) for c in meta.country],index=m.index); yr=meta.year
def cy_mean(series): return series.groupby([key.values,yr.values]).mean()
SW=pd.DataFrame({'module':cy_mean(LD[MODs].mean(1)),'acquired_other':cy_mean(LD[ACQs].mean(1)),
                 'total':cy_mean(pd.Series(np.log(m.sum(1).values)-np.log(bac),index=m.index))}).reset_index()
SW.columns=['key','Year','module','acquired_other','total']
OUT={'spec_sha256':hashlib.sha256(txt.encode()).hexdigest(),'test_P':{},'test_C':{}}
# ---------- TEST P ----------
def fe_model(df,pred,extra_fe=None):
    y=np.log(df.res.clip(0.01,0.99)/(1-df.res.clip(0.01,0.99)))
    Xd=pd.get_dummies(df.key,drop_first=True).astype(float)
    Xd=pd.concat([Xd,pd.get_dummies(df.Year,prefix='y',drop_first=True).astype(float)],axis=1)
    if extra_fe is not None: Xd=pd.concat([Xd,pd.get_dummies(df[extra_fe],prefix='ph',drop_first=True).astype(float)],axis=1)
    Xd[pred]=df[pred].values
    res=sm.OLS(y,sm.add_constant(Xd)).fit(cov_type='cluster',cov_kwds={'groups':df.key})
    ci=res.conf_int().loc[pred]
    return {'coef':float(res.params[pred]),'P':float(res.pvalues[pred]),'CI':[float(ci[0]),float(ci[1])],'n_obs':int(len(df)),'n_countries':int(df.key.nunique())}
prim=cy_table(*PH['ecoli_3GC']).merge(SW,on=['key','Year'])
keepc=prim.key.value_counts(); prim=prim[prim.key.isin(keepc[keepc>=2].index)]
OUT['test_P']['primary']={p:fe_model(prim,p) for p in ['module','acquired_other','total']}
pool=[]
for k,(pa,dr) in PH.items():
    t=cy_table(pa,dr).merge(SW,on=['key','Year']); t['ph']=k; pool.append(t)
pool=pd.concat(pool); kc=pool.key.value_counts(); pool=pool[pool.key.isin(kc[kc>=2].index)]
OUT['test_P']['pooled']={p:fe_model(pool,p,extra_fe='ph') for p in ['module','acquired_other','total']}
# ---------- TEST C ----------
carb=acq[acq.fa_name.astype(str).str.contains('ndm|kpc|oxa-48|oxa-181|oxa-232|vim-|imp-',case=False,regex=True)]
# class scores aggregate ALL genes of the class present in the data (no prevalence filter): individual carbapenemase
# genes are rare in sewage (none reach 50% prevalence), but they still contribute to the class total.
inmat=set(m.columns)
CLS={'beta_lactam':sorted(set(acq[acq['class']=='beta_lactam'].gene)&inmat),'quinolone':sorted(set(acq[acq['class']=='quinolone'].gene)&inmat),
     'carbapenemase':sorted(set(carb.gene)&inmat)}
PAIR={'ecoli_3GC':'beta_lactam','ecoli_cipro':'quinolone','kpneu_3GC':'beta_lactam','kpneu_carb':'carbapenemase'}
w=pd.read_csv('/tmp/wdi/wdi_subset.csv'); gdp=w[w['Indicator Code']=='NY.GDP.PCAP.PP.KD'].set_index('Country Name')
yrs=[c for c in gdp.columns if c.isdigit() and int(c)<=2021]; gdp=gdp[yrs].ffill(axis=1)[yrs[-1]]; gdp.index=[nm(i) for i in gdp.index]
cn=pd.Series([nm(c) for c in meta.country],index=m.index)
region=meta.groupby(cn.values).Region.first()
Ps=[]
for ph,cl in PAIR.items():
    genes=CLS[cl]
    if len(genes)<3: OUT['test_C'][ph]={'class':cl,'n_genes':len(genes),'status':'too few genes'}; Ps.append(1.0); continue
    sc=pd.Series(np.log(m[genes].sum(1).values+0.5)-np.log(bac),index=m.index).groupby(cn.values).mean()   # class total per bacterial fragment
    modsub=sorted(set(genes)&set(MODs))
    scm=pd.Series(np.log(m[modsub].sum(1).values+0.5)-np.log(bac),index=m.index).groupby(cn.values).mean() if len(modsub)>=5 else None
    t=cy_table(*PH[ph]).groupby('key')[['Resistant','InterpretableAST']].sum()
    res=(t.Resistant/t.InterpretableAST)
    df=pd.DataFrame({'res':res}).join(pd.DataFrame({'cls':sc,'gdp':gdp,'region':region}),how='inner').dropna()
    rho,p=spearmanr(df.cls,df.res); Ps.append(p)
    y=np.log(df.res.clip(0.01,0.99)/(1-df.res.clip(0.01,0.99)))
    Xd=pd.get_dummies(df.region,drop_first=True).astype(float); Xd['loggdp']=np.log(df.gdp); Xd['cls']=df.cls
    mfit=sm.OLS(y,sm.add_constant(Xd)).fit()
    ent={'class':cl,'n_genes':len(genes),'n_countries':int(len(df)),'spearman':[float(rho),float(p)],
         'adjusted_coef':float(mfit.params['cls']),'adjusted_P':float(mfit.pvalues['cls'])}
    if scm is not None:
        df2=df.join(pd.DataFrame({'mod_cls':scm}),how='inner').dropna()
        ent['module_restricted']={'n_genes':len(modsub),'spearman':list(map(float,spearmanr(df2.mod_cls,df2.res)))}
    OUT['test_C'][ph]=ent
adj=multipletests(Ps,method='holm')[1]
for (ph,_),a in zip(PAIR.items(),adj):
    if 'status' not in OUT['test_C'][ph]: OUT['test_C'][ph]['spearman_P_holm']=float(a)
json.dump(OUT,open('frozen_results5.json','w'),indent=1)
print('TEST P (within-country, E. coli 3GC): n_obs=%d, countries=%d'%(OUT['test_P']['primary']['module']['n_obs'],OUT['test_P']['primary']['module']['n_countries']))
for p in ['module','acquired_other','total']:
    r=OUT['test_P']['primary'][p]; print(f"   {p:15s} coef={r['coef']:+.3f} [{r['CI'][0]:+.2f},{r['CI'][1]:+.2f}] P={r['P']:.3f}")
print('TEST P pooled: n_obs=%d'%OUT['test_P']['pooled']['module']['n_obs'])
for p in ['module','acquired_other','total']:
    r=OUT['test_P']['pooled'][p]; print(f"   {p:15s} coef={r['coef']:+.3f} P={r['P']:.3f}")
print('\nTEST C (class-matched):')
for ph,v in OUT['test_C'].items():
    if 'status' in v: print(f"   {ph:12s} {v['class']:14s} {v['status']}"); continue
    mr=v.get('module_restricted',{}).get('spearman',[float('nan'),float('nan')])[0]
    print(f"   {ph:12s} {v['class']:14s} n={v['n_countries']:3d} genes={v['n_genes']:4d} rho={v['spearman'][0]:+.3f} Holm={v['spearman_P_holm']:.4f} income-adj P={v['adjusted_P']:.3f} module-only rho={mr:+.3f}")
