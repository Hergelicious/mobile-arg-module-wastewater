"""ADDED after frozen test C: is the class-matched association specific, or does any ARG class track any phenotype?
Every class x phenotype combination, unadjusted and adjusted for log GDP per capita and WHO region."""
import warnings; warnings.filterwarnings("ignore")
import json, re, numpy as np, pandas as pd, statsmodels.api as sm
from scipy.stats import spearmanr
U='/mnt/user-data/uploads/'; nm=lambda s: re.sub(r'[^a-z]','',str(s).lower())
d=pd.read_csv('/tmp/glass/GLASS2022-master/compiled_WHO_GLASS_2022.csv'); b=d[d.Specimen=='BLOOD']
PH={'ecoli_3GC':('Escherichia coli',['Cefotaxime','Ceftriaxone','Ceftazidime']),'ecoli_cipro':('Escherichia coli',['Ciprofloxacin']),
    'kpneu_3GC':('Klebsiella pneumoniae',['Cefotaxime','Ceftriaxone','Ceftazidime']),'kpneu_carb':('Klebsiella pneumoniae',['Meropenem','Imipenem'])}
m=pd.read_pickle('martiny_matrix.pkl')
meta=pd.read_excel(U+'41467_2025_66070_MOESM4_ESM.xlsx').drop_duplicates('complete_name').set_index('complete_name')
ref=pd.read_csv('bacterial_reference.csv',index_col=0).iloc[:,0]; gp=meta.genepid.to_dict()
keep=[i for i in m.index if i in meta.index and gp.get(i) in ref.index]; m=m.loc[keep]; meta=meta.loc[keep]
bac=np.array([ref[gp[i]] for i in keep],dtype=float)
sd2=pd.read_excel(U+'41467_2025_66070_MOESM5_ESM.xlsx'); acq=sd2[sd2.database=='resfinder']; inmat=set(m.columns)
carb=acq[acq.fa_name.astype(str).str.contains('ndm|kpc|oxa-48|oxa-181|oxa-232|vim-|imp-',case=False,regex=True)]
CLS={'beta_lactam':sorted(set(acq[acq['class']=='beta_lactam'].gene)&inmat),'quinolone':sorted(set(acq[acq['class']=='quinolone'].gene)&inmat),
     'carbapenemase':sorted(set(carb.gene)&inmat),'aminoglycoside':sorted(set(acq[acq['class']=='aminoglycoside'].gene)&inmat),
     'tetracycline':sorted(set(acq[acq['class']=='tetracycline'].gene)&inmat),'macrolide_MLS':sorted(set(acq[acq['class'].astype(str).str.contains('macrolide')].gene)&inmat)}
cn=pd.Series([nm(c) for c in meta.country],index=m.index); region=meta.groupby(cn.values).Region.first()
w=pd.read_csv('/tmp/wdi/wdi_subset.csv'); gdp=w[w['Indicator Code']=='NY.GDP.PCAP.PP.KD'].set_index('Country Name')
yrs=[c for c in gdp.columns if c.isdigit() and int(c)<=2021]; gdp=gdp[yrs].ffill(axis=1)[yrs[-1]]; gdp.index=[nm(i) for i in gdp.index]
score={k:pd.Series(np.log(m[v].sum(1).values+0.5)-np.log(bac),index=m.index).groupby(cn.values).mean() for k,v in CLS.items() if len(v)>=3}
OUT={'matched_pairs':{'ecoli_3GC':'beta_lactam','ecoli_cipro':'quinolone','kpneu_3GC':'beta_lactam','kpneu_carb':'carbapenemase'},'grid':{}}
for ph,(pa,dr) in PH.items():
    s=b[(b.PathogenName==pa)&(b.AbTargets.isin(dr))].groupby('CountryTerritoryArea')[['Resistant','InterpretableAST']].sum()
    s=s[s.InterpretableAST>=20]; res=(s.Resistant/s.InterpretableAST); res.index=[nm(i) for i in res.index]
    OUT['grid'][ph]={}
    for cl,sc in score.items():
        df=pd.DataFrame({'res':res}).join(pd.DataFrame({'cls':sc,'gdp':gdp,'region':region}),how='inner').dropna()
        rho,p=spearmanr(df.cls,df.res)
        y=np.log(df.res.clip(0.01,0.99)/(1-df.res.clip(0.01,0.99)))
        Xd=pd.get_dummies(df.region,drop_first=True).astype(float); Xd['loggdp']=np.log(df.gdp); Xd['cls']=df.cls
        f=sm.OLS(y,sm.add_constant(Xd)).fit()
        OUT['grid'][ph][cl]={'rho':round(float(rho),3),'P':float(p),'adj_coef':round(float(f.params['cls']),3),'adj_P':float(f.pvalues['cls']),'n':int(len(df))}
json.dump(OUT,open('class_placebo.json','w'),indent=1)
g=pd.DataFrame({ph:{cl:v['rho'] for cl,v in d_.items()} for ph,d_ in OUT['grid'].items()})
ga=pd.DataFrame({ph:{cl:v['adj_P'] for cl,v in d_.items()} for ph,d_ in OUT['grid'].items()})
print('Spearman rho (rows = ARG class, columns = phenotype):'); print(g.round(2).to_string())
print('\nincome-adjusted P:'); print(ga.round(3).to_string())
