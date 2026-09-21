import warnings; warnings.filterwarnings("ignore")
import json, re, numpy as np, pandas as pd
from scipy.stats import spearmanr
import statsmodels.formula.api as smf
U='/mnt/user-data/uploads/'
OUT=json.load(open('frozen_results3.json'))
S=json.load(open('frozen_spec3.json'))
m=pd.read_pickle('martiny_matrix.pkl')
meta=pd.read_excel(U+'41467_2025_66070_MOESM4_ESM.xlsx').drop_duplicates(subset='complete_name').set_index('complete_name')
keep=[i for i in m.index if i in meta.index]; m=m.loc[keep]; meta=meta.loc[keep]
MOD=sorted({p for v in S['gene_mapping']['module_mapped'].values() for p in v if p in m.columns})
prev=(m>0).mean(); MODc=[g for g in MOD if prev[g]>=0.5]
X=m.values.copy()
for i in range(X.shape[0]):
    row=X[i]; pos=row[row>0]
    if len(pos): row[row==0]=0.65*pos.min()
L=np.log(X); L=L-L.mean(1,keepdims=True); LD=pd.DataFrame(L,index=m.index,columns=m.columns)
df=pd.DataFrame({'score':LD[MODc].mean(1).values,'city':meta.city.astype(str).values,
                 'country':meta.country.astype(str).values,'year':meta.year.astype('Int64').astype(str).values}).dropna()
md=smf.mixedlm('score ~ 1',df,groups='country',vc_formula={'city':'0+C(city)','year':'0+C(year)'}).fit(reml=True)
comp=dict(zip(md.model.exog_vc.names,[float(v) for v in md.vcomp]))
cr=np.asarray(md.cov_re); comp['country']=float(cr.ravel()[0]) if cr.size else 0.0; comp['resid']=float(md.scale)
tot=sum(comp.values()); OUT['E_temporal']['variance_fractions']={k:round(v/tot,3) for k,v in comp.items()}
# G: cross-compartment, 12 shared cities
share_sew=(m[MODc].sum(1)/m.sum(1)).groupby(meta.city.astype(str).values).mean()
z=pd.read_excel(U+'41467_2025_59019_MOESM8_ESM.xlsx','Fig. S2b and Fig. S4',index_col=0)
lay=pd.read_csv('sample_layers.csv',index_col=0).loc[z.index]
share_sludge=lay.mobile_fraction.groupby(z.City.astype(str).values).mean()
norm=lambda s: re.sub(r'[^a-z]','',str(s).lower())
a={norm(k):v for k,v in share_sew.items()}; b={norm(k):v for k,v in share_sludge.items()}
common=sorted(set(a)&set(b))
if len(common)>=5:
    rho,p=spearmanr([a[c] for c in common],[b[c] for c in common])
    OUT['G_cross_compartment']={'n_cities':len(common),'cities':common,'rho':float(rho),'P':float(p),
        'status':'exploratory: different cohorts, years and pipelines'}
else: OUT['G_cross_compartment']={'n_cities':len(common),'status':'too few shared cities'}
json.dump(OUT,open('frozen_results3.json','w'),indent=1)
print(json.dumps({k:OUT[k] for k in ['E_temporal','G_cross_compartment']},indent=1))
