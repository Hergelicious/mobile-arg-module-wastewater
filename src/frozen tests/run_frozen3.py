"""Frozen Tier A tests D, E, G in the Martiny sewage cohort.
DOCUMENTED DEVIATION: the specification asks for additive log-ratio normalisation with bacterial (mOTU) fragments as
reference; that file was not available, so centred log-ratio (CLR) over the ARG composition is used instead. Zeros are
replaced by 65% of the per-sample detection limit, as in Martiny et al. Genes are restricted to those detected in >=50%
of samples (as for the 'core' genes in the discovery cohort) to keep the comparison sets well measured."""
import warnings; warnings.filterwarnings("ignore")
import json, hashlib, numpy as np, pandas as pd
from scipy.stats import spearmanr
txt=open('frozen_spec3.json').read()
assert hashlib.sha256(txt.encode()).hexdigest()=='b6d1fadaa45434f52ea5071119065cef7bbc48dd98d231c3365a59936ca3e996','spec changed!'
S=json.loads(txt); rng=np.random.default_rng(20260920)
U='/mnt/user-data/uploads/'
m=pd.read_pickle('martiny_matrix.pkl')
meta=pd.read_excel(U+'41467_2025_66070_MOESM4_ESM.xlsx')
meta=meta.drop_duplicates(subset='complete_name',keep='first').set_index('complete_name')  # 2 duplicated records (city-name spelling)
keep=[i for i in m.index if i in meta.index]; m=m.loc[keep]; meta=meta.loc[keep]
sd2=pd.read_excel(U+'41467_2025_66070_MOESM5_ESM.xlsx')
acq_all=set(sd2[sd2.database=='resfinder'].gene); fg_all=set(sd2[sd2.database!='resfinder'].gene)
MOD=sorted({p for v in S['gene_mapping']['module_mapped'].values() for p in v if p in m.columns})
prev=(m>0).mean()
core=set(prev[prev>=0.5].index)
ACQ=sorted((acq_all&core)-set(MOD)); FG=sorted((fg_all-acq_all)&core); MODc=[g for g in MOD if g in core]  # latent = not in ResFinder (IDs can be listed under several databases)
OUT={'spec_sha256':hashlib.sha256(txt.encode()).hexdigest(),'deviation':'CLR instead of ALR (mOTU file unavailable)',
     'n_samples':len(m),'n_cities':int(meta.city.nunique()),'n_countries':int(meta.country.nunique()),
     'n_module_genes_detected':len(MODc),'n_acquired_comparison':len(ACQ),'n_functional_comparison':len(FG)}
# CLR with 65% detection-limit zero replacement
X=m.values.copy()
for i in range(X.shape[0]):
    row=X[i]; pos=row[row>0]
    if len(pos): row[row==0]=0.65*pos.min()
L=np.log(X); L=L-L.mean(1,keepdims=True)
LD=pd.DataFrame(L,index=m.index,columns=m.columns)
city=LD.groupby(meta.city.values).mean()
ctry=meta.groupby('city').country.first().loc[city.index]
def r2(Y,H):
    Yc=Y-Y.mean(0); P=H@np.linalg.pinv(H.T@H)@H.T
    return 1-((Yc-P@Yc)**2).sum(0)/(Yc**2).sum(0)
dum=lambda l: pd.get_dummies(pd.Series(l)).values.astype(float)
cols=MODc+ACQ+FG; assert len(set(cols))==len(cols), "gene sets overlap"; Y=city[cols].values
mask_mod=np.array([c in set(MODc) for c in cols]); mask_acq=np.array([c in set(ACQ) for c in cols]); mask_fg=np.array([c in set(FG) for c in cols])
obs=r2(Y,dum(ctry.values))
NP=9999; nl_a=np.empty(NP); nl_f=np.empty(NP)
for i in range(NP):
    p=rng.permutation(ctry.values); o=r2(Y,dum(p))
    nl_a[i]=o[mask_mod].mean()-o[mask_acq].mean(); nl_f[i]=o[mask_mod].mean()-o[mask_fg].mean()
d_a=obs[mask_mod].mean()-obs[mask_acq].mean(); d_f=obs[mask_mod].mean()-obs[mask_fg].mean()
OUT['D_validation_primary']={'D_ext':float(d_a),'P':float((np.sum(nl_a>=d_a)+1)/(NP+1)),
    'mean_R2_module':float(obs[mask_mod].mean()),'mean_R2_acquired_other':float(obs[mask_acq].mean())}
OUT['D_validation_secondary_vs_FG']={'D_ext':float(d_f),'P':float((np.sum(nl_f>=d_f)+1)/(NP+1)),'mean_R2_FG':float(obs[mask_fg].mean())}
# within-country co-variation vs abundance-matched acquired sets
W=LD[cols]-LD[cols].groupby(meta.country.values).transform('mean')
C=np.corrcoef(W.values.T); np.fill_diagonal(C,np.nan)
def mpc(idx): 
    sub=C[np.ix_(idx,idx)]; return float(np.nanmean(sub))
mod_idx=np.where(mask_mod)[0]; acq_idx=np.where(mask_acq)[0]
ab=np.log10(m[cols].mean().values)
null=[]
for _ in range(1000):
    pick=[]
    for i in mod_idx:
        cand=acq_idx[np.abs(ab[acq_idx]-ab[i])<0.5]
        pick.append(rng.choice(cand) if len(cand) else rng.choice(acq_idx))
    null.append(mpc(np.array(sorted(set(pick)))))
obs_r=mpc(mod_idx)
OUT['covariation']={'mean_within_country_r_module':obs_r,'null_median':float(np.median(null)),
                    'P':float((np.sum(np.array(null)>=obs_r)+1)/1001)}
# E: temporal stability
score=LD[MODc].mean(1)
df=pd.DataFrame({'score':score.values,'city':meta.city.values,'country':meta.country.values,'year':meta.year.values}).dropna()
cy=df.groupby(['city','year']).score.mean().reset_index()
multi=cy.groupby('city').year.nunique(); multi=multi[multi>=2].index
first=cy[cy.city.isin(multi)].sort_values('year').groupby('city').first().score
last=cy[cy.city.isin(multi)].sort_values('year').groupby('city').last().score
rho,p=spearmanr(first,last)
OUT['E_temporal']={'n_cities_multi_year':int(len(multi)),'rho_first_vs_last_year':float(rho),'P':float(p)}
import statsmodels.formula.api as smf
try:
    md=smf.mixedlm('score ~ 1',df,groups=df.country,vc_formula={'city':'0+C(city)','year':'0+C(year)'}).fit(reml=True)
    comp=dict(zip(md.model.exog_vc.names,[float(v) for v in md.vcomp])); comp['country']=float(md.cov_re.iloc[0,0]); comp['resid']=float(md.scale)
    tot=sum(comp.values()); OUT['E_temporal']['variance_fractions']={k:round(v/tot,3) for k,v in comp.items()}
except Exception as e: OUT['E_temporal']['variance_fractions']=str(e)
json.dump(OUT,open('frozen_results3.json','w'),indent=1); print(json.dumps(OUT,indent=1))
