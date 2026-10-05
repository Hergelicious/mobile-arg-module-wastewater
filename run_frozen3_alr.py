"""Frozen spec 3 tests with the PRE-SPECIFIED normalisation (additive log-ratio against bacterial mOTU fragments),
now that kingdom_motus_pad_agg.csv is available; plus test F (drivers).
Documented details: the length-adjusted bacterial count is inf for 4 records, so raw bacterial fragment counts are used
as the reference. Test F uses WHO GLASS consumption (via Our World in Data; antibiotics + antituberculosis drugs,
2016-2023, 73 countries) because Browne et al. 2021 estimates were not obtainable; the year nearest to each country's
median sampling year is used (the 'capped at 2018' clause referred to Browne's range and no longer applies). The GDP
partial correlation is not run: no GDP data exist for the sewage cohort in the files available."""
import warnings; warnings.filterwarnings("ignore")
import json, hashlib, re, numpy as np, pandas as pd
from scipy.stats import spearmanr
import statsmodels.formula.api as smf
txt=open('frozen_spec3.json').read()
assert hashlib.sha256(txt.encode()).hexdigest()=='b6d1fadaa45434f52ea5071119065cef7bbc48dd98d231c3365a59936ca3e996','spec changed!'
S=json.loads(txt); rng=np.random.default_rng(20260920)
U='/mnt/user-data/uploads/'
m=pd.read_pickle('martiny_matrix.pkl')
meta=pd.read_excel(U+'41467_2025_66070_MOESM4_ESM.xlsx').drop_duplicates('complete_name').set_index('complete_name')
ref=pd.read_csv('bacterial_reference.csv',index_col=0).iloc[:,0]
gp=meta.genepid.to_dict()
keep=[i for i in m.index if i in meta.index and gp.get(i) in ref.index]; m=m.loc[keep]; meta=meta.loc[keep]
bac=np.array([ref[gp[i]] for i in keep],dtype=float)
sd2=pd.read_excel(U+'41467_2025_66070_MOESM5_ESM.xlsx')
acq_all=set(sd2[sd2.database=='resfinder'].gene); fg_all=set(sd2[sd2.database!='resfinder'].gene)
MOD=sorted({p for v in S['gene_mapping']['module_mapped'].values() for p in v if p in m.columns})
prev=(m>0).mean(); core=set(prev[prev>=0.5].index)
ACQ=sorted((acq_all&core)-set(MOD)); FG=sorted((fg_all-acq_all)&core); MODc=[g for g in MOD if g in core]  # latent = not in ResFinder (IDs can be listed under several databases)
OUT={'spec_sha256':hashlib.sha256(txt.encode()).hexdigest(),'normalisation':'ALR vs raw bacterial mOTU fragments (pre-specified)',
     'n_samples':len(m),'n_cities':int(meta.city.nunique()),'n_countries':int(meta.country.nunique()),
     'n_module_genes_detected':len(MODc),'n_acquired_comparison':len(ACQ),'n_functional_comparison':len(FG)}
X=m.values.copy()
for i in range(X.shape[0]):
    row=X[i]; pos=row[row>0]
    if len(pos): row[row==0]=0.65*pos.min()
LD=pd.DataFrame(np.log(X)-np.log(bac)[:,None],index=m.index,columns=m.columns)
city=LD.groupby(meta.city.values).mean(); ctry=meta.groupby('city').country.first().loc[city.index]
def r2(Y,H):
    Yc=Y-Y.mean(0); P=H@np.linalg.pinv(H.T@H)@H.T; return 1-((Yc-P@Yc)**2).sum(0)/(Yc**2).sum(0)
dum=lambda l: pd.get_dummies(pd.Series(l)).values.astype(float)
cols=MODc+ACQ+FG; assert len(set(cols))==len(cols), "gene sets overlap"; Y=city[cols].values
mm_=np.array([c in set(MODc) for c in cols]); ma=np.array([c in set(ACQ) for c in cols]); mf=np.array([c in set(FG) for c in cols])
obs=r2(Y,dum(ctry.values)); NP=9999; na=np.empty(NP); nf=np.empty(NP)
for i in range(NP):
    o=r2(Y,dum(rng.permutation(ctry.values))); na[i]=o[mm_].mean()-o[ma].mean(); nf[i]=o[mm_].mean()-o[mf].mean()
da=obs[mm_].mean()-obs[ma].mean(); df_=obs[mm_].mean()-obs[mf].mean()
OUT['D_validation_primary']={'D_ext':float(da),'P':float((np.sum(na>=da)+1)/(NP+1)),'mean_R2_module':float(obs[mm_].mean()),'mean_R2_acquired_other':float(obs[ma].mean())}
OUT['D_validation_secondary_vs_FG']={'D_ext':float(df_),'P':float((np.sum(nf>=df_)+1)/(NP+1)),'mean_R2_FG':float(obs[mf].mean())}
W=LD[cols]-LD[cols].groupby(meta.country.values).transform('mean'); Cm=np.corrcoef(W.values.T); np.fill_diagonal(Cm,np.nan)
mpc=lambda idx: float(np.nanmean(Cm[np.ix_(idx,idx)]))
mi=np.where(mm_)[0]; ai=np.where(ma)[0]; ab=np.log10(m[cols].mean().values); null=[]
for _ in range(1000):
    pick=[(rng.choice(ai[np.abs(ab[ai]-ab[i])<0.5]) if np.any(np.abs(ab[ai]-ab[i])<0.5) else rng.choice(ai)) for i in mi]
    null.append(mpc(np.array(sorted(set(pick)))))
orr=mpc(mi); OUT['covariation']={'mean_within_country_r_module':orr,'null_median':float(np.median(null)),'P':float((np.sum(np.array(null)>=orr)+1)/1001)}
score=LD[MODc].mean(1)
d=pd.DataFrame({'score':score.values,'city':meta.city.astype(str).values,'country':meta.country.astype(str).values,
                'region':meta.Region.astype(str).values,'year':meta.year.values}).dropna()
cy=d.groupby(['city','year']).score.mean().reset_index(); multi=cy.groupby('city').year.nunique(); multi=multi[multi>=2].index
srt=cy[cy.city.isin(multi)].sort_values('year'); first=srt.groupby('city').first().score; last=srt.groupby('city').last().score
rho,p=spearmanr(first,last); OUT['E_temporal']={'n_cities_multi_year':int(len(multi)),'rho_first_vs_last_year':float(rho),'P':float(p)}
d2=d.copy(); d2['year']=d2.year.astype(int).astype(str)
md=smf.mixedlm('score ~ 1',d2,groups='country',vc_formula={'city':'0+C(city)','year':'0+C(year)'}).fit(reml=True)
comp=dict(zip(md.model.exog_vc.names,[float(v) for v in md.vcomp])); cr=np.asarray(md.cov_re); comp['country']=float(cr.ravel()[0]) if cr.size else 0.0; comp['resid']=float(md.scale)
tot=sum(comp.values()); OUT['E_temporal']['variance_fractions']={k:round(v/tot,3) for k,v in comp.items()}
# G (share is a ratio of fragment counts, so it does not depend on the normalisation)
share=(m[MODc].sum(1)/m.sum(1)).groupby(meta.city.astype(str).values).mean()
z=pd.read_excel(U+'41467_2025_59019_MOESM8_ESM.xlsx','Fig. S2b and Fig. S4',index_col=0)
lay=pd.read_csv('sample_layers.csv',index_col=0).loc[z.index]; sl=lay.mobile_fraction.groupby(z.City.astype(str).values).mean()
nm=lambda s: re.sub(r'[^a-z]','',str(s).lower()); A={nm(k):v for k,v in share.items()}; B={nm(k):v for k,v in sl.items()}; com=sorted(set(A)&set(B))
rg,pg=spearmanr([A[c] for c in com],[B[c] for c in com]); OUT['G_cross_compartment']={'n_cities':len(com),'rho':float(rg),'P':float(pg)}
# F: drivers
cons=pd.read_csv(U+'antibiotic-consumption-rate.csv'); vcol=cons.columns[-1]
cs=d.groupby('country').score.mean(); yr=d.groupby('country').year.median()
alias={'unitedstates':'unitedstatesofamerica','russia':'russianfederation'}
cmap={}
for ent,g in cons.groupby('Entity'): cmap[nm(ent)]=g.set_index('Year')[vcol]
rows=[]
for c in cs.index:
    k=nm(c); k=k if k in cmap else alias.get(k,k)
    if k in cmap:
        s=cmap[k]; y=s.index[np.argmin(np.abs(s.index-yr[c]))]; rows.append((c,cs[c],float(s[y]),int(y),d[d.country==c].region.iloc[0]))
F=pd.DataFrame(rows,columns=['country','score','ddd','year_used','region'])
rf,pf_=spearmanr(F.score,F.ddd)
OUT['F_drivers']={'source':'WHO GLASS via Our World in Data (deviation: Browne et al. 2021 not obtainable)','n_countries_matched':len(F),
                  'rho':float(rf),'P':float(pf_)}
try:
    mf2=smf.mixedlm('score ~ ddd',F,groups='region').fit(reml=True)
    OUT['F_drivers']['mixed_region']={'coef':float(mf2.params['ddd']),'P':float(mf2.pvalues['ddd'])}
except Exception as e: OUT['F_drivers']['mixed_region']=str(e)
OUT['F_drivers']['gdp_partial']='not run: no GDP data for the sewage cohort in the files available'
F.to_csv('F_country_table.csv',index=False)
json.dump(OUT,open('frozen_results3_alr.json','w'),indent=1,allow_nan=False); print(json.dumps(OUT,indent=1))
