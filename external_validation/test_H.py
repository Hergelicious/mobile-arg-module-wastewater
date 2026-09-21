"""Test H (frozen spec 3) plus analyses added once collection dates proved available.
H primary (pre-specified): D after residualising each gene on run-level factors that vary in the SRA run table
  (read length 301 vs 302 bp; release batch). Instrument, centre and layout are constant.
Added: (1) D after residualising on collection date (uninformative by design: country explains 97% of date variance);
       (2) within-country date sensitivity in the USA and China: mean |rho| of mobile vs other genes with collection
           date, against 999 permutations of dates among plants within each country."""
import warnings; warnings.filterwarnings("ignore")
import json, numpy as np, pandas as pd
from scipy.stats import rankdata
rng=np.random.default_rng(424242); U='/mnt/user-data/uploads/'
z=pd.read_excel(U+'41467_2025_59019_MOESM8_ESM.xlsx','Fig. S2b and Fig. S4',index_col=0)
r=pd.read_csv(U+'SraRunTable.csv',low_memory=False); w=r[r['Assay Type']=='WGS'].set_index('Sample Name').loc[z.index]
g=pd.read_csv('gene_table_final.csv',index_col=0); core=list(g.index); M=g.M_prior.values.astype(bool)
G=z[core].astype(float); LP=np.log10(G+G.replace(0,np.nan).min()/2).groupby(z.WWTPID).mean()
pm=z.groupby('WWTPID').agg(country=('CountryRegion_full','first'),city=('City','first')).loc[LP.index]
K=json.load(open('canonical.json'))
def r2(Y,H):
    Yc=Y-Y.mean(0); P=H@np.linalg.pinv(H.T@H)@H.T; return 1-((Yc-P@Yc)**2).sum(0)/(Yc**2).sum(0)
dum=lambda l: pd.get_dummies(pd.Series(l)).values.astype(float)
assert abs((lambda o:o[M].mean()-o[~M].mean())(r2(LP.values,dum(pm.country)))-K['D'][0])<1e-12
def resid(Y,Xc):
    X=np.c_[np.ones(len(Xc)),Xc]; return Y-X@np.linalg.lstsq(X,Y,rcond=None)[0]
cc=pm.groupby('city').country.first()
def Dperm(Y,n=999):
    o=r2(Y,dum(pm.country)); d0=o[M].mean()-o[~M].mean(); nl=np.empty(n)
    for i in range(n):
        lab=pd.Series(rng.permutation(cc.values),index=cc.index); oo=r2(Y,dum(pm.city.map(lab))); nl[i]=oo[M].mean()-oo[~M].mean()
    return [float(d0),float((np.sum(nl>=d0)+1)/(n+1))]
OUT={}
run=pd.DataFrame({'len301':(w.AvgSpotLen.astype(str)=='301').astype(float),'rel2':(w.ReleaseDate.astype(str).str.startswith('2021-06-29')).astype(float)},index=z.index).groupby(z.WWTPID).mean().loc[LP.index]
OUT['run_factors']={'plants_with_301bp_runs':int((run.len301>0).sum()),'plants_in_second_release':int((run.rel2>0).sum())}
OUT['H_primary_run_level']=Dperm(resid(LP.values,run.values))
dt=pd.to_datetime(w.Collection_Date,errors='coerce'); day=(dt-dt.min()).dt.days
dayp=day.groupby(z.WWTPID).mean().loc[LP.index]
ok=dayp.notna().values
OUT['date_country_R2']=float(1-((dayp-dayp.groupby(pm.country).transform('mean'))**2).sum()/((dayp-dayp.mean())**2).sum())
Yd=LP.values.copy(); Yd[ok]=resid(LP.values[ok],dayp.values[ok][:,None])
OUT['added_D_after_date_adjustment']=Dperm(Yd)
res={}
for c in ['USA','China']:
    k=(pm.country==c).values & ok
    Yc=LP.values[k]; dc=dayp.values[k]
    Ry=np.column_stack([rankdata(Yc[:,j]) for j in range(Yc.shape[1])]); Ry=(Ry-Ry.mean(0))/Ry.std(0)
    def rho(dv):
        rd=rankdata(dv); rd=(rd-rd.mean())/rd.std(); return (Ry*rd[:,None]).mean(0)
    ob=np.abs(rho(dc)); st=ob[M].mean()-ob[~M].mean()
    nl=np.array([(lambda a:a[M].mean()-a[~M].mean())(np.abs(rho(rng.permutation(dc)))) for _ in range(999)])
    res[c]={'n_plants':int(k.sum()),'n_dates':int(len(set(dc))),'mean_abs_rho_mobile':float(ob[M].mean()),'mean_abs_rho_other':float(ob[~M].mean()),
            'difference':float(st),'P':float((np.sum(nl>=st)+1)/1000)}
OUT['added_within_country_date_sensitivity']=res
json.dump(OUT,open('test_H.json','w'),indent=1); print(json.dumps(OUT,indent=1))
