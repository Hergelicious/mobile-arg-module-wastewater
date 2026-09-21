"""Independent re-derivation of every headline number, written from scratch."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, json
from scipy.stats import spearmanr, mannwhitneyu
rng=np.random.default_rng(4242)
x=pd.ExcelFile('/mnt/user-data/uploads/41467_2025_59019_MOESM8_ESM.xlsx')
s=pd.read_excel(x,'Fig. S2b and Fig. S4',index_col=0)
g=pd.read_csv('gene_table2.csv',index_col=0); core=list(g.index); M=g.M_prior.values
G=s[core].astype(float)
V={}
# --- design ---
V['n_samples']=len(s); V['n_plants']=s.WWTPID.nunique(); V['n_cities']=s.City.nunique(); V['n_countries']=s.CountryRegion_full.nunique()
V['n_mobile']=int(M.sum()); V['n_other']=int((~M).sum())
# --- load ---
lt=np.log10(s.ARG_per_cell)
V['load_sd']=float(lt.std()); V['load_IQR_fold']=float(10**(lt.quantile(.75)-lt.quantile(.25)))
for lev,col in [('country','CountryRegion_full'),('city','City'),('plant','WWTPID')]:
    V[f'load_R2_{lev}']=float(1-((lt-lt.groupby(s[col]).transform('mean'))**2).sum()/((lt-lt.mean())**2).sum())
# --- layers ---
mob=list(g.index[M]); oth=list(g.index[~M])
lmob=np.log10(G[mob].sum(1)); loth=np.log10(G[oth].sum(1))
V['layer_sd_mob']=float(lmob.std()); V['layer_sd_oth']=float(loth.std())
V['share_oth_median']=float((G[oth].sum(1)/G[core].sum(1)).median())
V['mobile_share_min_med_max']=[float(v) for v in np.percentile(G[mob].sum(1)/G[core].sum(1),[0,50,100])]
top=G[oth].mean().sort_values(ascending=False)
V['oth_top5_share']=float(top.head(5).sum()/top.sum()); V['oth_top5']=list(top.index[:5])
V['oth_sd_wo_top5']=float(np.log10(G[list(top.index[5:])].sum(1)).std())
# per-gene SD at plant level
LP=np.log10(G+G.replace(0,np.nan).min()/2).groupby(s.WWTPID).mean()
V['per_gene_sd_mob']=float(LP[mob].std().median()); V['per_gene_sd_oth']=float(LP[oth].std().median())
V['per_gene_sd_p']=float(mannwhitneyu(LP[mob].std(),LP[oth].std()).pvalue)
def mpc(cols):
    c=LP[cols].corr().values; iu=np.triu_indices(len(cols),1); return float(np.nanmean(c[iu]))
V['r_mob']=mpc(mob); V['r_oth']=mpc(oth)
# --- country R2 machinery (independent implementation via group means, not projection) ---
pm=s.groupby('WWTPID').agg(country=('CountryRegion_full','first'),city=('City','first')).loc[LP.index]
def R2(df,lab):
    d=df-df.mean(); m=d.groupby(lab.values).transform('mean')
    return 1-(( (d-m)**2 ).sum()/ (d**2).sum())
obs=R2(LP,pm.country)
cc=pm.groupby('city').country.first(); nulls=[]
for _ in range(999):
    lab=pd.Series(rng.permutation(cc.values),index=cc.index)
    nulls.append(R2(LP,pm.city.map(lab)).values)
N=np.array(nulls); exc=obs.values-N.mean(0)
V['mean_R2_obs']=float(obs.mean()); V['mean_R2_null']=float(N.mean())
V['excess_mob']=float(exc[M].mean()); V['excess_oth']=float(exc[~M].mean())
V['D_obs']=float(obs.values[M].mean()-obs.values[~M].mean())
Dn=N[:,M].mean(1)-N[:,~M].mean(1); V['D_p']=float((np.sum(Dn>=V['D_obs'])+1)/1000)
V['D_excess_diff']=float(exc[M].mean()-exc[~M].mean())
from statsmodels.stats.multitest import multipletests
pg=((N>=obs.values).sum(0)+1)/1000
V['n_sig_fdr']=int(multipletests(pg,method='fdr_bh')[0].sum())
# top / bottom genes
ex=pd.Series(exc,index=core).sort_values(ascending=False)
V['top10']=[(k,round(v,3),bool(g.M_prior[k])) for k,v in ex.head(10).items()]
V['top10_n_mobile']=int(sum(g.M_prior[k] for k in ex.head(10).index))
# within-class
cls=g['class'].str.split(';').str[0]
V['within_class']={c:[float(exc[(cls==c).values&M].mean()),float(exc[(cls==c).values&~M].mean())]
                   for c in ['Aminoglycoside','Beta.lactam','MLS','Tetracycline']}
# SS shares in multi-city countries
multi=['Australia','Brazil','China','Germany','South Africa','Taiwan','USA']; k=pm.country.isin(multi)
d=LP[k]-LP[k].mean(); mc=d.groupby(pm.country[k].values).transform('mean'); mci=d.groupby(pm.city[k].values).transform('mean')
tot=(d**2).sum(); cS=(mc**2).sum()/tot; ciS=((mci-mc)**2).sum()/tot
V['SS_country']=[float(cS[M].mean()),float(cS[~M].mean())]; V['SS_city']=[float(ciS[M].mean()),float(ciS[~M].mean())]
# --- decoding, ONE consistent procedure ---
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import LeaveOneGroupOut, cross_val_predict
Xs=LP[k]; ys=pm.country[k].values; gs=pm.city[k].values
def dec(genes):
    clf=make_pipeline(StandardScaler(),LogisticRegression(C=0.5,max_iter=3000,class_weight='balanced'))
    return balanced_accuracy_score(ys,cross_val_predict(clf,Xs[list(genes)].values,ys,groups=gs,cv=LeaveOneGroupOut()))
DR=np.random.default_rng(7)
pa=[dec(DR.choice(mob,30,replace=False)) for _ in range(20)]
oa=[dec(DR.choice(oth,30,replace=False)) for _ in range(20)]
V['decode_mob']=float(np.mean(pa)); V['decode_oth']=float(np.mean(oa)); V['decode_diff']=float(np.mean(pa)-np.mean(oa))
nd=[]
for _ in range(200):
    p=DR.permutation(core); nd.append(dec(p[:30])-dec(p[30:60]))
nd=np.array(nd); V['decode_p']=float((np.sum(nd>=V['decode_diff'])+1)/201)
V['decode_chance']=1/len(multi)
rand=[dec(DR.choice(core,30,replace=False)) for _ in range(60)]
V['decode_random_mean']=float(np.mean(rand))
json.dump(V,open('verify.json','w'),indent=1,default=str)
for kk,vv in V.items(): print(kk,':',vv)
