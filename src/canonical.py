"""Single source of truth for every number in the manuscript. 9999 permutations."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, json
from scipy.stats import spearmanr, mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
rng=np.random.default_rng(20260919)
x=pd.ExcelFile('/mnt/user-data/uploads/41467_2025_59019_MOESM8_ESM.xlsx')
s=pd.read_excel(x,'Fig. S2b and Fig. S4',index_col=0)
f4=pd.read_excel(x,'Fig. 4c',index_col=0).loc[s.index]
g=pd.read_csv('gene_table2.csv',index_col=0); core=list(g.index); M=g.M_prior.values
mob=list(g.index[M]); oth=list(g.index[~M])
G=s[core].astype(float)
LP=np.log10(G+G.replace(0,np.nan).min()/2).groupby(s.WWTPID).mean()
pm=s.groupby('WWTPID').agg(country=('CountryRegion_full','first'),city=('City','first'),GC=('GC','mean')).loc[LP.index]
Y=LP.values; K={}
NP=9999
def r2(Yv,H):
    Yc=Yv-Yv.mean(0); P=H@np.linalg.pinv(H.T@H)@H.T
    return 1-((Yc-P@Yc)**2).sum(0)/(Yc**2).sum(0)
def dum(lab): return pd.get_dummies(lab).values.astype(float)
cc=pm.groupby('city').country.first()
obs=r2(Y,dum(pm.country))
N=np.empty((NP,len(core)))
for i in range(NP):
    lab=pd.Series(rng.permutation(cc.values),index=cc.index); N[i]=r2(Y,dum(pm.city.map(lab)))
exc=obs-N.mean(0)
Dn=N[:,M].mean(1)-N[:,~M].mean(1); Dobs=obs[M].mean()-obs[~M].mean()
K['n_perm']=NP
K['mean_R2_obs']=float(obs.mean()); K['mean_R2_null']=float(N.mean())
K['excess']=[float(exc[M].mean()),float(exc[~M].mean())]; K['D_excess_diff']=float(exc[M].mean()-exc[~M].mean())
K['D']=[float(Dobs),float((np.sum(Dn>=Dobs)+1)/(NP+1))]
pg=((N>=obs).sum(0)+1)/(NP+1); K['n_sig_fdr']=int(multipletests(pg,method='fdr_bh')[0].sum())
ex=pd.Series(exc,index=core).sort_values(ascending=False)
K['top10']=[[k,round(float(v),3),bool(g.M_prior[k])] for k,v in ex.head(10).items()]
K['top10_n_mobile']=int(sum(g.M_prior[k] for k in ex.head(10).index))
K['bottom5']=[[k,round(float(v),3)] for k,v in ex.tail(5).items()]
# across-gene tests, permutation-based
cpl=g.M_coupling.values
ro=spearmanr(cpl,exc)[0]; rn=np.array([spearmanr(cpl,N[i]-N.mean(0))[0] for i in range(2000)])
K['rho_coupling']=[float(ro),float((np.sum(rn>=ro)+1)/2001)]
X=np.c_[np.ones(len(core)),M.astype(float),g.mean_log,g.prevalence,g.sd_log]
bfun=lambda v: np.linalg.lstsq(X,v,rcond=None)[0][1]
bo=bfun(exc); bn=np.array([bfun(N[i]-N.mean(0)) for i in range(2000)])
K['beta_mob']=[float(bo),float((np.sum(bn>=bo)+1)/2001)]
K['coupling_median']=[float(np.median(cpl[M])),float(np.median(cpl[~M]))]
# D variants (999 perms each, adequate for P~0.001 resolution)
def Dperm(P2,Y2,mask,n=999,stat=None):
    H=dum(P2.country); o=r2(Y2,H); c2=P2.groupby('city').country.first()
    st=stat if stat else (lambda r: r[mask].mean()-r[~mask].mean())
    d0=st(o); nl=np.empty(n)
    for i in range(n):
        lab=pd.Series(rng.permutation(c2.values),index=c2.index); nl[i]=st(r2(Y2,dum(P2.city.map(lab))))
    return [float(d0),float((np.sum(nl>=d0)+1)/(n+1))]
allrows=pd.Series(True,index=pm.index)
ml=g.mean_log; avail=list(oth); pairs=[]
for m in g.loc[mob].sort_values('mean_log',ascending=False).index:
    d=(ml[avail]-ml[m]).abs(); j=d.idxmin(); pairs.append((m,j)); avail.remove(j)
im=[core.index(a) for a,_ in pairs]; io=[core.index(b) for _,b in pairs]
K['matched_absdiff_median']=float(np.median([abs(ml[a]-ml[b]) for a,b in pairs]))
K['D_matched']=Dperm(pm,Y,M,stat=lambda r: r[im].mean()-r[io].mean())
q=pd.qcut(g.mean_log,5,labels=False).values
K['D_strat']=Dperm(pm,Y,M,stat=lambda r: np.mean([r[(q==kk)&M].mean()-r[(q==kk)&~M].mean() for kk in range(5)]))
cls=g['class'].str.split(';').str[0].values; okc=['Aminoglycoside','Beta.lactam','MLS','Tetracycline']
K['D_class']=Dperm(pm,Y,M,stat=lambda r: np.mean([r[(cls==c)&M].mean()-r[(cls==c)&~M].mean() for c in okc]))
K['within_class']={c:[float(exc[(cls==c)&M].mean()),float(exc[(cls==c)&~M].mean()),int(((cls==c)&M).sum()),int(((cls==c)&~M).sum())] for c in okc}
lc=np.log(G+G.replace(0,np.nan).min()/2); CLR=lc.sub(lc.mean(1),axis=0).groupby(s.WWTPID).mean().loc[pm.index].values
K['D_CLR']=Dperm(pm,CLR,M)
def resid(y,xv):
    Xm=np.c_[np.ones(len(xv)),xv]; return y-Xm@np.linalg.lstsq(Xm,y,rcond=None)[0]
Ygc=np.column_stack([resid(Y[:,i],pm.GC.values) for i in range(Y.shape[1])])
K['D_GC']=Dperm(pm,Ygc,M)
hi=cpl>=np.quantile(cpl,2/3); lo=cpl<=np.quantile(cpl,1/3)
K['D_tertile']=Dperm(pm,Y,M,stat=lambda r: r[hi].mean()-r[lo].mean())
for nm,exc_c in [('noUSA',['USA']),('noCHN',['China']),('noBoth',['USA','China'])]:
    k=~pm.country.isin(exc_c); K['D_'+nm]=Dperm(pm[k],Y[k.values],M)
prev80=(G>0).mean()>0.8; m80=M[prev80.values]
K['D_prev80']=[float(r2(Y[:,prev80.values],dum(pm.country))[m80].mean()-r2(Y[:,prev80.values],dum(pm.country))[~m80].mean()),int(prev80.sum())]
# jackknife / LOFO
import re as _re
K['jackknife']={}
for c in pm.country.unique():
    k=(pm.country!=c).values; o=r2(Y[k],dum(pm.country[k])); K['jackknife'][c]=float(o[M].mean()-o[~M].mean())
fam=lambda n: _re.split(r'[\(\-0-9]',n)[0].lower() or n
fams={m:fam(m) for m in mob}; K['LOFO']={}
for f in sorted(set(fams.values())):
    keep=np.array([not (c in mob and fams[c]==f) for c in core])
    o=obs[keep]; mk=M[keep]; K['LOFO'][f]=float(o[mk].mean()-o[~mk].mean())
K['jk_range']=[min(K['jackknife'].values()),max(K['jackknife'].values())]
K['LOFO_range']=[min(K['LOFO'].values()),max(K['LOFO'].values())]
# bootstrap CI
bs=np.empty(1000); cb=pm.groupby('country').city.unique()
for i in range(1000):
    idx=[]; lab=[]
    for ctry,cl in cb.items():
        for ci in rng.choice(cl,len(cl),replace=True):
            w=np.where(pm.city.values==ci)[0]; idx+=list(w); lab+=[ctry]*len(w)
    o=r2(Y[idx],dum(pd.Series(lab))); bs[i]=o[M].mean()-o[~M].mean()
K['D_boot_CI']=[float(np.quantile(bs,.025)),float(np.quantile(bs,.975))]
# SS shares
multi=['Australia','Brazil','China','Germany','South Africa','Taiwan','USA']; k=pm.country.isin(multi)
def SS(lab_c,lab_ci):
    d=LP[k]-LP[k].mean(); mc=d.groupby(lab_c).transform('mean'); mci=d.groupby(lab_ci).transform('mean'); t=(d**2).sum()
    return ((mc**2).sum()/t).values, (((mci-mc)**2).sum()/t).values
cS,ciS=SS(pm.country[k].values,pm.city[k].values)
st0=(cS[M].mean()-cS[~M].mean())-(ciS[M].mean()-ciS[~M].mean())
c2=pm[k].groupby('city').country.first(); nl=np.empty(999)
for i in range(999):
    lab=pd.Series(rng.permutation(c2.values),index=c2.index)
    a,b=SS(pm.city[k].map(lab).values,pm.city[k].values); nl[i]=(a[M].mean()-a[~M].mean())-(b[M].mean()-b[~M].mean())
K['SS_country']=[float(cS[M].mean()),float(cS[~M].mean())]; K['SS_city']=[float(ciS[M].mean()),float(ciS[~M].mean())]
K['SS_contrast_p']=float((np.sum(nl>=st0)+1)/1000)
# layers
lmob=np.log10(G[mob].sum(1)); loth=np.log10(G[oth].sum(1)); share=G[mob].sum(1)/G[core].sum(1)
K['layer_sd']=[float(lmob.std()),float(loth.std())]
K['per_gene_sd']=[float(LP[mob].std().median()),float(LP[oth].std().median()),float(mannwhitneyu(LP[mob].std(),LP[oth].std()).pvalue)]
def mpc(cols):
    c=LP[cols].corr().values; iu=np.triu_indices(len(cols),1); return float(np.nanmean(c[iu]))
K['r_within']=[mpc(mob),mpc(oth)]
rr=np.array([mpc(list(rng.choice(core,45,replace=False))) for _ in range(500)])
K['r_random45']=[float(np.median(rr)),float(np.quantile(rr,.975))]
lsd=lambda cols: float(np.log10(G[cols].sum(1)).std())
rs=np.array([lsd(list(rng.choice(core,45,replace=False))) for _ in range(2000)])
K['sd_random45']=[float(np.median(rs)),float((np.sum(rs>=K['layer_sd'][0])+1)/2001)]
mm=[]
for _ in range(500):
    av=list(oth); pick=[]   # ordered list: set order varies between runs
    for m in rng.permutation(mob):
        cand=[o for o in av if abs(ml[o]-ml[m])<=0.25] or [min(av,key=lambda o:abs(ml[o]-ml[m]))]
        o=cand[rng.integers(len(cand))]; pick.append(o); av.remove(o)
    mm.append(lsd(pick))
mm=np.array(mm)
K['sd_matched']=[float(np.median(mm)),float(np.quantile(mm,.025)),float(np.quantile(mm,.975)),float(np.mean(mm>=K['layer_sd'][0]))]
K['sd_ratio_matched']=float(K['layer_sd'][0]/np.median(mm))
top=G[oth].mean().sort_values(ascending=False)
K['oth_top5_share']=float(top.head(5).sum()/top.sum()); K['oth_sd_wo_top5']=lsd(list(top.index[5:]))
K['share_oth_median']=float((G[oth].sum(1)/G[core].sum(1)).median())
K['mobile_share']=[float(share.min()),float(share.median()),float(share.max())]
# REML
L=pd.DataFrame({'mob':lmob,'oth':loth,'country':s.CountryRegion_full,'city':s.City,'plant':s.WWTPID})
K['REML']={}
for yv in ['mob','oth']:
    m=smf.mixedlm(f"{yv} ~ 1",L,groups=L.country,re_formula="1",vc_formula={"city":"0+C(city)","plant":"0+C(plant)"}).fit(reml=True)
    d=dict(zip(m.model.exog_vc.names,[float(v) for v in m.vcomp])); d['country']=float(m.cov_re.iloc[0,0]); d['resid']=float(m.scale)
    tot=sum(d.values()); K['REML'][yv]={'abs':d,'frac':{kk:vv/tot for kk,vv in d.items()}}
K['REML_ratio']=K['REML']['mob']['abs']['country']/K['REML']['oth']['abs']['country']
json.dump(K,open('canonical.json','w'),indent=1)
print('D',K['D'],'excess',K['excess'],'nsig',K['n_sig_fdr'],'CI',K['D_boot_CI'])
print('rho',K['rho_coupling'],'beta',K['beta_mob'],'SSp',K['SS_contrast_p'])
print('top10 mobile',K['top10_n_mobile'],K['top10'][:7])
