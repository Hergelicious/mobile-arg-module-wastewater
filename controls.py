import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, json, re
from scipy.stats import spearmanr, mannwhitneyu
import statsmodels.formula.api as smf
rng=np.random.default_rng(11)
x=pd.ExcelFile('/mnt/user-data/uploads/41467_2025_59019_MOESM8_ESM.xlsx')
s=pd.read_excel(x,'Fig. S2b and Fig. S4',index_col=0)
g=pd.read_csv('gene_table.csv',index_col=0); core=list(g.index); G=s[core].astype(float)
mob=list(g.index[g.M_prior]); oth=list(g.index[~g.M_prior]); M=g.M_prior.values
L=np.log10(G+G.replace(0,np.nan).min()/2)
LP=L.groupby(s.WWTPID).mean()
pm=s.groupby('WWTPID').agg(country=('CountryRegion_full','first'),city=('City','first')).loc[LP.index]
g['mean_log']=np.log10(G.mean()); g['sd_log']=LP.std(); g['prev']=(G>0).mean()
C={}

# ---------- 1. layer variance contrast ----------
lsd=lambda genes: float(np.log10(G[genes].sum(axis=1)).std())
C['layer_sd']={'mobile':lsd(mob),'other':lsd(oth)}
rs=np.array([lsd(list(rng.choice(core,45,replace=False))) for _ in range(2000)])
C['layer_sd_random45']={'median':float(np.median(rs)),'q025':float(np.quantile(rs,.025)),'q975':float(np.quantile(rs,.975)),'p_mobile':float((np.sum(rs>=lsd(mob))+1)/2001)}
# abundance-matched other sets (nearest neighbour on mean log abundance, random tie-breaking, 500 draws, with caliper 0.25)
ml=g.mean_log; mm=[]
for _ in range(500):
    avail=list(oth); pick=[]   # ordered list: set order varies between runs
    for m in rng.permutation(mob):
        cand=[o for o in avail if abs(ml[o]-ml[m])<=0.25]
        if not cand: cand=[min(avail,key=lambda o:abs(ml[o]-ml[m]))]
        o=cand[rng.integers(len(cand))]; pick.append(o); avail.remove(o)
    mm.append(lsd(pick))
mm=np.array(mm)
C['layer_sd_matched_other']={'median':float(np.median(mm)),'q025':float(np.quantile(mm,.025)),'q975':float(np.quantile(mm,.975)),'frac_ge_mobile':float(np.mean(mm>=lsd(mob)))}
C['per_gene_sd_median']={'mobile':float(g.sd_log[g.M_prior].median()),'other':float(g.sd_log[~g.M_prior].median()),'MWU_p':float(mannwhitneyu(g.sd_log[g.M_prior],g.sd_log[~g.M_prior]).pvalue)}
# covariance structure: mean pairwise correlation within sets (plant level)
def mpc(genes):
    c=LP[genes].corr().values; iu=np.triu_indices(len(genes),1); return float(np.nanmean(c[iu]))
C['mean_pairwise_r']={'mobile':mpc(mob),'other':mpc(oth)}
nullr=np.array([mpc(list(rng.choice(core,45,replace=False))) for _ in range(500)])
C['mean_pairwise_r_random45']={'median':float(np.median(nullr)),'q975':float(np.quantile(nullr,.975))}
top=G[oth].mean().sort_values(ascending=False)
C['other_top5']={'genes':list(top.index[:5]),'share':float(top.head(5).sum()/top.sum()),'sd_without_top5':lsd(list(top.index[5:]))}

# ---------- permutation machinery ----------
cc=pm.groupby('city').country.first()
def r2(Y,H):
    Yc=Y-Y.mean(0);P=H@np.linalg.pinv(H.T@H)@H.T;return 1-((Yc-P@Yc)**2).sum(0)/(Yc**2).sum(0)
def perms(P,Y,n=999):
    c2=P.groupby('city').country.first(); out=[]
    for _ in range(n):
        lab=pd.Series(rng.permutation(c2.values),index=c2.index)
        out.append(r2(Y,pd.get_dummies(P.city.map(lab)).values.astype(float)))
    return np.array(out)
Y=LP.values; obs=r2(Y,pd.get_dummies(pm.country).values.astype(float)); null=perms(pm,Y)
exc=obs-null.mean(0)
def pD(stat_fn):
    o=stat_fn(obs); nl=np.array([stat_fn(n) for n in null]); return [float(o),float((np.sum(nl>=o)+1)/(len(nl)+1))]

# ---------- 2. abundance-controlled geographic contrast ----------
g['excess']=exc; g['mob']=g.M_prior.astype(int)
fit=smf.ols('excess ~ mob + mean_log + prev + sd_log',g).fit()
C['regression_excess']={'coef_mob':float(fit.params['mob']),'p_mob':float(fit.pvalues['mob']),'ci':list(map(float,fit.conf_int().loc['mob'])),
   'coef_mean_log':float(fit.params['mean_log']),'p_mean_log':float(fit.pvalues['mean_log']),'coef_sd':float(fit.params['sd_log']),'p_sd':float(fit.pvalues['sd_log']),'r2':float(fit.rsquared)}
# matched-pair D (fixed greedy nearest match, without replacement), permutation tested
avail=list(oth); pairs=[]
for m in g.loc[mob].sort_values('mean_log',ascending=False).index:
    d=(g.loc[avail,'mean_log']-g.loc[m,'mean_log']).abs(); j=d.idxmin(); pairs.append((m,j)); avail.remove(j)
im=[core.index(a) for a,_ in pairs]; io=[core.index(b) for _,b in pairs]
C['D_matched']=pD(lambda r: r[im].mean()-r[io].mean())
C['matched_abundance_diff_median']=float(np.median([abs(ml[a]-ml[b]) for a,b in pairs]))
# abundance-quintile stratified D
q=pd.qcut(g.mean_log,5,labels=False).values
def Dstrat(r):
    v=[r[(q==k)&M].mean()-r[(q==k)&~M].mean() for k in range(5) if ((q==k)&M).sum()>0 and ((q==k)&~M).sum()>0]; return np.mean(v)
C['D_abundance_stratified']=pD(Dstrat)
# within-drug-class D
cls=g['class'].str.split(';').str[0].values
okc=[c for c in sorted(set(cls)) if ((cls==c)&M).sum()>=1 and ((cls==c)&~M).sum()>=1]
def Dcls(r): return np.mean([r[(cls==c)&M].mean()-r[(cls==c)&~M].mean() for c in okc])
C['D_within_class']=pD(Dcls); C['within_class_classes']=sorted(okc)
C['within_class_detail']={c:[float(exc[(cls==c)&M].mean()),float(exc[(cls==c)&~M].mean()),int(((cls==c)&M).sum()),int(((cls==c)&~M).sum())] for c in okc}

# ---------- 3. CLR composition ----------
lc=np.log(G+G.replace(0,np.nan).min()/2); clr=lc.sub(lc.mean(axis=1),axis=0).groupby(s.WWTPID).mean().loc[pm.index].values
oc=r2(clr,pd.get_dummies(pm.country).values.astype(float)); nc=perms(pm,clr,499)
d0=oc[M].mean()-oc[~M].mean(); nd=np.array([n[M].mean()-n[~M].mean() for n in nc])
C['D_CLR']=[float(d0),float((np.sum(nd>=d0)+1)/500)]

# ---------- 4. leave-one-family-out ----------
fam=lambda gname: re.split(r'[\(\-0-9]',gname)[0].lower() or gname
fams=pd.Series({m:fam(m) for m in mob})
LOFO={}
for f in fams.unique():
    keep=np.array([not (c in mob and fams.get(c)==f) for c in core])
    r=obs[keep]; mk=M[keep]; LOFO[f]=float(r[mk].mean()-r[~mk].mean())
C['LOFO_min_max']=[min(LOFO.values()),max(LOFO.values())]; C['LOFO']=LOFO

# ---------- 5. country jackknife and cluster bootstrap CI for D ----------
JK={}
for c in pm.country.unique():
    k=(pm.country!=c).values; o=r2(Y[k],pd.get_dummies(pm.country[k]).values.astype(float)); JK[c]=float(o[M].mean()-o[~M].mean())
C['jackknife_min_max']=[min(JK.values()),max(JK.values())]; C['jackknife']=JK
bs=[]
cities_by=pm.groupby('country').city.unique()
for _ in range(500):
    rows=[]
    for ctry,cl in cities_by.items():
        pick=rng.choice(cl,len(cl),replace=True)
        for i,ci in enumerate(pick):
            idx=np.where(pm.city.values==ci)[0]; rows.append((idx,ctry))
    ii=np.concatenate([r[0] for r in rows]); lab=np.concatenate([[r[1]]*len(r[0]) for r in rows])
    o=r2(Y[ii],pd.get_dummies(pd.Series(lab)).values.astype(float)); bs.append(o[M].mean()-o[~M].mean())
C['D_bootstrap_CI']=[float(np.quantile(bs,.025)),float(np.quantile(bs,.975))]

# ---------- 6. MAG reconciliation ----------
mg=pd.read_excel(x,'Fig. 3c and 3d')
mg['cls']=mg.drug_class.str.replace('-','.')
mag_share=mg.cls.value_counts(normalize=True)
b=pd.read_excel(x,'Fig. 1b and Fig. S2a',index_col=0)
clscols=['Aminoglycoside','Beta.lactam','Glycopeptide','MLS','Sulfonamide','Tetracycline','Phenicol','Nitroimidazole','Rifamycin']
meta_share=(b[clscols].mean()/b[clscols].mean().sum())
rep=pd.DataFrame({'metagenome_share':meta_share,'MAG_share':mag_share.reindex(clscols).fillna(0)})
rep['MAG_over_metagenome']=rep.MAG_share/rep.metagenome_share
C['MAG_representation']=rep.round(4).to_dict('index')
# genes sharing each MAG class x mechanism cell: how many mobile vs other
cellmix={}
for gname in core:
    key=(g.loc[gname,'class'].split(';')[0],g.loc[gname,'mechanism'].split(';')[0])
    cellmix.setdefault(f"{key[0]}|{key[1]}",[0,0])[0 if gname in mob else 1]+=1
C['MAG_cell_gene_mix']=cellmix
C['MAG_phylum_top']=mg.group.value_counts(normalize=True).head(3).round(3).to_dict()
# fraction of mobile-annotated core-gene copies that fall in cells mixing both types
json.dump(C,open('controls.json','w'),indent=1)
for k,v in C.items():
    if k not in ('LOFO','jackknife'): print(k,':',v)
