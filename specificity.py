"""Is the module special, or would any comparable gene set do?
SIZE CONSTRAINT: the module contains 77 of the 114 core genes, so equal-sized alternatives cannot avoid overlapping it.
Comparisons therefore use 30-gene subsets: 30 random module genes against 30-gene sets drawn from the non-module genes,
matched on abundance, prevalence or drug class, and against 30 random core genes.
Compares the 77-gene module with six alternative sets of the same size on four criteria.
CIRCULARITY IS EXPLICIT: the module was defined by within-country co-variation, so criterion (i) is not an
independent test for it; criteria (ii)-(iv) are. Random and matched sets are drawn 200 times; medians reported."""
import warnings; warnings.filterwarnings("ignore")
import json, numpy as np, pandas as pd
rng=np.random.default_rng(808)
U='/mnt/user-data/uploads/'
g=pd.read_csv('gene_table_final.csv',index_col=0); core=list(g.index)
mm=pd.read_csv('module_membership.csv').set_index('gene').loc[core]
MOD=[x for x in core if bool(mm.loc[x,'in_module'])]; MOB=[x for x in core if bool(g.M_prior[x])]
z=pd.read_excel(U+'41467_2025_59019_MOESM8_ESM.xlsx','Fig. S2b and Fig. S4',index_col=0)
G=z[core].astype(float); L=np.log10(G+G.replace(0,np.nan).min()/2)
LP=L.groupby(z.WWTPID).mean(); ctry=z.groupby('WWTPID').CountryRegion_full.first().loc[LP.index]
W=LP-LP.groupby(ctry.values).transform('mean'); C=np.corrcoef(W.values.T); np.fill_diagonal(C,np.nan)
idx={x:i for i,x in enumerate(core)}
R2=g.R2_country; rho_int=pd.read_csv('intI1_rho.csv',index_col=0).iloc[:,0]
def icc1(v,grp):
    d=pd.DataFrame({'v':v,'g':grp}); d=d[d.g.map(d.g.value_counts())>=2]
    k=d.g.value_counts(); n=len(k); N=len(d)
    msb=(k*(d.groupby('g').v.mean()-d.v.mean())**2).sum()/(n-1); msw=((d.v-d.groupby('g').v.transform('mean'))**2).sum()/(N-n)
    kbar=(N-(k**2).sum()/N)/(n-1); return float((msb-msw)/(msb+(kbar-1)*msw))
def profile(genes):
    ii=[idx[x] for x in genes]; sub=C[np.ix_(ii,ii)]
    return {'co_variation_r':float(np.nanmean(sub)),'country_R2':float(R2[genes].mean()),
            'reliability_ICC':float(icc1(L[genes].mean(1),z.WWTPID)),'rho_intI1':float(rho_int[genes].mean())}
ab=g.mean_log; pv=g.prevalence; cls=g['class'].str.split(';').str[0]
NONMOD=[x for x in core if x not in MOD]; K=30
def draw_matched(target,on):
    pick=[]; avail=list(NONMOD)
    for x in target:
        if on=='abundance': cand=[y for y in avail if abs(ab[y]-ab[x])<0.5]
        elif on=='prevalence': cand=[y for y in avail if abs(pv[y]-pv[x])<0.1]
        elif on=='class': cand=[y for y in avail if cls[y]==cls[x]]
        else: cand=list(avail)
        if not cand: cand=list(avail)
        y=cand[rng.integers(len(cand))]; pick.append(y); avail.remove(y)
    return pick
def med(sets):
    return {k:round(float(np.median([profile(s_)[k] for s_ in sets])),3) for k in ['co_variation_r','country_R2','reliability_ICC','rho_intI1']}
mod_draws=[list(rng.choice(MOD,K,replace=False)) for _ in range(200)]
OUT={'set_size':K,
 'module_subset':med(mod_draws),
 'non_module':med([list(rng.choice(NONMOD,K,replace=False)) for _ in range(200)]),
 'random_core':med([list(rng.choice(core,K,replace=False)) for _ in range(200)]),
 'abundance_matched':med([draw_matched(t,'abundance') for t in mod_draws]),
 'prevalence_matched':med([draw_matched(t,'prevalence') for t in mod_draws]),
 'class_matched':med([draw_matched(t,'class') for t in mod_draws]),
 'curated_mobile_45_full':{k:round(v,3) for k,v in profile(MOB).items()},
 'module_77_full':{k:round(v,3) for k,v in profile(MOD).items()},
 'note':'The module was defined by within-country co-variation, so co_variation_r is not an independent criterion for it; country_R2, reliability and the intI1 association are.'}
json.dump(OUT,open('specificity.json','w'),indent=1)
print(pd.DataFrame({k:v for k,v in OUT.items() if isinstance(v,dict)}).T.to_string())
