"""
Unsupervised module discovery. DESIGN FIXED BEFORE RUNNING (no tuning to results):

Input   : plant-level log10 per-cell abundance of the 114 core ARGs (identical matrix to canonical.py).
Primary : gene-gene Pearson correlations computed on COUNTRY-DEMEANED data (within-country co-variation).
          -> discovery is blind to (a) the mobility labels and (b) between-country differences,
             so the later between-country test is not circular.
Cluster : average-linkage hierarchical clustering on distance 1 - r.
k       : the k in 2..12 that maximises mean silhouette (precomputed distance). Chosen automatically.
Enrich  : hypergeometric test of mobile-associated over-representation for every cluster with >=5 genes,
          Benjamini-Hochberg across clusters.
Geo test: D (module genes vs all other genes), observed country R2, 999 city-block permutations.
Stable  : 200 plant bootstraps -> recluster with the same k -> co-clustering frequency, ARI vs original.
Sensitiv: complete linkage; k-1 and k+1; raw (not demeaned) correlations.
All outputs -> module.json. Nothing is hand-copied.
"""
import warnings; warnings.filterwarnings("ignore")
import json, numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import hypergeom
from sklearn.metrics import silhouette_score, adjusted_rand_score
from statsmodels.stats.multitest import multipletests

rng = np.random.default_rng(777)
x = pd.ExcelFile('/mnt/user-data/uploads/41467_2025_59019_MOESM8_ESM.xlsx')
s = pd.read_excel(x, 'Fig. S2b and Fig. S4', index_col=0)
g = pd.read_csv('gene_table_final.csv', index_col=0)
core = list(g.index); M = g.M_prior.values.astype(bool)
G = s[core].astype(float)
LP = np.log10(G + G.replace(0, np.nan).min() / 2).groupby(s.WWTPID).mean()
pm = s.groupby('WWTPID').agg(country=('CountryRegion_full', 'first'), city=('City', 'first')).loc[LP.index]

# sanity: same matrix as the canonical analysis
K = json.load(open('canonical.json'))
def r2(Y, H):
    Yc = Y - Y.mean(0); P = H @ np.linalg.pinv(H.T @ H) @ H.T
    return 1 - ((Yc - P @ Yc) ** 2).sum(0) / (Yc ** 2).sum(0)
dum = lambda l: pd.get_dummies(l).values.astype(float)
obs = r2(LP.values, dum(pm.country))
assert abs((obs[M].mean() - obs[~M].mean()) - K['D'][0]) < 1e-12, "matrix mismatch with canonical"

def corr_dist(df, demean=True):
    d = df - df.groupby(pm.loc[df.index, 'country'].values).transform('mean') if demean else df
    c = np.corrcoef(d.values.T)
    np.fill_diagonal(c, 1.0)
    return c, np.clip(1 - c, 0, 2)

def cluster(dist, k, method='average'):
    Z = linkage(squareform(dist, checks=False), method=method)
    return fcluster(Z, k, criterion='maxclust')

def choose_k(dist, method='average'):
    sil = {}
    for k in range(2, 13):
        lab = cluster(dist, k, method)
        if len(set(lab)) > 1: sil[k] = silhouette_score(dist, lab, metric='precomputed')
    return max(sil, key=sil.get), sil

def enrichment(lab):
    rows = []
    N, Kmob = len(core), int(M.sum())
    for c in sorted(set(lab)):
        idx = lab == c; n = int(idx.sum()); m = int((idx & M).sum())
        p = hypergeom.sf(m - 1, N, Kmob, n) if n >= 5 else np.nan
        rows.append(dict(cluster=int(c), size=n, n_mobile=m, frac_mobile=m / n, p=p))
    df = pd.DataFrame(rows)
    ok = df.p.notna()
    df.loc[ok, 'q'] = multipletests(df.loc[ok, 'p'], method='fdr_bh')[1]
    return df

OUT = {}
C, Dm = corr_dist(LP, demean=True)
k, sil = choose_k(Dm)
lab = cluster(Dm, k)
enr = enrichment(lab)
OUT['k'] = int(k); OUT['silhouette'] = {int(a): float(b) for a, b in sil.items()}
OUT['clusters'] = enr.to_dict('records')
best = enr.loc[enr.q.idxmin()] if enr.q.notna().any() else None
mod_id = int(best.cluster)
mod = lab == mod_id
OUT['module'] = dict(cluster=mod_id, size=int(mod.sum()), n_mobile=int((mod & M).sum()),
                     p=float(best.p), q=float(best.q),
                     precision=float((mod & M).sum() / mod.sum()),
                     recall=float((mod & M).sum() / M.sum()),
                     genes=[core[i] for i in np.where(mod)[0]],
                     non_mobile_members=[core[i] for i in np.where(mod & ~M)[0]],
                     mobile_outside=[core[i] for i in np.where(~mod & M)[0]])
# within-module correlation (within-country) vs rest
iu = lambda mask: C[np.ix_(mask, mask)][np.triu_indices(mask.sum(), 1)]
OUT['module_mean_within_r'] = float(iu(mod).mean()); OUT['rest_mean_within_r'] = float(iu(~mod).mean())

# ---- geographic test of the DISCOVERED module (not circular: discovery used demeaned data) ----
cc = pm.groupby('city').country.first()
Dobs = obs[mod].mean() - obs[~mod].mean()
nl = np.empty(999)
for i in range(999):
    labp = pd.Series(rng.permutation(cc.values), index=cc.index)
    o = r2(LP.values, dum(pm.city.map(labp))); nl[i] = o[mod].mean() - o[~mod].mean()
OUT['D_module'] = [float(Dobs), float((np.sum(nl >= Dobs) + 1) / 1000)]
# same test restricted to NON-mobile module members vs non-mobile rest (does the module add beyond the labels?)
nm_in = mod & ~M; nm_out = ~mod & ~M
if nm_in.sum() >= 3:
    d2 = obs[nm_in].mean() - obs[nm_out].mean()
    nl2 = np.empty(999)
    for i in range(999):
        labp = pd.Series(rng.permutation(cc.values), index=cc.index)
        o = r2(LP.values, dum(pm.city.map(labp))); nl2[i] = o[nm_in].mean() - o[nm_out].mean()
    OUT['D_nonmobile_in_vs_out'] = [float(d2), float((np.sum(nl2 >= d2) + 1) / 1000), int(nm_in.sum())]
# module score
score = np.log10(G[[core[i] for i in np.where(mod)[0]]].sum(axis=1))
tot = np.log10(s.ARG_per_cell)
def r2s(v, grp): return float(1 - ((v - v.groupby(grp).transform('mean')) ** 2).sum() / ((v - v.mean()) ** 2).sum())
OUT['module_score'] = dict(sd=float(score.std()), R2_country=r2s(score, s.CountryRegion_full),
                           R2_country_total=r2s(tot, s.CountryRegion_full))

# ---- stability: plant bootstrap ----
ari, co = [], np.zeros((len(core), len(core)))
plants = LP.index.values
for b in range(200):
    pick = rng.choice(plants, len(plants), replace=True)
    df = LP.loc[pick]; df.index = range(len(pick))
    ctry = pm.loc[pick, 'country'].values
    d = df - df.groupby(ctry).transform('mean')
    cb = np.corrcoef(d.values.T); np.fill_diagonal(cb, 1)
    lb = cluster(np.clip(1 - cb, 0, 2), k)
    ari.append(adjusted_rand_score(lab, lb))
    co += (lb[:, None] == lb[None, :])
co /= 200
OUT['bootstrap_ARI'] = [float(np.median(ari)), float(np.quantile(ari, .05)), float(np.quantile(ari, .95))]
OUT['module_cocluster'] = float(co[np.ix_(mod, mod)][np.triu_indices(mod.sum(), 1)].mean())
OUT['module_member_stability'] = {core[i]: float(co[i, mod].mean()) for i in np.where(mod)[0]}

# ---- sensitivity ----
sens = {}
for name, dist, meth, kk in [('complete_linkage', Dm, 'complete', None),
                             ('k_minus_1', Dm, 'average', k - 1), ('k_plus_1', Dm, 'average', k + 1),
                             ('raw_correlations', corr_dist(LP, demean=False)[1], 'average', None)]:
    kk = kk if kk else choose_k(dist, meth)[0]
    lb = cluster(dist, kk, meth); e = enrichment(lb)
    b2 = e.loc[e.q.idxmin()]
    m2 = lb == int(b2.cluster)
    sens[name] = dict(k=int(kk), size=int(b2['size']), n_mobile=int(b2.n_mobile), q=float(b2.q),
                      jaccard_with_primary=float((m2 & mod).sum() / (m2 | mod).sum()))
OUT['sensitivity'] = sens

pd.DataFrame({'gene': core, 'cluster': lab, 'in_module': mod, 'mobile_associated': M}).to_csv('module_membership.csv', index=False)
np.save('module_Dnull.npy', nl)
def _clean(o):
    if isinstance(o, dict): return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, list): return [_clean(v) for v in o]
    if isinstance(o, float) and o != o: return None   # NaN -> null (strict JSON, readable by the Node build)
    return o
json.dump(_clean(OUT), open('module.json', 'w'), indent=1, allow_nan=False)
print(json.dumps({k2: v for k2, v in OUT.items() if k2 not in ('module_member_stability', 'silhouette')}, indent=1, default=str))
