"""
Two analyses using Zhu et al. Supplementary Data 1 (sequencing output) and 3 (MGE abundances).
DESIGN FIXED BEFORE RUNNING.

A. Direct integron coupling (the mobility test the paper lacked)
   Marker   : intI1 (class 1 integron integrase; Gillings et al. 2015) = sum of the 7 intI1 entries per sample.
              Secondary marker: qacEdelta1 (4 entries), part of the class 1 integron 3' conserved segment.
   Unit     : plant means of log10(abundance + half the smallest non-zero value).
   Statistic: for each of the 114 core genes, the WITHIN-COUNTRY partial Spearman correlation with intI1,
              controlling for total ARG load (same construction as the MGE-coupling index).
              D_int = mean correlation (mobile-associated) - mean correlation (other genes).
   Null     : intI1 values permuted among plants WITHIN each country (1,999 permutations). This keeps
              every country's intI1 distribution and all gene-gene correlations intact; only the
              plant-level pairing of intI1 with the resistome is broken.
   Also     : (i) curated cassette-type genes vs other genes; (ii) non-annotated members of the
              label-free module vs non-mobile genes outside it (tests the module beyond the labels);
              (iii) sensitivity excluding sul1 (physically part of the integron 3' segment).
B. Sequencing depth
   Covariate: plant mean log10(HQ bases). D recomputed on gene abundances residualised on depth,
              999 city-block permutations (same as the primary analysis).
"""
import warnings; warnings.filterwarnings("ignore")
import json, re, numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr
rng = np.random.default_rng(1729)
U = '/mnt/user-data/uploads/'
s = pd.read_excel(U + '41467_2025_59019_MOESM8_ESM.xlsx', 'Fig. S2b and Fig. S4', index_col=0)
q = pd.read_excel(U + '41467_2025_59019_MOESM3_ESM.xlsx', header=2).set_index('Sample').loc[s.index]
m = pd.read_excel(U + '41467_2025_59019_MOESM5_ESM.xlsx', header=2); m = m.rename(columns={m.columns[0]: 'id'})
g = pd.read_csv('gene_table_final.csv', index_col=0); core = list(g.index); M = g.M_prior.values.astype(bool)
mm = pd.read_csv('module_membership.csv').set_index('gene').loc[core]
inmod = mm.in_module.values.astype(bool)
K = json.load(open('canonical.json'))

G = s[core].astype(float)
LP = np.log10(G + G.replace(0, np.nan).min() / 2).groupby(s.WWTPID).mean()
pm = s.groupby('WWTPID').agg(country=('CountryRegion_full', 'first'), city=('City', 'first')).loc[LP.index]
tot = np.log10(s.ARG_per_cell).groupby(s.WWTPID).mean().loc[LP.index]

# sanity: identical matrix to the canonical analysis
def r2(Y, H):
    Yc = Y - Y.mean(0); P = H @ np.linalg.pinv(H.T @ H) @ H.T
    return 1 - ((Yc - P @ Yc) ** 2).sum(0) / (Yc ** 2).sum(0)
dum = lambda l: pd.get_dummies(l).values.astype(float)
o = r2(LP.values, dum(pm.country)); assert abs((o[M].mean() - o[~M].mean()) - K['D'][0]) < 1e-12

def marker(name):
    rows = m[m.MGE_Gene_name.astype(str) == name]
    v = rows[list(s.index)].sum(axis=0).astype(float)
    pc = v[v > 0].min() / 2
    return np.log10(v + pc).groupby(s.WWTPID).mean().loc[LP.index], int(len(rows)), float((v > 0).mean())
intI1, n_int, prev_int = marker('intI1')
qac, n_qac, prev_qac = marker('qacEdelta')
OUT = {'intI1_entries': n_int, 'intI1_prevalence': prev_int, 'qac_entries': n_qac, 'qac_prevalence': prev_qac}

ctry = pm.country.values
def demean(v): return v - pd.Series(v).groupby(ctry).transform('mean').values
def resid(y, x):
    X = np.c_[np.ones(len(x)), x]; return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
rt = rankdata(demean(tot.values))
Yr = np.column_stack([resid(rankdata(demean(LP.values[:, i])), rt) for i in range(len(core))])
Yr = (Yr - Yr.mean(0)) / Yr.std(0)
def corr_with(marker_vals):
    mr = resid(rankdata(demean(marker_vals)), rt); mr = (mr - mr.mean()) / mr.std()
    return (Yr * mr[:, None]).mean(0)

groups = {}
for c in np.unique(ctry): groups[c] = np.where(ctry == c)[0]
def within_perm(v):
    w = v.copy()
    for idx in groups.values():
        if len(idx) > 1: w[idx] = v[rng.permutation(idx)]
    return w

cassette_re = re.compile(r"^(sul1|aadA|ANT\(3''\)|catB|cmlA|OXA-(1|2|5|9|10|17|46)$|GES|VEB|AAC\(6'\)-Ib|EreA)")
CAS = np.array([bool(cassette_re.match(x)) for x in core])
OUT['n_cassette_genes'] = int(CAS.sum())
nonmob_in = inmod & ~M; nonmob_out = ~inmod & ~M
nos1 = np.array([x != 'sul1' for x in core])
stats = {
  'mobile_vs_other':      lambda r: r[M].mean() - r[~M].mean(),
  'cassette_vs_other':    lambda r: r[CAS].mean() - r[~M].mean(),
  'module_nonannot_vs_outside': lambda r: r[nonmob_in].mean() - r[nonmob_out].mean(),
  'mobile_vs_other_no_sul1':    lambda r: r[M & nos1].mean() - r[~M].mean(),
}
for mk, vals in [('intI1', intI1.values), ('qacEdelta1', qac.values)]:
    r = corr_with(vals)
    NP = 1999; null = {k: np.empty(NP) for k in stats}
    for i in range(NP):
        rp = corr_with(within_perm(vals))
        for k, f in stats.items(): null[k][i] = f(rp)
    res = {}
    for k, f in stats.items():
        ob = f(r); res[k] = [float(ob), float((np.sum(null[k] >= ob) + 1) / (NP + 1))]
    res['mean_r_mobile'] = float(r[M].mean()); res['mean_r_other'] = float(r[~M].mean())
    res['mean_r_cassette'] = float(r[CAS].mean())
    res['mean_r_module_nonannot'] = float(r[nonmob_in].mean()); res['mean_r_outside_nonmobile'] = float(r[nonmob_out].mean())
    res['top10'] = [[core[i], round(float(r[i]), 3), bool(M[i])] for i in np.argsort(-r)[:10]]
    OUT[mk] = res
    if mk == 'intI1': np.save('intI1_corr.npy', r); np.save('intI1_null.npy', null['mobile_vs_other'])

# ---- B. sequencing depth ----
depth = np.log10(q['HQ bases'].astype(float)).groupby(s.WWTPID).mean().loc[LP.index]
share = np.log10((G[g.index[g.M_prior]].sum(1) / G[core].sum(1))).groupby(s.WWTPID).mean().loc[LP.index]
OUT['depth'] = {'HQ_Gb_median': float(np.median(q['HQ bases']) / 1e9), 'HQ_Gb_range': [float(q['HQ bases'].min() / 1e9), float(q['HQ bases'].max() / 1e9)],
                'country_R2': float(r2(depth.values[:, None], dum(pm.country))[0]),
                'rho_mobile_share': float(spearmanr(depth, share)[0])}
Yd = np.column_stack([resid(LP.values[:, i], depth.values) for i in range(len(core))])
cc = pm.groupby('city').country.first()
od = r2(Yd, dum(pm.country)); d0 = od[M].mean() - od[~M].mean(); nl = np.empty(999)
for i in range(999):
    lab = pd.Series(rng.permutation(cc.values), index=cc.index); oo = r2(Yd, dum(pm.city.map(lab))); nl[i] = oo[M].mean() - oo[~M].mean()
OUT['depth']['D_adjusted'] = [float(d0), float((np.sum(nl >= d0) + 1) / 1000)]
json.dump(OUT, open('integron_depth.json', 'w'), indent=1)
print(json.dumps(OUT, indent=1))
