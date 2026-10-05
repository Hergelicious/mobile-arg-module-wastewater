"""
Two-layer resistome analysis of Zhu et al. (Nat Commun 2025) source data.
Uses ONLY the published source-data workbook (41467_2025_59019_MOESM8_ESM.xlsx).
"""
import re, json, warnings
import numpy as np, pandas as pd
from scipy.optimize import nnls
from scipy.stats import spearmanr, mannwhitneyu
warnings.filterwarnings("ignore")
rng = np.random.default_rng(2026)
F = "/mnt/user-data/uploads/41467_2025_59019_MOESM8_ESM.xlsx"
x = pd.ExcelFile(F)
R = {}  # results

# ---------------- data ----------------
s = pd.read_excel(x, "Fig. S2b and Fig. S4", index_col=0)
b = pd.read_excel(x, "Fig. 1b and Fig. S2a", index_col=0).loc[s.index]
f4 = pd.read_excel(x, "Fig. 4c", index_col=0).loc[s.index]
gcols = list(s.columns[s.columns.get_loc("ARG_per_cell") + 1:])
G = s[gcols].astype(float)
classes = ['Aminoglycoside','Beta.lactam','Chloramphenicol','Diaminopyrimidine','Fluoroquinolone',
           'Glycopeptide','MLS','Multidrug','Nitroimidazole','Penam','Peptide','Phenicol',
           'Rifamycin','Sulfonamide','Tetracycline']
mechs = ['antibiotic.efflux','antibiotic.inactivation','antibiotic.target.alteration',
         'antibiotic.target.protection','antibiotic.target.replacement']

# gene -> class / mechanism recovered from the workbook (class totals = sums of subtypes)
def assign(cols):
    out = {}
    A = b[gcols].values
    for c in cols:
        coef, _ = nnls(A, b[c].values)
        for g, w in zip(gcols, coef):
            if w > 0.5: out.setdefault(g, []).append(c)
    return out
gcls, gmech = assign(classes), assign(mechs)
R["n_genes_total"] = len(gcols)
R["n_genes_class_assigned"] = sum(1 for g in gcols if g in gcls)

meta = s[["Continent","CountryRegion_full","City","WWTPID","GC","GDP.per.capita..dollars","population"]].copy()
meta.columns = ["continent","country","city","plant","GC","gdp","pop"]
meta["total"] = s["ARG_per_cell"]
meta["mge"] = f4["MGE.nearARG.percell"]

R["n_samples"], R["n_plants"], R["n_cities"], R["n_countries"] = len(s), meta.plant.nunique(), meta.city.nunique(), meta.country.nunique()
R["samples_per_country"] = meta.country.value_counts().to_dict()
R["cities_per_country"] = meta.groupby("country").city.nunique().to_dict()

# core gene set: prevalence > 50%
prev = (G > 0).mean()
core = prev[prev > 0.5].index.tolist()
R["n_core_genes"] = len(core)
pseudo = G[core].replace(0, np.nan).min() / 2
LG = np.log10(G[core] + pseudo)

# plant-level aggregation
pm = meta.groupby("plant").agg(country=("country","first"), city=("city","first"), GC=("GC","mean"),
                               gdp=("gdp","mean"), pop=("pop","mean"),
                               total=("total", lambda v: np.log10(v).mean()),
                               mge=("mge", lambda v: np.log10(v + 0.005).mean()))
LP = LG.groupby(meta.plant).mean().loc[pm.index]

# ---------------- mobility indices ----------------
# M_prior: a-priori list of genes classically carried on class 1 integrons / plasmids / transposons
prior_re = re.compile(r'^(sul1|sul2|sul3|aadA|dfrA|cmlA|catB|floR|OXA-(1|2|5|9|10|17|20|46|129|198)$|arr-3|GES|VEB|PER|IMP|VIM|'
                      r'ANT\(3|AAC\(6.\)-Ib|APH\(6\)|APH\(3..\)-Ib|ErmB|mph|msrE|EreA|lnuF|tet\(M\)|tet\(Q\)|tet\(X\)|tet\(C\)|tet\(A\)|QnrS)')
M_prior = pd.Series([bool(prior_re.match(g)) for g in core], index=core)

# M_cpl: data-derived MGE coupling. Within-country (country-demeaned) partial Spearman between
# gene abundance and MGE-adjacent-ARG density, controlling for total ARG load, at plant level.
def demean(df, by):
    return df - df.groupby(by).transform("mean")
Z = demean(LP, pm.country); mz = demean(pm[["mge"]], pm.country)["mge"]; tz = demean(pm[["total"]], pm.country)["total"]
def resid(y, xv):
    X = np.c_[np.ones(len(xv)), xv]; return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
from scipy.stats import rankdata
mr = resid(rankdata(mz), rankdata(tz))
M_cpl = pd.Series({g: spearmanr(resid(rankdata(Z[g]), rankdata(tz)), mr)[0] for g in core})

# M_mag: MAG-derived mobile fraction per class x mechanism (only cells with >=10 ORFs)
mg = pd.read_excel(x, "Fig. 3c and 3d")
mg["cls"] = mg.drug_class.str.replace("-", ".")
mg["mech"] = "antibiotic." + mg.resistance_mechanism.str.replace("antibiotic ", "").str.replace(" ", ".")
cell = mg.groupby(["cls","mech"]).mobile.agg(["mean","size"])
cell = cell[cell["size"] >= 10]
M_mag = pd.Series({g: cell.loc[(gcls[g][0], gmech[g][0]), "mean"]
                   if g in gcls and g in gmech and (gcls[g][0], gmech[g][0]) in cell.index else np.nan for g in core})
R["mag_cells"] = {f"{a}|{m}": [round(v,3), int(n)] for (a,m),(v,n) in cell.iterrows()}
R["mobility_index_agreement"] = {
    "Mcpl_prior_vs_other_median": [float(M_cpl[M_prior].median()), float(M_cpl[~M_prior].median())],
    "Mcpl_prior_MWU_p": float(mannwhitneyu(M_cpl[M_prior], M_cpl[~M_prior]).pvalue),
    "Mcpl_vs_Mmag_rho": list(map(float, spearmanr(M_cpl, M_mag, nan_policy="omit")[:2])),
    "n_genes_with_Mmag": int(M_mag.notna().sum())}

# ---------------- geographic signal ----------------
cities = pm.city.values
city_country = pm.groupby("city").country.first()
Hc = pd.get_dummies(pm.country).values.astype(float)
def r2_country(Y, H):
    Yc = Y - Y.mean(0); P = H @ np.linalg.pinv(H.T @ H) @ H.T
    return 1 - ((Yc - P @ Yc) ** 2).sum(0) / (Yc ** 2).sum(0)
def perm_H():
    # shuffle the city->country assignment (keeps cities intact, preserves nesting)
    lab = pd.Series(rng.permutation(city_country.values), index=city_country.index)
    return pd.get_dummies(pm.city.map(lab)).values.astype(float)
Y = LP.values
obs = r2_country(Y, Hc)
NPERM = 999
null = np.array([r2_country(Y, perm_H()) for _ in range(NPERM)])
excess = pd.Series(obs - null.mean(0), index=core)
pgene = pd.Series(((null >= obs).sum(0) + 1) / (NPERM + 1), index=core)
R2c = pd.Series(obs, index=core)
R["gene_R2_country_mean"] = float(obs.mean()); R["gene_R2_null_mean"] = float(null.mean())
from statsmodels.stats.multitest import multipletests
R["genes_sig_fdr05"] = int((multipletests(pgene, method="fdr_bh")[0]).sum())

# primary gene-level association
rho, p = spearmanr(M_cpl, excess); R["assoc_Mcpl_excess"] = [float(rho), float(p)]
rho2, p2 = spearmanr(M_mag, excess, nan_policy="omit"); R["assoc_Mmag_excess"] = [float(rho2), float(p2)]
R["prior_excess_means"] = [float(excess[M_prior].mean()), float(excess[~M_prior].mean())]

# set-level D statistic with permutation null that preserves gene-gene correlation
def Dstat(r2v, mask): return r2v[mask].mean() - r2v[~mask].mean()
hi = M_cpl >= M_cpl.quantile(2/3); lo = M_cpl <= M_cpl.quantile(1/3)
def D_hilo(r2v): return r2v[hi.values].mean() - r2v[lo.values].mean()
D_obs_prior = Dstat(obs, M_prior.values); D_obs_cpl = D_hilo(obs)
D_null_prior = np.array([Dstat(n, M_prior.values) for n in null]); D_null_cpl = np.array([D_hilo(n) for n in null])
R["D_prior"] = [float(D_obs_prior), float(((D_null_prior >= D_obs_prior).sum()+1)/(NPERM+1))]
R["D_cpl"] = [float(D_obs_cpl), float(((D_null_cpl >= D_obs_cpl).sum()+1)/(NPERM+1))]

# ---------------- variance components in multi-city countries ----------------
multi = [c for c, n in pm.groupby("country").city.nunique().items() if n >= 2]
R["multi_city_countries"] = multi
sub = pm.country.isin(multi)
def vc(v, cc, ci):
    v = v - v.mean(); tot = (v**2).sum()
    mc = v.groupby(cc).transform("mean"); mci = v.groupby(ci).transform("mean")
    return ((mc**2).sum()/tot, ((mci-mc)**2).sum()/tot, ((v-mci)**2).sum()/tot)
VC = pd.DataFrame({g: vc(LP.loc[sub, g], pm.country[sub], pm.city[sub]) for g in core}, index=["country","city_within","plant_within"]).T
R["VC_prior"] = VC[M_prior].mean().round(3).to_dict(); R["VC_other"] = VC[~M_prior].mean().round(3).to_dict()
R["VC_hiCpl"] = VC[hi].mean().round(3).to_dict(); R["VC_loCpl"] = VC[lo].mean().round(3).to_dict()
# ratio country / city-within per set
R["VC_ratio_prior_vs_other_MWU"] = float(mannwhitneyu((VC.country/(VC.city_within+1e-9))[M_prior],(VC.country/(VC.city_within+1e-9))[~M_prior]).pvalue)

# ---------------- country decoding (leave-one-city-out) ----------------
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import LeaveOneGroupOut, cross_val_predict
Xs, ys, gs = LP.loc[sub], pm.country[sub].values, pm.city[sub].values
def decode(genes):
    clf = make_pipeline(StandardScaler(), LogisticRegression(C=0.5, max_iter=3000, class_weight="balanced"))
    pr = cross_val_predict(clf, Xs[genes].values, ys, groups=gs, cv=LeaveOneGroupOut())
    return balanced_accuracy_score(ys, pr)
k = int(min(M_prior.sum(), (~M_prior).sum(), 30))
acc_prior = np.mean([decode(list(rng.choice(M_prior[M_prior].index, k, replace=False))) for _ in range(20)])
acc_other = np.mean([decode(list(rng.choice(M_prior[~M_prior].index, k, replace=False))) for _ in range(20)])
kc = int(min(hi.sum(), lo.sum()))
acc_hi = decode(list(hi[hi].index)); acc_lo = decode(list(lo[lo].index))
acc_rand = np.array([decode(list(rng.choice(core, k, replace=False))) for _ in range(100)])
R["decode"] = {"k": k, "prior": float(acc_prior), "other": float(acc_other), "hiCpl": float(acc_hi), "loCpl": float(acc_lo),
               "random_mean": float(acc_rand.mean()), "random_95": float(np.quantile(acc_rand, .95)),
               "chance": 1/len(multi), "n_plants": int(sub.sum())}

# ---------------- total-load invariance & transferability ----------------
lt = np.log10(meta.total)
R["load_sd_log10"] = float(lt.std()); R["load_fold_IQR"] = float(10**(lt.quantile(.75)-lt.quantile(.25)))
for lev in ["country","city","plant"]:
    R[f"load_R2_{lev}"] = float(1 - ((lt - lt.groupby(meta[lev]).transform("mean"))**2).sum()/((lt-lt.mean())**2).sum())
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold
mcols = [c for c in s.columns[:s.columns.get_loc("ARG_per_cell")] if pd.api.types.is_numeric_dtype(s[c])
         and s[c].isna().mean() < 0.3 and c not in ["GC","rrn.copy.mean","comm.mean.copy"]]
Xm = s[mcols].fillna(s[mcols].median())
prior_genes = M_prior[M_prior].index; other_genes = M_prior[~M_prior].index
targets = {"total": lt, "mobile_layer": np.log10(G[prior_genes].sum(1)), "intrinsic_layer": np.log10(G[other_genes].sum(1)),
           "mobile_fraction": G[prior_genes].sum(1)/G[core].sum(1)}
TR = {}
for tn, yv in targets.items():
    for lev in ["plant","city","country"]:
        pr = cross_val_predict(RandomForestRegressor(500, min_samples_leaf=3, random_state=0, n_jobs=-1),
                               Xm, yv, groups=meta[lev], cv=GroupKFold(10))
        TR[f"{tn}|{lev}"] = float(1 - ((yv-pr)**2).sum()/((yv-yv.mean())**2).sum())
R["transfer_R2"] = TR
R["layer_sd_log10"] = {k2: float(v.std()) for k2, v in targets.items() if k2 != "mobile_fraction"}
R["mobile_fraction_range"] = [float(targets["mobile_fraction"].min()), float(targets["mobile_fraction"].median()), float(targets["mobile_fraction"].max())]

# ---------------- replicate stability (ICC among plants with >=2 samples) ----------------
def icc1(v, grp):
    d = pd.DataFrame({"v": v, "g": grp}); d = d[d.g.map(d.g.value_counts()) >= 2]
    k0 = d.g.value_counts(); n = len(k0); N = len(d)
    ms_b = (k0 * (d.groupby("g").v.mean() - d.v.mean())**2).sum()/(n-1)
    ms_w = ((d.v - d.groupby("g").v.transform("mean"))**2).sum()/(N-n)
    kbar = (N - (k0**2).sum()/N)/(n-1)
    return (ms_b - ms_w)/(ms_b + (kbar-1)*ms_w), n
R["ICC"] = {k2: [float(icc1(v, meta.plant)[0]), icc1(v, meta.plant)[1]] for k2, v in targets.items()}

# ---------------- socio-economic gradient ----------------
import statsmodels.formula.api as smf
cm = pd.DataFrame({"mf": targets["mobile_fraction"], "intr": targets["intrinsic_layer"], "mob": targets["mobile_layer"],
                   "gdp": np.log10(meta.gdp), "pop": np.log10(meta["pop"]), "country": meta.country, "city": meta.city})
cty = cm.groupby("city").agg({"mf":"mean","intr":"mean","mob":"mean","gdp":"mean","pop":"mean","country":"first"}).dropna()
cty["lmf"] = np.log10(cty.mf)
SE = {}
for yv in ["lmf","intr","mob"]:
    SE[f"{yv}_city_rho_gdp"] = list(map(float, spearmanr(cty.gdp, cty[yv])[:2]))
    m = smf.mixedlm(f"{yv} ~ gdp + pop", cty, groups=cty.country).fit(reml=True)
    SE[f"{yv}_mixed_gdp"] = [float(m.params["gdp"]), float(m.pvalues["gdp"])]
    SE[f"{yv}_mixed_pop"] = [float(m.params["pop"]), float(m.pvalues["pop"])]
cn = cty.groupby("country").agg({"lmf":"mean","gdp":"mean","intr":"mean"})
SE["country_rho_gdp_lmf"] = list(map(float, spearmanr(cn.gdp, cn.lmf)[:2]))
SE["country_rho_gdp_intr"] = list(map(float, spearmanr(cn.gdp, cn.intr)[:2]))
R["socioeconomic"] = SE

# ---------------- host breadth of mobile vs non-mobile MAG-borne ARGs ----------------
HB = {}
for c0 in ["Beta-lactam","Glycopeptide","Tetracycline","Aminoglycoside"]:
    d = mg[mg.drug_class == c0]; m1 = d[d.mobile == 1].group.values; m0 = d[d.mobile == 0].group.values
    n = min(len(m1), len(m0))
    if n < 20: continue
    def rich(v): return len(set(rng.choice(v, n, replace=False)))
    def shan(v):
        p = pd.Series(rng.choice(v, n, replace=False)).value_counts(normalize=True).values; return -(p*np.log(p)).sum()
    r1 = [rich(m1) for _ in range(500)]; r0 = [rich(m0) for _ in range(500)]
    h1 = [shan(m1) for _ in range(500)]; h0 = [shan(m0) for _ in range(500)]
    HB[c0] = {"n_mobile": len(m1), "n_nonmobile": len(m0), "n_rarefied": n,
              "rich_mobile": float(np.mean(r1)), "rich_non": float(np.mean(r0)),
              "shannon_mobile": float(np.mean(h1)), "shannon_non": float(np.mean(h0)),
              "p_shannon_mobile_gt": float(np.mean(np.array(h1) <= np.array(h0)))}
R["host_breadth"] = HB

# ---------------- robustness of D (prior set) ----------------
def D_for(mask_rows, Yl=None):
    P2 = pm[mask_rows]; Y2 = (LP if Yl is None else Yl)[mask_rows].values
    H2 = pd.get_dummies(P2.country).values.astype(float)
    o = r2_country(Y2, H2); cc2 = P2.groupby("city").country.first()
    nl = []
    for _ in range(499):
        lab = pd.Series(rng.permutation(cc2.values), index=cc2.index)
        nl.append(Dstat(r2_country(Y2, pd.get_dummies(P2.city.map(lab)).values.astype(float)), M_prior.values))
    d0 = Dstat(o, M_prior.values); nl = np.array(nl)
    return [float(d0), float(((nl >= d0).sum()+1)/500)]
RB = {}
RB["exclude_USA_China"] = D_for(~pm.country.isin(["USA","China"]))
LPgc = LP.apply(lambda col: pd.Series(resid(col.values, pm.GC.values), index=col.index))
RB["GC_adjusted"] = D_for(pd.Series(True, index=pm.index), LPgc)
caps = []
for _ in range(50):
    keep = pm.groupby("country", group_keys=False).apply(lambda d: d.sample(min(len(d), 4), random_state=int(rng.integers(1e9)))).index
    Y2 = LP.loc[keep].values; H2 = pd.get_dummies(pm.loc[keep, "country"]).values.astype(float)
    caps.append(Dstat(r2_country(Y2, H2), M_prior.values))
RB["cap4_per_country_D_mean_and_frac_pos"] = [float(np.mean(caps)), float(np.mean(np.array(caps) > 0))]
s5 = pd.read_excel(x, "Fig. S5b", index_col=0)
RB["contig_vs_read_rho"] = float(spearmanr(s5.contig_based, s5.reads_based)[0])
# alternative pseudo-count / prevalence threshold
core2 = prev[prev > 0.8].index.tolist(); m2 = M_prior.reindex(core2).fillna(False).values
if len(core2) > 20:
    LP2 = np.log10(G[core2] + G[core2].replace(0,np.nan).min()/2).groupby(meta.plant).mean().loc[pm.index]
    o2 = r2_country(LP2.values, Hc); RB["prev80_D"] = [float(Dstat(o2, m2)), len(core2)]
R["robustness"] = RB

# ---------------- save ----------------
gene_tab = pd.DataFrame({"prevalence": prev[core], "class": [";".join(gcls.get(g, [])) for g in core],
                         "mechanism": [";".join(gmech.get(g, [])) for g in core], "M_prior": M_prior,
                         "M_coupling": M_cpl, "M_MAG": M_mag, "R2_country": R2c, "excess_R2": excess, "perm_p": pgene,
                         "VC_country": VC.country, "VC_city_within": VC.city_within, "VC_plant_within": VC.plant_within})
gene_tab.to_csv("gene_table.csv")
np.save("null_Dprior.npy", D_null_prior); np.save("null_Dcpl.npy", D_null_cpl)
np.save("acc_rand.npy", acc_rand)
cty.to_csv("city_table.csv")
pd.DataFrame(targets).assign(country=meta.country, plant=meta.plant).to_csv("sample_layers.csv")
json.dump(R, open("results.json", "w"), indent=1, default=str)
print(json.dumps(R, indent=1, default=str))
