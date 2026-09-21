"""Assembles Supplementary Tables S1 and S2 from the analysis outputs -> S1_tmp.csv, S2_tmp.csv (read by make_xlsx.py)."""
import json, pandas as pd
K=json.load(open('canonical.json')); D=json.load(open('decode.json')); R=json.load(open('results.json'))
C=json.load(open('controls.json')); MD=json.load(open('module.json')); EF=json.load(open('effort.json')); FR=json.load(open('frozen_results.json')); IJ=json.load(open('integron_depth.json')); F2=json.load(open('frozen_results2.json')); F3=json.load(open('frozen_results3.json')); FA=json.load(open('frozen_results3_alr.json')); FS3=json.load(open('f_specificity.json')); TH=json.load(open('test_H.json')); SP=json.load(open('frozen_spec.json'))
# ---- S1 ----
d=pd.read_csv('gene_table_final.csv').rename(columns={'Unnamed: 0':'gene'})
d=d[['gene','class','mechanism','prevalence','mean_log','sd_log','M_prior','M_coupling','M_MAG','R2_country','excess_R2','perm_p','VC_country','VC_city_within','VC_plant_within']]
d.columns=['gene','drug_class','mechanism','prevalence','mean_log10_abundance','plant_level_SD','mobile_associated','MGE_coupling_index','MAG_mobility_index','country_R2','excess_country_R2','permutation_P','SS_share_country','SS_share_city_within','SS_share_plant_within']
d=d.round(4)
mm=pd.read_csv('module_membership.csv')
s1=d.merge(mm[['gene','cluster','in_module']].rename(columns={'in_module':'in_label_free_module'}),on='gene')
ri=pd.read_csv('intI1_rho.csv',index_col=0)['rho_intI1_within_country'].round(4)
s1['rho_intI1_within_country']=s1.gene.map(ri)
assert s1.rho_intI1_within_country.notna().all()
# ---- S2 ----
T=K['transfer_R2']
rows=[("Primary D (observed country R2, mobile - other)",K['D'][0],K['D'][1],"9999 city-block permutations"),
("Difference in excess country R2",K['D_excess_diff'],None,""),
("Bootstrap 95% CI for D",f"{K['D_boot_CI'][0]:.3f}, {K['D_boot_CI'][1]:.3f}",None,"1000 city-cluster resamples"),
("Mean country R2, observed / permutation null",f"{K['mean_R2_obs']:.3f} / {K['mean_R2_null']:.3f}",None,"marginal R2 of country"),
("Genes with significant country structure (FDR<0.05)",K['n_sig_fdr'],None,"of 114 core genes"),
("beta(mobile) adjusted for abundance, prevalence, SD",K['beta_mob'][0],K['beta_mob'][1],"permutation P"),
("D, abundance-matched gene pairs",K['D_matched'][0],K['D_matched'][1],"median |delta log10| = %.2f"%K['matched_absdiff_median']),
("D, abundance-stratified (quintiles)",K['D_strat'][0],K['D_strat'][1],""),
("D, within drug class (4 classes)",K['D_class'][0],K['D_class'][1],""),
("D, centred log-ratio data",K['D_CLR'][0],K['D_CLR'][1],""),
("D, GC-content adjusted",K['D_GC'][0],K['D_GC'][1],""),
("D, MGE-coupling tertiles (high - low)",K['D_tertile'][0],K['D_tertile'][1],""),
("D, genes with prevalence >80%",K['D_prev80'][0],None,"%d genes"%K['D_prev80'][1]),
("D, excluding USA",K['D_noUSA'][0],K['D_noUSA'][1],""),
("D, excluding China",K['D_noCHN'][0],K['D_noCHN'][1],""),
("D, excluding USA and China",K['D_noBoth'][0],K['D_noBoth'][1],"not significant"),
("D, country jackknife range",f"{K['jk_range'][0]:.3f} to {K['jk_range'][1]:.3f}",None,"min = China excluded"),
("D, leave-one-gene-family-out range",f"{K['LOFO_range'][0]:.3f} to {K['LOFO_range'][1]:.3f}",None,""),
("rho(MGE coupling, excess R2)",K['rho_coupling'][0],K['rho_coupling'][1],"permutation P; not significant"),
("Median MGE coupling, mobile / other",f"{K['coupling_median'][0]:.2f} / {K['coupling_median'][1]:.2f}",None,""),
("SS share between countries, mobile / other",f"{K['SS_country'][0]:.2f} / {K['SS_country'][1]:.2f}",K['SS_contrast_p'],"7 multi-city countries; contrast P"),
("SS share between cities, mobile / other",f"{K['SS_city'][0]:.2f} / {K['SS_city'][1]:.2f}",None,""),
("Country decoding accuracy, mobile / other",f"{D['mobile']:.2f} / {D['other']:.2f}",D['p'],"50 draws of 30 genes; chance %.2f"%D['chance']),
("Country decoding, random 30-gene draws",D['random_mean'],None,"contains both gene types"),
("Layer SD, mobile / non-mobile (log10)",f"{K['layer_sd'][0]:.3f} / {K['layer_sd'][1]:.3f}",None,""),
("Layer SD, abundance-matched non-mobile sets",K['sd_matched'][0],None,"95%% range %.3f-%.3f"%(K['sd_matched'][1],K['sd_matched'][2])),
("Mobile layer SD / matched-set median",K['sd_ratio_matched'],None,""),
("Layer SD, random 45-gene sets (median)",K['sd_random45'][0],K['sd_random45'][1],"P for mobile layer"),
("Per-gene SD, mobile / other (median)",f"{K['per_gene_sd'][0]:.2f} / {K['per_gene_sd'][1]:.2f}",K['per_gene_sd'][2],"Mann-Whitney"),
("Mean pairwise r, mobile / other",f"{K['r_within'][0]:.2f} / {K['r_within'][1]:.2f}",None,"random 45-gene median %.2f"%K['r_random45'][0]),
("Non-mobile layer: top-5 gene share",K['oth_top5_share'],None,"SD without top 5 = %.3f"%K['oth_sd_wo_top5']),
("Non-mobile share of core ARG copies (median)",K['share_oth_median'],None,""),
("Mobile share min / median / max",f"{100*K['mobile_share'][0]:.1f}% / {100*K['mobile_share'][1]:.1f}% / {100*K['mobile_share'][2]:.1f}%",None,""),
("REML country variance fraction, mobile / non-mobile",f"{K['REML']['mob']['frac']['country']:.2f} / {K['REML']['oth']['frac']['country']:.2f}",None,""),
("REML between-country variance ratio",K['REML_ratio'],None,"mobile / non-mobile, absolute components"),
("ARG load SD (log10) / IQR fold-range",f"{R['load_sd_log10']:.3f} / {R['load_fold_IQR']:.2f}",None,""),
("Transferability R2, total (plant/city/country)",f"{T['total|plant']:.3f} / {T['total|city']:.3f} / {T['total|country']:.3f}",None,"imputation fitted within folds"),
("Transferability R2, mobile layer (plant/city/country)",f"{T['mobile_layer|plant']:.3f} / {T['mobile_layer|city']:.3f} / {T['mobile_layer|country']:.3f}",None,"imputation fitted within folds"),
("Transferability R2, non-mobile layer (plant/city/country)",f"{T['intrinsic_layer|plant']:.3f} / {T['intrinsic_layer|city']:.3f} / {T['intrinsic_layer|country']:.3f}",None,"imputation fitted within folds"),
("Contig- vs read-based ARG totals (Spearman)",R['robustness']['contig_vs_read_rho'],None,""),
("Label-free clustering: number of clusters (max silhouette)",MD['k'],None,"within-country correlations, average linkage"),
("Label-free module size / mobile-associated members",f"{MD['module']['size']} / {MD['module']['n_mobile']}",f"q={MD['module']['q']:.1e}","hypergeometric enrichment, BH"),
("Label-free module precision / recall",f"{MD['module']['precision']:.2f} / {MD['module']['recall']:.2f}",None,"vs curated 45-gene set"),
("Label-free module D vs other genes",round(MD['D_module'][0],4),MD['D_module'][1],"999 city-block permutations"),
("Label-free: non-annotated module members vs non-mobile outside",round(MD['D_nonmobile_in_vs_out'][0],4),MD['D_nonmobile_in_vs_out'][1],f"{MD['D_nonmobile_in_vs_out'][2]} genes"),
("Label-free: within-country mean r, module / rest",f"{MD['module_mean_within_r']:.3f} / {MD['rest_mean_within_r']:.3f}",None,""),
("Label-free: module co-clustering in bootstraps",round(MD['module_cocluster'],3),None,"200 plant bootstraps"),
("Label-free: bootstrap ARI median (5-95%)",f"{MD['bootstrap_ARI'][0]:.2f} ({MD['bootstrap_ARI'][1]:.2f}-{MD['bootstrap_ARI'][2]:.2f})",None,""),
("Label-free: min Jaccard with primary across sensitivity runs",round(min(v['jaccard_with_primary'] for v in MD['sensitivity'].values()),3),None,"complete linkage, k-1, k+1, raw correlations"),
("D, adjusted for sampling effort (proxy A: max sequence number, source Fig. S1a)",round(EF['D_adj_effA'][0],4),EF['D_adj_effA'][1],"999 city-block permutations"),
("D, adjusted for sampling effort (proxy C: max sequence number, source Fig. S1c)",round(EF['D_adj_effC'][0],4),EF['D_adj_effC'][1],"999 city-block permutations"),
("D, adjusted for sampling effort (both proxies)",round(EF['D_adj_both'][0],4),EF['D_adj_both'][1],"999 city-block permutations"),
("Sampling effort: country marginal R2 (proxy A / C)",f"{EF['effA_country_R2']:.2f} / {EF['effC_country_R2']:.2f}",None,"plant level"),
("Sampling effort vs mobile share: Spearman rho (proxy A / C)",f"{EF['effA_rho_mobile_share']:.2f} / {EF['effC_rho_mobile_share']:.2f}",None,"plant level"),
("D, adjusted for sequencing depth (log10 HQ bases)",round(FR['C_depth']['D_depth_adjusted'],4),FR['C_depth']['P'],"999 city-block permutations"),
("Sequencing depth: median (range) HQ Gb per sample",f"{FR['C_depth']['HQ_Gb_median']:.1f} ({FR['C_depth']['HQ_Gb_range'][0]:.1f}-{FR['C_depth']['HQ_Gb_range'][1]:.1f})",None,f"country marginal R2 = {FR['C_depth']['depth_country_R2']:.2f}"),
("PRE-SPECIFIED A: module co-membership in Martiny et al. sewage network (share of pairs)",f"{FR['A_primary']['observed_share']:.3f} vs null {FR['A_primary']['null_mean']:.3f}",FR['A_primary']['P'],f"{len(SP['A_network_test']['module_genes_mapped'])} genes, {FR['A_n_module_nodes']} nodes; not supported"),
("PRE-SPECIFIED A (secondary): acquired-acquired link density",f"{FR['A_secondary']['observed_density']:.3f} vs null {FR['A_secondary']['null_mean']:.3f}",FR['A_secondary']['P'],"10,000 random node sets"),
("PRE-SPECIFIED B: within-country rho with intI1, mobile / other",f"{FR['B_primary']['mean_rho_mobile']:.3f} / {FR['B_primary']['mean_rho_other']:.3f}",FR['B_primary']['P'],"9,999 within-country permutations; supported"),
("PRE-SPECIFIED B (secondary): within-country rho with intI1, module / outside",f"{FR['B_secondary']['mean_rho_module']:.3f} / {FR['B_secondary']['mean_rho_outside']:.3f}",FR['B_secondary']['P'],"9,999 within-country permutations"),
("Frozen specification SHA-256",FR['spec_sha256'],None,"frozen_spec.json in Supplementary Code"),
("EXTERNAL, CLR sensitivity (sewage): module vs other acquired genes, mean country R2",f"{F3['D_validation_primary']['mean_R2_module']:.3f} / {F3['D_validation_primary']['mean_R2_acquired_other']:.3f}",F3['D_validation_primary']['P'],f"D = {F3['D_validation_primary']['D_ext']:.4f}; NOT SUPPORTED; {F3['n_samples']} samples, {F3['n_countries']} countries"),
("EXTERNAL, CLR sensitivity (sewage): module vs latent (functional) genes",round(F3['D_validation_secondary_vs_FG']['D_ext'],4),F3['D_validation_secondary_vs_FG']['P'],"largely restates Martiny et al. acquired vs latent"),
("EXTERNAL, CLR sensitivity (sewage): within-country co-variation of module vs matched acquired sets",f"{F3['covariation']['mean_within_country_r_module']:.3f} / {F3['covariation']['null_median']:.3f}",F3['covariation']['P'],"SUPPORTED; 1,000 abundance-matched sets"),
("EXTERNAL, CLR sensitivity (sewage): temporal stability, first vs last sampling year",round(F3['E_temporal']['rho_first_vs_last_year'],3),F3['E_temporal']['P'],f"SUPPORTED; {F3['E_temporal']['n_cities_multi_year']} cities"),
("EXTERNAL, CLR sensitivity (sewage): variance fractions of module score (city/year/country/resid)",str(F3['E_temporal']['variance_fractions']),None,"REML"),
("EXTERNAL: module share, sewage vs activated sludge across shared cities",round(F3['G_cross_compartment']['rho'],3),F3['G_cross_compartment']['P'],f"exploratory; {F3['G_cross_compartment']['n_cities']} cities"),
("Frozen specification 3 SHA-256",F3['spec_sha256'],None,"frozen_spec3.json in Supplementary Code"),
("EXTERNAL, PRE-SPECIFIED ALR (sewage): module vs other acquired genes, mean country R2",f"{FA['D_validation_primary']['mean_R2_module']:.3f} / {FA['D_validation_primary']['mean_R2_acquired_other']:.3f}",FA['D_validation_primary']['P'],f"D = {FA['D_validation_primary']['D_ext']:.4f}; NOT SUPPORTED; {FA['n_samples']} samples, {FA['n_countries']} countries"),
("EXTERNAL, ALR: module vs latent (functional) genes",round(FA['D_validation_secondary_vs_FG']['D_ext'],4),FA['D_validation_secondary_vs_FG']['P'],"largely restates Martiny et al. acquired vs latent"),
("EXTERNAL, ALR: within-country co-variation of module vs matched acquired sets",f"{FA['covariation']['mean_within_country_r_module']:.3f} / {FA['covariation']['null_median']:.3f}",FA['covariation']['P'],"SUPPORTED; 1,000 abundance-matched sets"),
("EXTERNAL, ALR: temporal stability, first vs last sampling year",round(FA['E_temporal']['rho_first_vs_last_year'],3),FA['E_temporal']['P'],f"SUPPORTED; {FA['E_temporal']['n_cities_multi_year']} cities"),
("EXTERNAL, ALR: variance fractions of module score (city/year/country/resid)",str(FA['E_temporal']['variance_fractions']),None,"REML"),
("PRE-SPECIFIED F: module score vs national antibiotic consumption",round(FA['F_drivers']['rho'],3),FA['F_drivers']['P'],f"SUPPORTED; {FA['F_drivers']['n_countries_matched']} countries; WHO GLASS via Our World in Data (deviation from Browne et al.)"),
("F secondary: mixed model with world region as random intercept (coefficient)",round(FA['F_drivers']['mixed_region']['coef'],4),FA['F_drivers']['mixed_region']['P'],"weaker after region adjustment"),
("F secondary: partial correlation controlling GDP",None,None,"not run: no GDP data for the sewage cohort"),
("EXPLORATORY (after F): consumption vs acquired genes outside module",round(FS3['acquired_not_module']['rho'],3),FS3['acquired_not_module']['P'],"association not module-specific"),
("EXPLORATORY (after F): consumption vs latent genes",round(FS3['functional_latent']['rho'],3),FS3['functional_latent']['P'],""),
("PRE-SPECIFIED H: D adjusted for run-level fields (read length, release batch)",round(TH['H_primary_run_level'][0],4),TH['H_primary_run_level'][1],f"SUPPORTED; instrument, centre and layout constant; {TH['run_factors']['plants_with_301bp_runs']} plants with 301 bp runs, {TH['run_factors']['plants_in_second_release']} in second release"),
("ADDED: country share of variance in collection date (discovery cohort)",round(TH['date_country_R2'],3),None,"country and sampling season nearly aliased"),
("ADDED: D adjusted for collection date",round(TH['added_D_after_date_adjustment'][0],4),TH['added_D_after_date_adjustment'][1],"halved"),
("ADDED: date sensitivity, mobile minus other, USA",round(TH['added_within_country_date_sensitivity']['USA']['difference'],4),TH['added_within_country_date_sensitivity']['USA']['P'],"999 within-country permutations"),
("ADDED: date sensitivity, mobile minus other, China",round(TH['added_within_country_date_sensitivity']['China']['difference'],4),TH['added_within_country_date_sensitivity']['China']['P'],"999 within-country permutations"),
("PRE-SPECIFIED C: ResFinder-acquired vs other genes, mean country R2",f"{F2['mean_R2']['resfinder']:.3f} / {F2['mean_R2']['other']:.3f}",F2['primary_D_RF'][1],f"D = {F2['primary_D_RF'][0]:.3f}; {F2['n_acquired']} vs {F2['n_other']} genes; 9,999 city-block permutations"),
("PRE-SPECIFIED C (secondary): agreement with curated list, Cohen kappa",round(F2['agreement']['kappa'],3),None,f"both {F2['agreement']['both']}, ResFinder only {F2['agreement']['resfinder_only']}, curated only {F2['agreement']['curated_only']}, neither {F2['agreement']['neither']}"),
("PRE-SPECIFIED C (secondary): ResFinder-acquired vs other among non-curated genes (D)",round(F2['D_RF_noncurated'][0],4),F2['D_RF_noncurated'][1],f"{F2['n_noncurated_split'][0]} vs {F2['n_noncurated_split'][1]} genes"),
("PRE-SPECIFIED C (secondary): ResFinder-acquired genes in label-free module",f"{F2['module_enrichment']['resfinder_in_module']} of {F2['n_acquired']}",F2['module_enrichment']['P'],"hypergeometric"),
("Frozen specification 2 SHA-256",F2['spec_sha256'],None,"frozen_spec2.json in Supplementary Code"),
("EXPLORATORY (post hoc): intI1 partial rho controlling total load, mobile / other",f"{IJ['intI1']['mean_r_mobile']:.3f} / {IJ['intI1']['mean_r_other']:.3f}",IJ['intI1']['mobile_vs_other'][1],"1,999 within-country permutations"),
("EXPLORATORY: intI1 partial rho, excluding sul1 (difference)",round(IJ['intI1']['mobile_vs_other_no_sul1'][0],4),IJ['intI1']['mobile_vs_other_no_sul1'][1],"1,999 within-country permutations"),
("EXPLORATORY: intI1 partial rho, non-annotated module members vs non-mobile outside (difference)",round(IJ['intI1']['module_nonannot_vs_outside'][0],4),IJ['intI1']['module_nonannot_vs_outside'][1],"1,999 within-country permutations"),
("EXPLORATORY: qacEdelta1 partial rho, mobile / other",f"{IJ['qacEdelta1']['mean_r_mobile']:.3f} / {IJ['qacEdelta1']['mean_r_other']:.3f}",IJ['qacEdelta1']['mobile_vs_other'][1],"1,999 within-country permutations"),
("EXPLORATORY: qacEdelta1 partial rho, non-annotated module members vs non-mobile outside (difference)",round(IJ['qacEdelta1']['module_nonannot_vs_outside'][0],4),IJ['qacEdelta1']['module_nonannot_vs_outside'][1],"1,999 within-country permutations")]
s2=pd.DataFrame(rows,columns=['analysis','estimate','P','note'])
s2['estimate']=s2.estimate.apply(lambda v: round(v,4) if isinstance(v,float) else v)
s2['P']=s2.P.apply(lambda v: ('<0.001' if isinstance(v,float) and v<0.001 else (round(v,4) if isinstance(v,float) else v)))
s1.to_csv('S1_tmp.csv',index=False); s2.to_csv('S2_tmp.csv',index=False)
print('built',s1.shape,s2.shape)
