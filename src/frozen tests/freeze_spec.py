"""FROZEN SPECIFICATION — written before reading any outcome values:
   - Martiny SD3 community labels (comm_id) and SD4 link weights have NOT been read.
   - Zhu Supplementary Data 3 (MGE abundances) values have NOT been read.
Only gene/node NAMES were inspected to build the mapping below."""
import json, re, hashlib, pandas as pd
g=pd.read_csv('gene_table_final.csv',index_col=0); mm=pd.read_csv('module_membership.csv').set_index('gene')
core=list(g.index)
nodes=pd.read_excel('/mnt/user-data/uploads/41467_2025_66070_MOESM6_ESM.xlsx',usecols=['taxa_id','source','type','name_extracted_amr'])
acq=nodes[(nodes.source=='amr')&(nodes.type=='acquired')]
def norm_card(x):
    """CARD subtype name -> ResFinder-style gene token (exact-gene matching, no family collapsing)."""
    y=x.strip()
    m=re.fullmatch(r'(Erm|Ere|Mph|Msr|Lnu|Mef|Lsa|Vat|Vga)([A-Z])(\d*)',y,re.I)
    if m: return f"{m.group(1).lower()}({m.group(2).lower()}){m.group(3)}"
    m=re.fullmatch(r'(OXA|TEM|GES|VEB|PER|CARB|CTX-M|SHV|IMP|VIM|NDM|KPC|CMY|ACT|LCR|NPS)-?(\d+)',y,re.I)
    if m: return f"bla{m.group(1).lower()}-{m.group(2)}"
    return y.lower()
def node_members(name):
    return {re.sub(r'_\d+_[^_]+$','',p) for p in str(name).split('/')}
mapping={}
for gene in core:
    tok=norm_card(gene); hits=[r.taxa_id for r in acq.itertuples() if tok in node_members(r.name_extracted_amr)]
    if len(hits)==1: mapping[gene]=hits[0]
spec={
 'version':'frozen-1',
 'A_network_test':{
   'data':'Martiny et al. 2025 SD3 (nodes) and SD4 (links); acquired ARG nodes only (n=%d)'%len(acq),
   'mapping_rule':'exact gene match after CARD->ResFinder name normalisation; a gene maps only if it matches exactly one node; functional nodes excluded (named by clone IDs, not mappable)',
   'mapping':mapping,
   'module_genes_mapped':sorted([k for k in mapping if mm.loc[k,'in_module']]),
   'nonmodule_genes_mapped':sorted([k for k in mapping if not mm.loc[k,'in_module']]),
   'primary_statistic':'fraction of pairs of distinct mapped module nodes sharing the same comm_id',
   'null':'10,000 random draws of the same number of distinct nodes from all acquired ARG nodes; one-sided P=(#null>=obs+1)/10001',
   'secondary':'density of acquired-acquired links among mapped module nodes vs the same null',
   'success_criterion':'primary P<0.05',
   'uninformative_if':'fewer than 6 distinct module nodes mapped',
   'reporting':'reported in the manuscript whatever the outcome'},
 'B_intI1_test':{
   'data':'Zhu et al. 2025 Supplementary Data 3 (MGE abundances per sample)',
   'intI1_definition':"rows whose MGE_Gene_name equals 'intI1' (case-insensitive), summed per sample; if absent the test is not run",
   'transform':'plant means of log10(x + half the smallest positive value); country means removed from genes and intI1',
   'primary_statistic':'mean within-country Spearman correlation with intI1 for the 45 mobile-associated genes minus that for the 69 other core genes',
   'null':'9,999 permutations of intI1 values among plants within country; one-sided P',
   'secondary':'same statistic for the 77 label-free module genes versus the 37 genes outside the module',
   'success_criterion':'primary P<0.05 in the predicted (positive) direction',
   'reporting':'reported in the manuscript whatever the outcome'},
 'C_depth_adjustment':'Covariate analysis, not an outcome test: D recomputed after regressing each gene (plant level) on log10 HQ bases from Zhu Supplementary Data 1'}
txt=json.dumps(spec,indent=1,sort_keys=True)
open('frozen_spec.json','w').write(txt)
print(txt[:1400]); print('\nSHA-256 of frozen_spec.json:',hashlib.sha256(txt.encode()).hexdigest())
