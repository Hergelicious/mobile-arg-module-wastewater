"""FROZEN SPECIFICATION 2 (external classification test). Written and hashed before the test statistic is computed.
Honest status: per-gene country R2 values were already known from the primary analysis; what this test adds is a
classification defined by an external database (ResFinder, via Martiny et al. 2025 Supplementary Data 2), applied by a
mechanical rule, instead of our own curated list."""
import json, re, hashlib, pandas as pd
g=pd.read_csv('gene_table_final.csv',index_col=0); core=list(g.index)
sd2=pd.read_excel('/mnt/user-data/uploads/41467_2025_66070_MOESM5_ESM.xlsx')
rf=sd2[sd2.database=='resfinder'].fa_name.astype(str)
stems={re.sub(r'_\d+_[^_]+$','',n.split('|',1)[1]).lower() for n in rf}
def norm_card(x):   # identical to freeze_spec.py (frozen test A)
    y=x.strip()
    m=re.fullmatch(r'(Erm|Ere|Mph|Msr|Lnu|Mef|Lsa|Vat|Vga)([A-Z])(\d*)',y,re.I)
    if m: return f"{m.group(1).lower()}({m.group(2).lower()}){m.group(3)}"
    m=re.fullmatch(r'(OXA|TEM|GES|VEB|PER|CARB|CTX-M|SHV|IMP|VIM|NDM|KPC|CMY|ACT|LCR|NPS)-?(\d+)',y,re.I)
    if m: return f"bla{m.group(1).lower()}-{m.group(2)}"
    return y.lower()
present=sorted([x for x in core if norm_card(x) in stems])
spec={'version':'frozen-2',
 'D_external_classification':{
  'data':'Martiny et al. 2025 Supplementary Data 2 (ARG metadata): entries with database == resfinder (%d entries, %d gene stems)'%(len(rf),len(stems)),
  'rule':'a core gene is ResFinder-acquired if its name, normalised by the function used in frozen test A, exactly equals the stem of any ResFinder entry (allele/accession suffix removed); no manual edits',
  'resfinder_acquired_genes':present,
  'n_acquired':len(present),'n_other':len(core)-len(present),
  'primary_statistic':'D_RF = mean observed plant-level country R2 of ResFinder-acquired genes minus that of all other core genes',
  'null':'9,999 permutations reassigning whole cities among countries (as for the primary D); one-sided P',
  'secondary':['Cohen kappa between ResFinder classification and the curated mobile-associated list',
               'D_RF restricted to genes NOT on the curated list (ResFinder-acquired vs non-acquired among the 69 non-curated genes), same null',
               'hypergeometric enrichment of ResFinder-acquired genes in the label-free module'],
  'success_criterion':'primary P<0.05 in the predicted (positive) direction',
  'uninformative_if':'fewer than 10 core genes classified as ResFinder-acquired',
  'reporting':'reported in the manuscript whatever the outcome'}}
txt=json.dumps(spec,indent=1,sort_keys=True)
open('frozen_spec2.json','w').write(txt)
print('n ResFinder-acquired core genes:',len(present)); print(present)
print('SHA-256:',hashlib.sha256(txt.encode()).hexdigest())
