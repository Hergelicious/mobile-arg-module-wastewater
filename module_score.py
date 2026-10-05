#!/usr/bin/env python3
"""Compute the wastewater resistome module metrics from any ARG abundance table.

Usage:  python3 module_score.py abundance.csv [--panel resistome_module_panel.csv] [--out scores.csv]

Input : CSV with genes as rows (first column = gene name, in CARD, ResFinder or PanRes nomenclature) and one
        column per sample holding abundances in any consistent unit (per cell, FPKM, counts, ALR input).
Output: per sample - module_score (mean log abundance of matched module genes), module_share (module / total),
        mobile_share (curated mobile-associated / total), and the number of panel genes matched.
Notes : the module is fixed (frozen_spec3.json); this tool only matches names and aggregates. A sample with
        fewer than 10 matched module genes is reported with a warning, because the score is then unreliable.
"""
import argparse, re, sys, numpy as np, pandas as pd
def norm(x):
    y=str(x).strip()
    m=re.fullmatch(r'(Erm|Ere|Mph|Msr|Lnu|Mef|Lsa|Vat|Vga)([A-Z])(\d*)',y,re.I)
    if m: return f"{m.group(1).lower()}({m.group(2).lower()}){m.group(3)}"
    m=re.fullmatch(r'(OXA|TEM|GES|VEB|PER|CARB|CTX-M|SHV|IMP|VIM|NDM|KPC|CMY|ACT|LCR|NPS)-?(\d+)',y,re.I)
    if m: return f"bla{m.group(1).lower()}-{m.group(2)}"
    return y.lower()
def score(table,panel):
    key={}
    for _,r in panel.iterrows():
        for k in [norm(r.gene_CARD), r.name_normalised]+([p for p in str(r.panres_ids).split(';') if p] if pd.notna(r.panres_ids) else []):
            key[str(k).lower()]=r.gene_CARD
    idx=[key.get(norm(i),key.get(str(i).lower())) for i in table.index]
    t=table.copy(); t['_panel']=idx; t=t[t._panel.notna()].groupby('_panel').sum()
    pm=panel.set_index('gene_CARD').loc[t.index]
    mod=t[pm.in_module.values]; mob=t[pm.mobile_associated_curated.values]
    pc=t[t>0].min().min()/2
    out=pd.DataFrame({'module_score':np.log10(mod+pc).mean(),'module_share':mod.sum()/t.sum(),
                      'mobile_share':mob.sum()/t.sum(),'panel_genes_matched':int(len(t)),'module_genes_matched':int(len(mod))})
    if len(mod)<10: print(f'WARNING: only {len(mod)} module genes matched; scores are unreliable',file=sys.stderr)
    return out
if __name__=='__main__':
    a=argparse.ArgumentParser(); a.add_argument('table'); a.add_argument('--panel',default='resistome_module_panel.csv'); a.add_argument('--out',default='module_scores.csv')
    n=a.parse_args()
    tab=pd.read_csv(n.table,index_col=0)
    s=score(tab,pd.read_csv(n.panel))
    s.to_csv(n.out); print(s.head().to_string()); print(f'\nwritten: {n.out}  ({len(s)} samples)')
