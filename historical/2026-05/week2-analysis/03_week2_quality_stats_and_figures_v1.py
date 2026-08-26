#!/usr/bin/env python3
from __future__ import annotations
import csv, math, random, statistics
from pathlib import Path
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
try:
    from scipy.stats import kruskal, mannwhitneyu
except Exception:
    kruskal = mannwhitneyu = None

ROOT=Path('.')
OUT=ROOT/'data_processed'
FIG=ROOT/'figures'
FIG.mkdir(exist_ok=True)

def read_csv(path):
    with open(path, newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))

def write_csv(path, rows, fields=None):
    if fields is None:
        fields=[]
        for r in rows:
            for k in r:
                if k not in fields: fields.append(k)
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

def fnum(x):
    try:
        return float(x)
    except Exception:
        return None

def mean_sd(vals):
    vals=[v for v in vals if v is not None and math.isfinite(v)]
    if not vals: return ('','','','')
    return (len(vals), statistics.mean(vals), statistics.stdev(vals) if len(vals)>1 else 0.0, statistics.median(vals))

def cohen_d(a,b):
    if len(a)<2 or len(b)<2: return None
    va=statistics.variance(a); vb=statistics.variance(b)
    pooled=math.sqrt(((len(a)-1)*va+(len(b)-1)*vb)/(len(a)+len(b)-2))
    return (statistics.mean(a)-statistics.mean(b))/pooled if pooled else None

def cliffs_delta(a,b):
    if not a or not b: return None
    gt=lt=0
    for x in a:
        for y in b:
            if x>y: gt+=1
            elif x<y: lt+=1
    return (gt-lt)/(len(a)*len(b))

steps=read_csv(OUT/'normal_steps.csv')
lesion=read_csv(OUT/'lesion_steps_scored.csv')
bp=read_csv(OUT/'base_pair_internal_scored.csv')
manifest=read_csv(OUT/'week1_expanded_pdb_manifest.csv')

# Data quality summaries
quality=[]
for label, subset in [
    ('core_original', [r for r in steps if r['baseline_subset']=='core_original']),
    ('expanded_main', [r for r in steps if r['baseline_subset'] in {'core_original','expanded_main'}]),
    ('expanded_main_nonterminal', [r for r in steps if r['baseline_subset'] in {'core_original','expanded_main'} and r['terminal_step']=='no']),
    ('all_normal_sensitivity', steps),
    ('lesion_all', lesion),
]:
    c=Counter(r['geometry_class'] for r in subset)
    quality.append({'dataset':label,'n_total':len(subset),'n_plausible':c.get('plausible',0),'n_outside':c.get('outside',0),'n_extreme':c.get('extreme',0),'unique_pdb':len(set(r['pdb_id'] for r in subset))})
write_csv(OUT/'week2_data_quality_summary.csv',quality,['dataset','n_total','n_plausible','n_outside','n_extreme','unique_pdb'])

# Kruskal and effect sizes for local bend by category
stats_rows=[]
for label, filt in [
    ('core_original', lambda r: r['baseline_subset']=='core_original' and r['geometry_class']=='plausible'),
    ('expanded_main', lambda r: r['baseline_subset'] in {'core_original','expanded_main'} and r['geometry_class']=='plausible'),
    ('expanded_main_nonterminal', lambda r: r['baseline_subset'] in {'core_original','expanded_main'} and r['geometry_class']=='plausible' and r['terminal_step']=='no'),
    ('all_normal_sensitivity', lambda r: r['geometry_class']=='plausible'),
]:
    rows=[r for r in steps if filt(r)]
    groups={cat:[fnum(r['local_bend']) for r in rows if r['step_category']==cat and fnum(r['local_bend']) is not None] for cat in ['AT-rich','GC-rich','mixed']}
    p=''
    H=''
    if kruskal and all(len(v)>0 for v in groups.values()):
        try:
            res=kruskal(*groups.values())
            H=f'{res.statistic:.6g}'; p=f'{res.pvalue:.6g}'
        except Exception: pass
    for cat, vals in groups.items():
        n,m,sd,med=mean_sd(vals)
        stats_rows.append({'baseline_subset':label,'analysis':'local_bend_by_category','group':cat,'n':n,'mean':f'{m:.6g}' if m!='' else '', 'sd':f'{sd:.6g}' if sd!='' else '', 'median':f'{med:.6g}' if med!='' else '', 'kruskal_H_all_categories':H,'kruskal_p_all_categories':p,'effect_vs':'','cohen_d':'','cliffs_delta':'','mannwhitney_p':''})
    pairs=[('GC-rich','AT-rich'),('mixed','AT-rich'),('GC-rich','mixed')]
    for a,b in pairs:
        av=groups[a]; bv=groups[b]
        mw=''
        if mannwhitneyu and av and bv:
            try: mw=f'{mannwhitneyu(av,bv,alternative="two-sided").pvalue:.6g}'
            except Exception: pass
        d=cohen_d(av,bv); cd=cliffs_delta(av,bv)
        stats_rows.append({'baseline_subset':label,'analysis':'local_bend_pairwise','group':a,'n':len(av),'mean':f'{statistics.mean(av):.6g}' if av else '', 'sd':f'{statistics.stdev(av):.6g}' if len(av)>1 else '', 'median':f'{statistics.median(av):.6g}' if av else '', 'kruskal_H_all_categories':H,'kruskal_p_all_categories':p,'effect_vs':b,'cohen_d':f'{d:.6g}' if d is not None else '', 'cliffs_delta':f'{cd:.6g}' if cd is not None else '', 'mannwhitney_p':mw})
write_csv(OUT/'week2_local_bend_stats.csv',stats_rows)

# Lesion summaries
lesion_summary=[]
for group_label, filt in [
    ('near_8OG_plausible_all_lesion', lambda r: r['near_8OG']=='yes' and r['geometry_class']=='plausible'),
    ('near_8OG_plausible_primary_free', lambda r: r['near_8OG']=='yes' and r['geometry_class']=='plausible' and r['bucket']=='8OG_lesion_primary'),
    ('near_8OG_plausible_repair_bound', lambda r: r['near_8OG']=='yes' and r['geometry_class']=='plausible' and 'repair' in r['bucket']),
    ('near_8OG_plausible_polymerase_bound', lambda r: r['near_8OG']=='yes' and r['geometry_class']=='plausible' and 'polymerase' in r['bucket']),
    ('near_8OG_outside_extreme', lambda r: r['near_8OG']=='yes' and r['geometry_class']!='plausible'),
]:
    rows=[r for r in lesion if filt(r) and fnum(r.get('D_step')) is not None]
    vals=[fnum(r['D_step']) for r in rows]
    n,m,sd,med=mean_sd(vals)
    lesion_summary.append({'group':group_label,'n':n,'mean_D_step':f'{m:.6g}' if m!='' else '', 'sd_D_step':f'{sd:.6g}' if sd!='' else '', 'median_D_step':f'{med:.6g}' if med!='' else '', 'min_D_step':f'{min(vals):.6g}' if vals else '', 'max_D_step':f'{max(vals):.6g}' if vals else '', 'unique_pdb':len(set(r['pdb_id'] for r in rows))})
write_csv(OUT/'week2_lesion_Dstep_summary.csv',lesion_summary)

# Top lesion rows plausible near 8OG and extreme near 8OG
plaus=[r for r in lesion if r['near_8OG']=='yes' and r['geometry_class']=='plausible' and fnum(r.get('D_step')) is not None]
plaus=sorted(plaus,key=lambda r:fnum(r['D_step']), reverse=True)[:25]
extreme=[r for r in lesion if r['near_8OG']=='yes' and r['geometry_class']!='plausible' and fnum(r.get('D_step')) is not None]
extreme=sorted(extreme,key=lambda r:fnum(r['D_step']), reverse=True)[:25]
top_fields=['pdb_id','bucket','step_number','step_seq','contains_8OG','flanks_8OG_sequence','geometry_class','D_step','shift','slide','rise','tilt','roll','twist','local_bend','z_shift','z_slide','z_rise','z_tilt','z_roll','z_twist']
write_csv(OUT/'week2_top_plausible_near8OG_Dstep.csv',plaus,top_fields)
write_csv(OUT/'week2_top_outside_extreme_near8OG_Dstep.csv',extreme,top_fields)

# 8OG direct base-pair summary
og_bp=[r for r in bp if r['contains_8OG']=='yes']
og_summary=[]
for label, filt in [
    ('8OG:A direct pairs', lambda r: r['opposite_base_to_8OG']=='A'),
    ('8OG:C direct pairs', lambda r: r['opposite_base_to_8OG']=='C'),
    ('8OG:A direct free only', lambda r: r['opposite_base_to_8OG']=='A' and r['bucket']=='8OG_lesion_primary'),
    ('8OG:C direct free only', lambda r: r['opposite_base_to_8OG']=='C' and r['bucket']=='8OG_lesion_primary'),
]:
    rows=[r for r in og_bp if filt(r)]
    rec={'group':label,'n':len(rows),'unique_pdb':len(set(r['pdb_id'] for r in rows))}
    for p in ['stretch','opening','propeller','D_pair']:
        vals=[abs(fnum(r[p])) if p!='D_pair' else fnum(r[p]) for r in rows if fnum(r[p]) is not None]
        n,m,sd,med=mean_sd(vals)
        rec[f'mean_abs_{p}' if p!='D_pair' else 'mean_D_pair']=f'{m:.6g}' if m!='' else ''
        rec[f'median_abs_{p}' if p!='D_pair' else 'median_D_pair']=f'{med:.6g}' if med!='' else ''
    og_summary.append(rec)
write_csv(OUT/'week2_8OG_direct_basepair_summary.csv',og_summary)
write_csv(OUT/'week2_8OG_direct_basepairs_scored.csv',og_bp,[k for k in og_bp[0].keys()] if og_bp else None)

# Figures
# 1 local bend boxplot expanded main
cats=['AT-rich','GC-rich','mixed']
fig, ax=plt.subplots(figsize=(7,4.5))
data=[]
for cat in cats:
    data.append([fnum(r['local_bend']) for r in steps if r['baseline_subset'] in {'core_original','expanded_main'} and r['geometry_class']=='plausible' and r['step_category']==cat and fnum(r['local_bend']) is not None])
ax.boxplot(data, tick_labels=cats, showmeans=True)
ax.set_title('Expanded normal B-DNA baseline: local bend by step category')
ax.set_ylabel('local bend = sqrt(tilt² + roll²) (deg)')
ax.grid(True, axis='y', alpha=0.3)
fig.tight_layout(); fig.savefig(FIG/'week2_expanded_local_bend_boxplot.png', dpi=180); plt.close(fig)

# 2 D_pair direct 8OG pairs scatter/bar style
fig, ax=plt.subplots(figsize=(7.5,4.5))
labels=[]; vals=[]
for r in og_bp:
    if fnum(r.get('D_pair')) is None: continue
    labels.append(f"{r['pdb_id']} {r['base_pair_class']}")
    vals.append(fnum(r['D_pair']))
ax.bar(range(len(vals)), vals)
ax.set_xticks(range(len(vals)))
ax.set_xticklabels(labels, rotation=45, ha='right')
ax.set_title('Direct 8OG base-pair internal distortion score')
ax.set_ylabel('D_pair')
ax.grid(True, axis='y', alpha=0.3)
fig.tight_layout(); fig.savefig(FIG/'week2_direct_8OG_Dpair_bar.png', dpi=180); plt.close(fig)

# 3 top plausible near 8OG D_step
fig, ax=plt.subplots(figsize=(8,5))
rows=plaus[:12]
labels=[f"{r['pdb_id']}:{r['step_number']} {r['step_seq']}" for r in reversed(rows)]
vals=[fnum(r['D_step']) for r in reversed(rows)]
ax.barh(range(len(vals)), vals)
ax.set_yticks(range(len(vals))); ax.set_yticklabels(labels)
ax.set_title('Top plausible near-8OG base-pair steps by D_step')
ax.set_xlabel('D_step')
ax.grid(True, axis='x', alpha=0.3)
fig.tight_layout(); fig.savefig(FIG/'week2_top_plausible_near8OG_Dstep.png', dpi=180); plt.close(fig)

# 4 D_step vs D_pair overview for direct 8OG structures by PDB
# match D_pair by pdb to max direct D_pair; D_step by max plausible near 8OG in same pdb
pdb_Dpair=defaultdict(float)
for r in og_bp:
    v=fnum(r.get('D_pair'))
    if v is not None: pdb_Dpair[r['pdb_id']]=max(pdb_Dpair[r['pdb_id']],v)
pdb_Dstep=defaultdict(float)
for r in lesion:
    if r['near_8OG']=='yes' and fnum(r.get('D_step')) is not None and r['geometry_class']=='plausible':
        pdb_Dstep[r['pdb_id']]=max(pdb_Dstep[r['pdb_id']],fnum(r['D_step']))
fig, ax=plt.subplots(figsize=(6,4.5))
for pdbid in sorted(set(pdb_Dpair)|set(pdb_Dstep)):
    if pdbid not in pdb_Dpair or pdbid not in pdb_Dstep: continue
    ax.scatter([pdb_Dstep[pdbid]],[pdb_Dpair[pdbid]])
    ax.annotate(pdbid,(pdb_Dstep[pdbid],pdb_Dpair[pdbid]),xytext=(3,3),textcoords='offset points',fontsize=8)
ax.set_title('Step-level vs base-pair internal distortion')
ax.set_xlabel('max plausible near-8OG D_step')
ax.set_ylabel('max direct 8OG D_pair')
ax.grid(True, alpha=0.3)
fig.tight_layout(); fig.savefig(FIG/'week2_Dstep_vs_Dpair_overview.png', dpi=180); plt.close(fig)

# Short markdown report
md=[]
md.append('# Week 2 parsing and first-pass analysis log\n')
md.append('## Completed\n')
for q in quality:
    md.append(f"- {q['dataset']}: total {q['n_total']}, plausible {q['n_plausible']}, outside {q['n_outside']}, extreme {q['n_extreme']}, PDB {q['unique_pdb']}\n")
md.append('\n## Key baseline result\n')
for r in stats_rows:
    if r['baseline_subset']=='expanded_main' and r['analysis']=='local_bend_by_category':
        md.append(f"- {r['group']}: n={r['n']}, mean local bend={r['mean']} deg, SD={r['sd']}.\n")
kw=[r for r in stats_rows if r['baseline_subset']=='expanded_main' and r['analysis']=='local_bend_by_category' and r['group']=='AT-rich']
if kw:
    md.append(f"- Expanded main Kruskal-Wallis p-value for local bend categories: {kw[0]['kruskal_p_all_categories']}.\n")
md.append('\n## Key 8OG result\n')
for r in og_summary:
    md.append(f"- {r['group']}: n={r['n']}, mean D_pair={r.get('mean_D_pair','')}, mean abs(opening)={r.get('mean_abs_opening','')}, mean abs(stretch)={r.get('mean_abs_stretch','')}.\n")
md.append('\n## Caveat\n')
md.append('- The first week screening script was intentionally conservative and falsely rejected normal B-DNA entries when generic text contained words like protein/enzyme. For week 2, inclusion is based on the manifest category plus direct NDB parsing from each mmCIF file.\n')
(ROOT/'week2_processing_log.md').write_text(''.join(md),encoding='utf-8')
print('Wrote week2 stats, figures, log')
