#!/usr/bin/env python3
from __future__ import annotations
import json, math, os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DP = ROOT/'data_processed'
FIG = ROOT/'figures'
DOC = ROOT/'docs'
FIG.mkdir(exist_ok=True)
DOC.mkdir(exist_ok=True)

PARAMS_STEP = ['shift','slide','rise','tilt','roll','twist']
PARAMS_PAIR = ['abs_stretch','abs_opening','abs_propeller']
ALLOWED_STANDARD = {'DA','DC','DG','DT','HOH','WAT','H2O','MG','NA','K','CA','CL','CO','CD','MN','ZN','NI','SR','BA','TL','NH4','BR'}


def clean_list(s):
    if pd.isna(s): return []
    return [x for x in str(s).split(';') if x and x!='nan']

def classify_pair_group(cls: str) -> str:
    s = str(cls).replace(' ', '').upper()
    if s in {'A:T','T:A'}: return 'AT_pair'
    if s in {'G:C','C:G'}: return 'GC_pair'
    return 'other'

def cliffs_delta(x, y):
    x=np.asarray(x); y=np.asarray(y)
    # positive means x > y
    if len(x)==0 or len(y)==0: return np.nan
    gt=0; lt=0
    for xi in x:
        gt += np.sum(xi > y)
        lt += np.sum(xi < y)
    return (gt-lt)/(len(x)*len(y))

def bootstrap_cluster_ci(df, value_col, group_col='step_category', cluster_col='pdb_id', n_boot=2000, seed=3409):
    rng=np.random.default_rng(seed)
    rows=[]
    for group, g in df.groupby(group_col):
        pids=np.array(sorted(g[cluster_col].unique()))
        vals=[]
        for _ in range(n_boot):
            sampled=rng.choice(pids, size=len(pids), replace=True)
            parts=[]
            for pid in sampled:
                parts.append(g[g[cluster_col]==pid][value_col].values)
            arr=np.concatenate(parts) if parts else np.array([])
            if len(arr): vals.append(np.mean(arr))
        vals=np.array(vals)
        rows.append({
            group_col: group,
            'mean': g[value_col].mean(),
            'cluster_boot_ci_low': float(np.quantile(vals,0.025)),
            'cluster_boot_ci_high': float(np.quantile(vals,0.975)),
            'n_steps': len(g),
            'n_structures': len(pids),
        })
    return pd.DataFrame(rows)

def zscore_rows(df, baseline, params, category_col='step_category'):
    rows=[]
    base_all={p:(baseline[p].mean(), baseline[p].std(ddof=1)) for p in params}
    base_by_cat={}
    for cat,g in baseline.groupby(category_col):
        base_by_cat[cat]={p:(g[p].mean(), g[p].std(ddof=1)) for p in params}
    for _,r in df.iterrows():
        cat=r.get(category_col)
        use_cat = cat if cat in base_by_cat and cat in ['AT-rich','GC-rich','mixed'] else 'all'
        d=r.to_dict(); ss=0.0
        for p in params:
            mu,sd=(base_by_cat[use_cat][p] if use_cat!='all' else base_all[p])
            if sd is None or not math.isfinite(sd) or sd==0 or pd.isna(r[p]):
                z=np.nan
            else:
                z=(r[p]-mu)/sd
                ss += z*z
            d[f'z_{p}']=z
        d['baseline_used_for_D_step']=use_cat
        d['D_step']=math.sqrt(ss) if ss>0 else np.nan
        rows.append(d)
    return pd.DataFrame(rows)


def main():
    normal = pd.read_csv(DP/'normal_steps_plausible.csv')
    screen = pd.read_csv(DP/'week2_screening_corrected.csv')
    lesion_near = pd.read_csv(DP/'lesion_steps_8og_near.csv')
    bp_all = pd.read_csv(DP/'base_pair_internal_all.csv')

    # Audit normal candidates: exclude nonstandard base / strong bound molecules from strict baseline.
    audit_rows=[]
    for _,r in screen[screen['bucket'].astype(str).str.startswith('normal')].iterrows():
        comps=clean_list(r.get('all_chem_comp_ids'))
        nonstandard=[c for c in comps if c not in ALLOWED_STANDARD]
        subset = normal.loc[normal['pdb_id'].eq(r['pdb_id']),'baseline_subset']
        subset = subset.iloc[0] if len(subset) else 'none'
        # strict: core or tier1 AND no nonstandard/strong ligand/mod base flags
        strict = subset in ['core_original','expanded_tier1'] and len(nonstandard)==0
        if strict:
            status='strict_baseline'
        elif subset in ['core_original','expanded_tier1'] and len(nonstandard)>0:
            status='sensitivity_or_exclude_nonstandard'
        elif subset=='sensitivity_tier2':
            status='sensitivity_only'
        elif subset=='screen_only':
            status='screen_only'
        else:
            status='exclude_or_manual_review'
        audit_rows.append({**r.to_dict(), 'parsed_baseline_subset':subset,
                           'nonstandard_comp_flags':';'.join(nonstandard),
                           'final_week2_status':status,
                           'strict_baseline_use': strict})
    audit=pd.DataFrame(audit_rows)
    strict_ids=set(audit.loc[audit['strict_baseline_use'],'pdb_id'])
    strict_normal=normal[normal['pdb_id'].isin(strict_ids) & normal['step_category'].isin(['AT-rich','GC-rich','mixed'])].copy()
    strict_normal.to_csv(DP/'week2_strict_normal_steps.csv', index=False)
    audit.to_csv(DP/'week2_baseline_screening_audit.csv', index=False)

    # Baseline summary and stats.
    baseline_summary=[]
    for cat,g in strict_normal.groupby('step_category'):
        rec={'step_category':cat,'n_steps':len(g),'n_structures':g['pdb_id'].nunique()}
        for p in ['shift','slide','rise','tilt','roll','twist','local_bend']:
            rec[f'{p}_mean']=g[p].mean(); rec[f'{p}_sd']=g[p].std(ddof=1); rec[f'{p}_median']=g[p].median()
        baseline_summary.append(rec)
    strict_summary=pd.DataFrame(baseline_summary)
    strict_summary.to_csv(DP/'week2_strict_baseline_summary.csv', index=False)

    # Nonterminal sensitivity.
    nt = strict_normal[~strict_normal['is_terminal_step']]
    nt_summary=[]
    for cat,g in nt.groupby('step_category'):
        nt_summary.append({'step_category':cat,'n_steps':len(g),'n_structures':g['pdb_id'].nunique(),
                           'local_bend_mean':g['local_bend'].mean(),'local_bend_sd':g['local_bend'].std(ddof=1),
                           'rise_mean':g['rise'].mean(),'twist_mean':g['twist'].mean()})
    pd.DataFrame(nt_summary).to_csv(DP/'week2_strict_nonterminal_sensitivity_summary.csv', index=False)

    # Tests and effect sizes for local_bend.
    groups={cat:g['local_bend'].dropna().values for cat,g in strict_normal.groupby('step_category')}
    kw = stats.kruskal(*(groups[c] for c in ['AT-rich','GC-rich','mixed'] if c in groups))
    gc=groups.get('GC-rich', np.array([])); at=groups.get('AT-rich', np.array([])); mixed=groups.get('mixed', np.array([]))
    def cohens_d(x,y):
        nx,ny=len(x),len(y)
        sx=np.var(x, ddof=1); sy=np.var(y, ddof=1)
        pooled=math.sqrt(((nx-1)*sx+(ny-1)*sy)/(nx+ny-2)) if nx+ny>2 else np.nan
        return (np.mean(x)-np.mean(y))/pooled if pooled else np.nan
    tests=pd.DataFrame([
        {'comparison':'Kruskal-Wallis AT/GC/mixed local_bend', 'statistic':kw.statistic, 'p_value':kw.pvalue, 'effect':'-', 'effect_value':np.nan},
        {'comparison':'GC-rich vs AT-rich local_bend', 'statistic':stats.mannwhitneyu(gc, at, alternative='two-sided').statistic, 'p_value':stats.mannwhitneyu(gc, at, alternative='two-sided').pvalue, 'effect':'Cohen_d_GC_minus_AT', 'effect_value':cohens_d(gc,at)},
        {'comparison':'GC-rich vs AT-rich local_bend', 'statistic':np.nan, 'p_value':np.nan, 'effect':'Cliff_delta_GC_minus_AT', 'effect_value':cliffs_delta(gc,at)},
        {'comparison':'GC-rich vs mixed local_bend', 'statistic':stats.mannwhitneyu(gc, mixed, alternative='two-sided').statistic, 'p_value':stats.mannwhitneyu(gc, mixed, alternative='two-sided').pvalue, 'effect':'Cohen_d_GC_minus_mixed', 'effect_value':cohens_d(gc,mixed)},
    ])
    tests.to_csv(DP/'week2_strict_baseline_tests.csv', index=False)
    boot=bootstrap_cluster_ci(strict_normal, 'local_bend')
    boot.to_csv(DP/'week2_structure_level_bootstrap_local_bend.csv', index=False)

    # Score 8OG-near steps against strict baseline.
    lesion_scored = zscore_rows(lesion_near, strict_normal, PARAMS_STEP)
    lesion_scored.to_csv(DP/'week2_lesion_8og_near_scored_strict.csv', index=False)
    lesion_plausible = lesion_scored[lesion_scored['geometry_quality']=='plausible'].copy()
    lesion_summary = lesion_plausible.groupby(['bucket','manifest_protein_bound']).agg(
        n_steps=('pdb_id','count'), mean_D_step=('D_step','mean'), median_D_step=('D_step','median'), max_D_step=('D_step','max')
    ).reset_index()
    lesion_summary.to_csv(DP/'week2_lesion_Dstep_summary_strict.csv', index=False)

    # Pair internal scoring.
    bp_all['pair_group']=bp_all['base_pair_class'].map(classify_pair_group)
    bp_all['abs_stretch']=bp_all['stretch'].abs()
    bp_all['abs_opening']=bp_all['opening'].abs()
    bp_all['abs_propeller']=bp_all['propeller'].abs()
    bp_baseline=bp_all[(bp_all['pdb_id'].isin(strict_ids)) & (bp_all['pair_group'].isin(['AT_pair','GC_pair'])) & (bp_all['contains_8OG'].astype(str).str.lower().eq('no'))].copy()
    pair_base=bp_baseline.groupby('pair_group')[PARAMS_PAIR].agg(['mean','std','count'])
    # flatten for export
    pair_summary=[]
    for group,g in bp_baseline.groupby('pair_group'):
        rec={'pair_group':group,'n_pairs':len(g),'n_structures':g['pdb_id'].nunique()}
        for p in PARAMS_PAIR:
            rec[f'{p}_mean']=g[p].mean(); rec[f'{p}_sd']=g[p].std(ddof=1); rec[f'{p}_median']=g[p].median()
        pair_summary.append(rec)
    pd.DataFrame(pair_summary).to_csv(DP/'week2_strict_pair_internal_baseline_summary.csv', index=False)

    bp8 = bp_all[bp_all['contains_8OG'].astype(str).str.lower().eq('yes')].copy()
    scored=[]
    base_stats={}
    for group,g in bp_baseline.groupby('pair_group'):
        base_stats[group]={p:(g[p].mean(), g[p].std(ddof=1)) for p in PARAMS_PAIR}
    for _,r in bp8.iterrows():
        use='AT_pair' if str(r['opposite_base_to_8OG']).upper()=='A' else ('GC_pair' if str(r['opposite_base_to_8OG']).upper()=='C' else 'unknown')
        d=r.to_dict(); ss=0.0
        for p in PARAMS_PAIR:
            if use in base_stats:
                mu,sd=base_stats[use][p]
                z=(r[p]-mu)/sd if sd and not pd.isna(sd) else np.nan
            else:
                z=np.nan
            d[f'z_{p}']=z
            if pd.notna(z): ss += z*z
        d['pair_baseline_used']=use
        d['D_pair']=math.sqrt(ss) if ss>0 else np.nan
        scored.append(d)
    bp8_scored=pd.DataFrame(scored)
    bp8_scored.to_csv(DP/'week2_8og_pair_internal_scored_strict.csv', index=False)

    # Figures.
    plt.figure(figsize=(7,4.5))
    labels=['기존 core','strict 확장','전체 normal 후보']
    counts=[55, len(strict_normal), len(normal[normal['geometry_quality'].eq('plausible')])]
    plt.bar(labels, counts)
    plt.ylabel('plausible base-pair steps')
    plt.title('정상 B-DNA baseline 규모 확장')
    for i,c in enumerate(counts): plt.text(i,c+3,str(c),ha='center')
    plt.tight_layout(); plt.savefig(FIG/'fig1_baseline_size_expansion.png', dpi=200); plt.close()

    order=['AT-rich','GC-rich','mixed']
    data=[strict_normal.loc[strict_normal.step_category.eq(c),'local_bend'] for c in order]
    plt.figure(figsize=(7,4.5))
    plt.boxplot(data, labels=order, showmeans=True)
    plt.ylabel('local bend (degree)')
    plt.title('Strict expanded baseline: local bend by step category')
    plt.tight_layout(); plt.savefig(FIG/'fig2_local_bend_by_category_strict.png', dpi=200); plt.close()

    plt.figure(figsize=(6,5))
    for cat,g in strict_normal.groupby('step_category'):
        plt.scatter(g['twist'], g['rise'], label=cat, alpha=0.75)
    plt.xlabel('twist (degree)'); plt.ylabel('rise (Å)')
    plt.title('Strict baseline: twist-rise distribution')
    plt.legend()
    plt.tight_layout(); plt.savefig(FIG/'fig3_twist_rise_strict.png', dpi=200); plt.close()

    top=lesion_scored.sort_values('D_step', ascending=False).head(15).copy()
    top['label']=top['pdb_id'].astype(str)+' step '+top['step_number'].astype(int).astype(str)+' '+top['relation_to_8OG_pair'].astype(str)
    plt.figure(figsize=(8,6))
    plt.barh(top['label'][::-1], top['D_step'][::-1])
    plt.xlabel('D_step')
    plt.title('8OG-near step distortion score (strict baseline)')
    plt.tight_layout(); plt.savefig(FIG/'fig4_D_step_8og_near_top15.png', dpi=220); plt.close()

    heat=lesion_scored[lesion_scored['geometry_quality'].eq('plausible')].sort_values('D_step', ascending=False).head(12)
    zcols=[f'z_{p}' for p in PARAMS_STEP]
    mat=heat[zcols].to_numpy(dtype=float)
    plt.figure(figsize=(8,5))
    im=plt.imshow(mat, aspect='auto')
    plt.colorbar(im, label='z-score')
    plt.xticks(range(len(zcols)), [c[2:] for c in zcols], rotation=45, ha='right')
    plt.yticks(range(len(heat)), heat['pdb_id'].astype(str)+' '+heat['step_number'].astype(int).astype(str))
    plt.title('8OG-near plausible steps: parameter-wise z-score')
    plt.tight_layout(); plt.savefig(FIG/'fig5_zscore_heatmap_8og_near.png', dpi=220); plt.close()

    # opening-stretch: baseline AT/GC and 8OG pairs.
    plt.figure(figsize=(7,5))
    for group,g in bp_baseline.groupby('pair_group'):
        plt.scatter(g['abs_stretch'], g['abs_opening'], alpha=0.35, label=group)
    for _,r in bp8_scored.iterrows():
        plt.scatter(r['abs_stretch'], r['abs_opening'], marker='x', s=80)
        plt.text(r['abs_stretch'], r['abs_opening'], f"{r['pdb_id']} {r['base_pair_class']}", fontsize=8)
    plt.xlabel('|stretch| (Å)'); plt.ylabel('|opening| (degree)')
    plt.title('Base-pair internal geometry: normal vs 8OG')
    plt.legend()
    plt.tight_layout(); plt.savefig(FIG/'fig6_opening_stretch_8og_pairs.png', dpi=220); plt.close()

    # D_step vs D_pair: pair D for 8OG direct, step score max/mean by PDB around 8OG.
    agg_step = lesion_scored[lesion_scored['relation_to_8OG_pair'].isin(['contains_8OG_pair','flanks_8OG_pair'])].groupby('pdb_id').agg(
        mean_D_step=('D_step','mean'), max_D_step=('D_step','max')
    ).reset_index()
    pair_direct=bp8_scored[['pdb_id','base_pair_class','opposite_base_to_8OG','D_pair','bucket','protein_bound_manifest']].copy()
    merged=pair_direct.merge(agg_step, on='pdb_id', how='left')
    merged.to_csv(DP/'week2_Dstep_Dpair_joined.csv', index=False)
    plt.figure(figsize=(6.5,5))
    for _,r in merged.iterrows():
        plt.scatter(r['mean_D_step'], r['D_pair'], marker='o')
        plt.text(r['mean_D_step'], r['D_pair'], f"{r['pdb_id']} {r['base_pair_class']}", fontsize=8)
    plt.xlabel('mean D_step near 8OG')
    plt.ylabel('D_pair for 8OG pair')
    plt.title('Step distortion vs pair internal distortion')
    plt.tight_layout(); plt.savefig(FIG/'fig7_Dstep_vs_Dpair.png', dpi=220); plt.close()

    key={
        'strict_baseline_n_steps': int(len(strict_normal)),
        'strict_baseline_n_structures': int(strict_normal['pdb_id'].nunique()),
        'original_core_n_steps': int((normal['baseline_subset'].eq('core_original') & normal['geometry_quality'].eq('plausible')).sum()),
        'all_normal_candidate_plausible_steps': int((normal['geometry_quality'].eq('plausible')).sum()),
        'excluded_from_strict_due_to_nonstandard_or_bound_compounds': audit.loc[audit['final_week2_status'].eq('sensitivity_or_exclude_nonstandard'), ['pdb_id','nonstandard_comp_flags']].to_dict('records'),
        'strict_local_bend_means': strict_summary[['step_category','n_steps','local_bend_mean','local_bend_sd']].to_dict('records'),
        'local_bend_kruskal_p': float(kw.pvalue),
        'local_bend_gc_vs_at_cohen_d': float(tests.loc[tests['effect'].eq('Cohen_d_GC_minus_AT'),'effect_value'].iloc[0]),
        'n_8og_near_steps_total': int(len(lesion_scored)),
        'n_8og_near_plausible_steps': int(len(lesion_plausible)),
        'top_8og_D_step_cases': lesion_scored.sort_values('D_step', ascending=False).head(8)[['pdb_id','bucket','step_number','step_seq','relation_to_8OG_pair','geometry_quality','D_step']].to_dict('records'),
        '8og_pair_internal_scored': bp8_scored[['pdb_id','bucket','protein_bound_manifest','base_pair_class','opposite_base_to_8OG','stretch','opening','propeller','D_pair']].to_dict('records')
    }
    (DP/'week2_key_results.json').write_text(json.dumps(key, indent=2, ensure_ascii=False), encoding='utf-8')

    md=f"""# Week 2 mmCIF parsing and baseline expansion result\n\n## Completed in session\n\n- Downloaded and parsed 43 RCSB mmCIF structures from the expanded Week 1 manifest.\n- Extracted `_ndb_struct_na_base_pair` and `_ndb_struct_na_base_pair_step` geometry tables.\n- Reproduced the original core baseline of 55 plausible normal B-DNA steps.\n- Built a stricter expanded baseline by excluding nonstandard/base-modified or strongly bound-compound candidates.\n\n## Final strict baseline\n\n- Original core baseline: **55 plausible steps**.\n- Strict expanded baseline: **{len(strict_normal)} plausible steps** from **{strict_normal['pdb_id'].nunique()} structures**.\n- All parsed normal candidates before strict filtering: **{int((normal['geometry_quality'].eq('plausible')).sum())} plausible steps**.\n\nThe strict baseline excludes candidates flagged by nonstandard components such as TAF, SPM, or NT from the main baseline and keeps them for sensitivity/manual review only.\n\n## Strict baseline local bend means\n\n{strict_summary[['step_category','n_steps','n_structures','local_bend_mean','local_bend_sd','rise_mean','twist_mean']].round(3).to_markdown(index=False)}\n\nKruskal-Wallis test for local bend across AT-rich, GC-rich, and mixed: p = {kw.pvalue:.3g}.\nGC-rich vs AT-rich Cohen's d = {tests.loc[tests['effect'].eq('Cohen_d_GC_minus_AT'),'effect_value'].iloc[0]:.3f}.\n\n## Week 2 output tables\n\n- `week2_baseline_screening_audit.csv`\n- `week2_strict_normal_steps.csv`\n- `week2_strict_baseline_summary.csv`\n- `week2_strict_baseline_tests.csv`\n- `week2_structure_level_bootstrap_local_bend.csv`\n- `week2_lesion_8og_near_scored_strict.csv`\n- `week2_8og_pair_internal_scored_strict.csv`\n- `week2_key_results.json`\n\n## Figures\n\n- `fig1_baseline_size_expansion.png`\n- `fig2_local_bend_by_category_strict.png`\n- `fig3_twist_rise_strict.png`\n- `fig4_D_step_8og_near_top15.png`\n- `fig5_zscore_heatmap_8og_near.png`\n- `fig6_opening_stretch_8og_pairs.png`\n- `fig7_Dstep_vs_Dpair.png`\n"""
    (DOC/'week2_final_progress_report.md').write_text(md, encoding='utf-8')
    print(json.dumps(key, indent=2, ensure_ascii=False))

if __name__ == '__main__':
    main()
