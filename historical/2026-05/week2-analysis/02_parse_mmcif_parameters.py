#!/usr/bin/env python3
from __future__ import annotations
import math
import json
from pathlib import Path
import pandas as pd
import numpy as np
import gemmi

ROOT = Path('.')
MANIFEST = ROOT / 'data_processed' / 'week1_expanded_pdb_manifest.csv'
MMCIF_DIR = ROOT / 'data_raw' / 'mmcif'
OUT_DIR = ROOT / 'data_processed'
OUT_DIR.mkdir(exist_ok=True, parents=True)

STD_BASE = {'DA':'A','A':'A','ADE':'A','DG':'G','G':'G','GUA':'G','DC':'C','C':'C','CYT':'C','DT':'T','T':'T','THY':'T'}
LESION_COMPS = {'8OG','O8G','8OX','OG','A1C6T'}  # A1C6T appears in some product-bound repair entries; not always in base-pair table
DNA_COMP_PREFIXES = {'DA','DG','DC','DT','A','G','C','T','8OG','O8G','A1C6T'}
COMMON_NONPOLY = {'HOH','DOD','WAT','MG','CA','K','NA','CL','ZN','MN','CO','TL','SPM','SO4','GOL','EDO','MOO'}
STRONG_LIGANDS = {'NDP','DTP','0KX','AAB'}  # keep as flags; not automatic baseline rejects unless normal

def clean(s):
    if s is None:
        return ''
    s = str(s)
    if s in ['?','.']:
        return ''
    return s.strip().strip("'").strip('"')

def to_float(x):
    x = clean(x)
    if x == '': return np.nan
    try: return float(x)
    except Exception: return np.nan

def to_int(x):
    x = clean(x)
    if x == '': return np.nan
    try: return int(float(x))
    except Exception: return np.nan

def table_df(block, prefix: str) -> pd.DataFrame:
    try:
        tab = block.find_mmcif_category(prefix)
    except Exception:
        return pd.DataFrame()
    if tab is None or len(tab) == 0:
        return pd.DataFrame()
    cols = [tag.split('.')[-1] for tag in tab.tags]
    rows = []
    for r in tab:
        rows.append([clean(r[i]) for i in range(len(cols))])
    return pd.DataFrame(rows, columns=cols)

def base_letter(comp: str) -> str:
    comp = clean(comp).upper()
    return STD_BASE.get(comp, 'O' if comp in LESION_COMPS else comp)

def is_lesion(comp: str) -> bool:
    return clean(comp).upper() in LESION_COMPS

def step_category(seq: str) -> str:
    if not isinstance(seq, str) or len(seq) != 2:
        return 'other'
    if any(x not in 'ATGC' for x in seq):
        return 'lesion/modified'
    if seq[0] in 'AT' and seq[1] in 'AT':
        return 'AT-rich'
    if seq[0] in 'GC' and seq[1] in 'GC':
        return 'GC-rich'
    return 'mixed'

def bp_class(i_comp, j_comp):
    i = base_letter(i_comp); j = base_letter(j_comp)
    comps = {i,j}
    if 'O' in comps:
        opp = (comps - {'O'})
        opp = next(iter(opp)) if opp else ''
        if opp == 'A': return '8OG:A'
        if opp == 'C': return '8OG:C'
        return '8OG:other'
    if comps == {'A','T'}: return 'normal_A:T/T:A'
    if comps == {'G','C'}: return 'normal_G:C/C:G'
    return 'other/mismatch'

def structure_info(block, pdb_id: str, manifest_row: dict) -> dict:
    ent = table_df(block, '_entity.')
    entpoly = table_df(block, '_entity_poly.')
    nonpoly = table_df(block, '_pdbx_entity_nonpoly.')
    chem = table_df(block, '_chem_comp.')
    exptl = table_df(block, '_exptl.')
    refine = table_df(block, '_refine.')

    poly_types = [] if entpoly.empty else entpoly.get('type', pd.Series(dtype=str)).map(clean).tolist()
    has_protein = any('polypeptide' in t.lower() or 'protein' in t.lower() for t in poly_types)
    dna_poly_count = sum(('deoxyribo' in t.lower()) or ('polydeoxyribonucleotide' in t.lower()) or ('dna' in t.lower()) for t in poly_types)
    nonpoly_ids = [] if nonpoly.empty else nonpoly.get('comp_id', pd.Series(dtype=str)).map(lambda x: clean(x).upper()).tolist()
    nonpoly_nonwater = sorted([x for x in set(nonpoly_ids) if x not in {'HOH','DOD','WAT'}])
    ion_or_solvent_flags = sorted([x for x in set(nonpoly_ids) if x in COMMON_NONPOLY and x not in {'HOH','DOD','WAT'}])
    strong_ligand_flags = sorted([x for x in set(nonpoly_ids) if x not in COMMON_NONPOLY])
    chem_ids = [] if chem.empty else chem.get('id', pd.Series(dtype=str)).map(lambda x: clean(x).upper()).tolist()
    all_comps = sorted(set(chem_ids))
    contains_8og = any(x in LESION_COMPS for x in all_comps)
    # method/resolution
    method = ''
    if not exptl.empty and 'method' in exptl.columns:
        method = clean(exptl['method'].iloc[0])
    res = np.nan
    for c in ['ls_d_res_high','pdbx_diffrn_id']:
        if not refine.empty and c in refine.columns:
            res = to_float(refine[c].iloc[0]); break
    if np.isnan(res):
        # try reflns or entry in manifest
        res = to_float(manifest_row.get('resolution_A',''))
    return {
        'pdb_id': pdb_id,
        'poly_types': ';'.join(poly_types),
        'has_protein_by_entity_poly': has_protein,
        'dna_poly_count': dna_poly_count,
        'nonpoly_comp_ids': ';'.join(nonpoly_nonwater),
        'ion_or_solvent_flags': ';'.join(ion_or_solvent_flags),
        'strong_ligand_flags': ';'.join(strong_ligand_flags),
        'contains_8og_comp': contains_8og,
        'method_parsed': method,
        'resolution_A_parsed': res,
        'all_chem_comp_ids': ';'.join(all_comps),
    }

def geometry_quality(row) -> str:
    vals = [row.get(k, np.nan) for k in ['shift','slide','rise','tilt','roll','twist']]
    if any(pd.isna(v) for v in vals):
        return 'missing'
    shift, slide, rise, tilt, roll, twist = vals
    if rise < 0.5 or rise > 6.0 or twist < 0 or twist > 80 or abs(tilt) > 60 or abs(roll) > 60 or abs(shift) > 8 or abs(slide) > 8:
        return 'extreme'
    if not (2.0 <= rise <= 4.2 and 15 <= twist <= 55 and abs(tilt) <= 30 and abs(roll) <= 30 and abs(shift) <= 4 and abs(slide) <= 4):
        return 'outside'
    return 'plausible'

def collect_rows():
    manifest = pd.read_csv(MANIFEST).fillna('')
    manifest['pdb_id'] = manifest['pdb_id'].astype(str).str.upper()
    info_rows=[]; bp_rows=[]; step_rows=[]
    for _, mrow in manifest.iterrows():
        pdb_id = mrow['pdb_id']
        cif_path = MMCIF_DIR / f'{pdb_id}.cif'
        if not cif_path.exists():
            info_rows.append({'pdb_id':pdb_id,'file_exists':False, **mrow.to_dict()})
            continue
        try:
            block = gemmi.cif.read_file(str(cif_path))[0]
        except Exception as e:
            info_rows.append({'pdb_id':pdb_id,'file_exists':True,'parse_error':str(e), **mrow.to_dict()})
            continue
        sinfo = structure_info(block, pdb_id, mrow.to_dict())
        sinfo.update(mrow.to_dict())
        sinfo['file_exists'] = True
        info_rows.append(sinfo)

        bp = table_df(block, '_ndb_struct_na_base_pair.')
        if not bp.empty:
            for num_col in ['shear','stretch','stagger','buckle','propeller','opening']:
                if num_col in bp.columns: bp[num_col] = bp[num_col].map(to_float)
            for num_col in ['i_label_seq_id','j_label_seq_id','pair_number']:
                if num_col in bp.columns: bp[num_col] = bp[num_col].map(to_int)
            bp['pdb_id'] = pdb_id
            bp['bucket'] = mrow['bucket']
            bp['priority'] = mrow['priority']
            bp['structure_type'] = 'normal' if str(mrow['bucket']).startswith('normal') else 'lesion_or_support'
            bp['manifest_protein_bound'] = mrow['protein_bound']
            bp['has_protein_by_entity_poly'] = sinfo['has_protein_by_entity_poly']
            bp['i_base'] = bp['i_label_comp_id'].map(base_letter) if 'i_label_comp_id' in bp.columns else ''
            bp['j_base'] = bp['j_label_comp_id'].map(base_letter) if 'j_label_comp_id' in bp.columns else ''
            bp['pair_class'] = bp.apply(lambda r: bp_class(r.get('i_label_comp_id',''), r.get('j_label_comp_id','')), axis=1)
            bp['contains_8OG'] = bp.apply(lambda r: is_lesion(r.get('i_label_comp_id','')) or is_lesion(r.get('j_label_comp_id','')), axis=1)
            bp['opposite_base'] = bp.apply(lambda r: r['j_base'] if is_lesion(r.get('i_label_comp_id','')) else (r['i_base'] if is_lesion(r.get('j_label_comp_id','')) else ''), axis=1)
            bp_rows.append(bp)

        st = table_df(block, '_ndb_struct_na_base_pair_step.')
        if not st.empty:
            for num_col in ['shift','slide','rise','tilt','roll','twist','x_displacement','y_displacement','helical_rise','inclination','tip','helical_twist']:
                if num_col in st.columns: st[num_col] = st[num_col].map(to_float)
            for num_col in ['i_label_seq_id_1','j_label_seq_id_1','i_label_seq_id_2','j_label_seq_id_2','step_number']:
                if num_col in st.columns: st[num_col] = st[num_col].map(to_int)
            st['pdb_id'] = pdb_id
            st['bucket'] = mrow['bucket']
            st['priority'] = mrow['priority']
            st['structure_type'] = 'normal' if str(mrow['bucket']).startswith('normal') else 'lesion_or_support'
            st['manifest_protein_bound'] = mrow['protein_bound']
            st['has_protein_by_entity_poly'] = sinfo['has_protein_by_entity_poly']
            st['nonpoly_comp_ids'] = sinfo['nonpoly_comp_ids']
            st['ion_or_solvent_flags'] = sinfo['ion_or_solvent_flags']
            st['strong_ligand_flags'] = sinfo['strong_ligand_flags']
            st['resolution_A_parsed'] = sinfo['resolution_A_parsed']
            # top strand sequence from i base of pair 1 and pair 2
            st['base1_top'] = st['i_label_comp_id_1'].map(base_letter) if 'i_label_comp_id_1' in st.columns else ''
            st['base2_top'] = st['i_label_comp_id_2'].map(base_letter) if 'i_label_comp_id_2' in st.columns else ''
            st['step_seq'] = st['base1_top'].fillna('') + st['base2_top'].fillna('')
            st['step_category'] = st['step_seq'].map(step_category)
            st['contains_8OG'] = st.apply(lambda r: any(is_lesion(r.get(c,'')) for c in ['i_label_comp_id_1','j_label_comp_id_1','i_label_comp_id_2','j_label_comp_id_2']), axis=1)
            # relation to lesion pair numbers: identify 8OG pair numbers from bp table if available
            lesion_pair_nums = []
            if not bp.empty and 'pair_number' in bp.columns:
                lesion_pair_nums = [int(x) for x in bp.loc[bp['contains_8OG'],'pair_number'].dropna().tolist()]
            def rel_to_8og(step_no):
                if pd.isna(step_no) or not lesion_pair_nums:
                    return ''
                s=int(step_no)
                # step s connects pair s and s+1 in most NDB outputs
                if any(k in {s, s+1} for k in lesion_pair_nums):
                    return 'contains_8OG_pair'
                if any(k in {s-1, s+2} for k in lesion_pair_nums):
                    return 'flanks_8OG_pair'
                return 'other_step'
            st['relation_to_8OG_pair'] = st['step_number'].map(rel_to_8og) if 'step_number' in st.columns else ''
            st['local_bend'] = np.sqrt(st['tilt']**2 + st['roll']**2)
            st['geometry_quality'] = st.apply(geometry_quality, axis=1)
            # terminal flag
            max_step = int(st['step_number'].max()) if 'step_number' in st.columns and not st['step_number'].isna().all() else np.nan
            st['is_terminal_step'] = st['step_number'].isin([1, max_step]) if not pd.isna(max_step) else False
            # baseline bucket labels
            if str(mrow['bucket']).startswith('normal'):
                if 'Core' in str(mrow['priority']):
                    baseline_subset = 'core_original'
                elif 'Tier 1' in str(mrow['priority']):
                    baseline_subset = 'expanded_tier1'
                elif 'Tier 2' in str(mrow['priority']):
                    baseline_subset = 'sensitivity_tier2'
                else:
                    baseline_subset = 'screen_only'
            else:
                baseline_subset = ''
            st['baseline_subset'] = baseline_subset
            step_rows.append(st)
    info = pd.DataFrame(info_rows)
    bp_all = pd.concat(bp_rows, ignore_index=True) if bp_rows else pd.DataFrame()
    st_all = pd.concat(step_rows, ignore_index=True) if step_rows else pd.DataFrame()
    return manifest, info, bp_all, st_all

def corrected_screen(info: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for _, r in info.iterrows():
        bucket = r.get('bucket','')
        is_normal = str(bucket).startswith('normal')
        file_exists = bool(r.get('file_exists', False))
        has_bp = False; has_step=False
        if file_exists:
            try:
                block=gemmi.cif.read_file(str(MMCIF_DIR / f"{r['pdb_id']}.cif"))[0]
                has_bp=len(block.find_mmcif_category('_ndb_struct_na_base_pair.'))>0
                has_step=len(block.find_mmcif_category('_ndb_struct_na_base_pair_step.'))>0
            except Exception:
                pass
        reason=[]
        if not file_exists: decision='file_missing'; reason.append('mmCIF not downloaded')
        elif not has_bp or not has_step: decision='manual_review'; reason.append('NDB base pair/step category missing')
        elif is_normal and bool(r.get('has_protein_by_entity_poly', False)):
            decision='reject_baseline'; reason.append('entity_poly contains polypeptide')
        elif is_normal and bool(r.get('contains_8og_comp', False)):
            decision='reject_baseline'; reason.append('normal candidate contains 8OG/modification flag')
        elif is_normal and str(r.get('strong_ligand_flags','')):
            decision='manual_review'; reason.append('non-common ligand present: '+str(r.get('strong_ligand_flags','')))
        elif is_normal and 'Screen' in str(r.get('priority','')):
            decision='candidate_screen_only'; reason.append('manifest screen-only')
        elif is_normal and 'Tier 2' in str(r.get('priority','')):
            decision='candidate_sensitivity'; reason.append('tier-2 cation/condition sensitivity subset')
        elif is_normal:
            decision='candidate_expanded_baseline'; reason.append('DNA-only/polydeoxyribonucleotide; NDB step category present')
        else:
            decision='lesion_or_support'; reason.append('not used as normal baseline; bucket-specific analysis')
        d = r.to_dict()
        d.update({'has_ndb_base_pair':has_bp,'has_ndb_base_pair_step':has_step,'screen_decision_corrected':decision,'screen_reason_corrected':'; '.join(reason)})
        rows.append(d)
    return pd.DataFrame(rows)

def summarize_baselines(steps: pd.DataFrame) -> pd.DataFrame:
    normal = steps[(steps['structure_type']=='normal') & (steps['geometry_quality']=='plausible')].copy()
    rows=[]
    versions = {
        'core_original': normal[normal['baseline_subset']=='core_original'],
        'expanded_core_plus_tier1': normal[normal['baseline_subset'].isin(['core_original','expanded_tier1'])],
        'expanded_core_plus_tier1_nonterminal': normal[normal['baseline_subset'].isin(['core_original','expanded_tier1']) & (~normal['is_terminal_step'])],
        'all_normal_candidates': normal[normal['baseline_subset'].isin(['core_original','expanded_tier1','sensitivity_tier2','screen_only'])],
        'sensitivity_tier2_only': normal[normal['baseline_subset']=='sensitivity_tier2'],
    }
    params = ['shift','slide','rise','tilt','roll','twist','local_bend']
    for ver, df in versions.items():
        if df.empty: continue
        for cat, g in df.groupby('step_category'):
            rec={'baseline_version':ver,'step_category':cat,'n_steps':len(g),'n_structures':g['pdb_id'].nunique()}
            for p in params:
                rec[f'{p}_mean']=g[p].mean()
                rec[f'{p}_sd']=g[p].std(ddof=1)
                rec[f'{p}_median']=g[p].median()
            rows.append(rec)
        # all categories total
        rec={'baseline_version':ver,'step_category':'all','n_steps':len(df),'n_structures':df['pdb_id'].nunique()}
        for p in params:
            rec[f'{p}_mean']=df[p].mean(); rec[f'{p}_sd']=df[p].std(ddof=1); rec[f'{p}_median']=df[p].median()
        rows.append(rec)
    return pd.DataFrame(rows)

def main():
    manifest, info, bp, steps = collect_rows()
    screen = corrected_screen(info)
    # output raw/parsing tables
    info.to_csv(OUT_DIR/'week2_structure_info.csv', index=False)
    screen.to_csv(OUT_DIR/'week2_screening_corrected.csv', index=False)
    bp.to_csv(OUT_DIR/'base_pair_internal_all.csv', index=False)
    steps.to_csv(OUT_DIR/'base_pair_steps_all.csv', index=False)
    normal_steps = steps[steps['structure_type']=='normal'].copy()
    lesion_steps = steps[steps['structure_type']!='normal'].copy()
    normal_steps.to_csv(OUT_DIR/'normal_steps_raw.csv', index=False)
    normal_steps[normal_steps['geometry_quality']=='plausible'].to_csv(OUT_DIR/'normal_steps_plausible.csv', index=False)
    lesion_steps.to_csv(OUT_DIR/'lesion_steps_raw.csv', index=False)
    lesion_near = lesion_steps[lesion_steps['relation_to_8OG_pair'].isin(['contains_8OG_pair','flanks_8OG_pair']) | lesion_steps['contains_8OG']].copy()
    lesion_near.to_csv(OUT_DIR/'lesion_steps_8og_near.csv', index=False)
    bp.to_csv(OUT_DIR/'base_pair_internal.csv', index=False)
    bp_8og = bp[bp['contains_8OG']].copy()
    bp_8og.to_csv(OUT_DIR/'base_pair_internal_8og.csv', index=False)
    summary = summarize_baselines(steps)
    summary.to_csv(OUT_DIR/'baseline_summary.csv', index=False)
    # manifest counts
    report = {
        'n_mmcif_files': int(info['file_exists'].sum()),
        'n_structures_total': int(len(info)),
        'n_normal_candidates': int(info['bucket'].astype(str).str.startswith('normal').sum()),
        'n_lesion_support': int((~info['bucket'].astype(str).str.startswith('normal')).sum()),
        'n_base_pair_rows': int(len(bp)),
        'n_step_rows_total': int(len(steps)),
        'n_normal_steps_raw': int(len(normal_steps)),
        'n_normal_steps_plausible': int((normal_steps['geometry_quality']=='plausible').sum()),
        'n_lesion_steps_raw': int(len(lesion_steps)),
        'n_lesion_steps_8og_near': int(len(lesion_near)),
        'baseline_counts_by_subset': normal_steps[normal_steps['geometry_quality']=='plausible'].groupby('baseline_subset')['pdb_id'].count().to_dict(),
        'baseline_structures_by_subset': normal_steps[normal_steps['geometry_quality']=='plausible'].groupby('baseline_subset')['pdb_id'].nunique().to_dict(),
        'screen_decision_counts': screen['screen_decision_corrected'].value_counts().to_dict(),
        'geometry_quality_normal': normal_steps['geometry_quality'].value_counts().to_dict(),
        'geometry_quality_lesion': lesion_steps['geometry_quality'].value_counts().to_dict(),
    }
    (OUT_DIR/'week2_parse_summary.json').write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(report, indent=2, ensure_ascii=False))

if __name__ == '__main__':
    main()
