#!/usr/bin/env python3
"""Week 2 parser: extract NDB nucleic-acid geometry from expanded mmCIF set.

Outputs:
  data_processed/base_pair_steps_all.csv
  data_processed/base_pair_internal_all.csv
  data_processed/normal_steps.csv
  data_processed/lesion_steps.csv
  data_processed/base_pair_internal.csv
  data_processed/baseline_summary.csv
  data_processed/week2_parse_log.csv
"""
from __future__ import annotations

import csv
import math
import re
import shlex
from collections import defaultdict, Counter
from pathlib import Path
from statistics import mean, pstdev, stdev
from typing import Any

ROOT = Path('.')
MANIFEST = ROOT / 'data_processed' / 'week1_expanded_pdb_manifest.csv'
MMCIF_DIR = ROOT / 'data_raw' / 'mmcif'
OUT_DIR = ROOT / 'data_processed'

MISSING = {'?', '.', ''}
EIGHT_OG_ALIASES = {'8OG', 'O8G', 'OG', '8OHG', 'OHG'}
BASE_MAP = {
    'DA': 'A', 'A': 'A', 'ADE': 'A',
    'DT': 'T', 'T': 'T', 'THY': 'T',
    'DC': 'C', 'C': 'C', 'CYT': 'C',
    'DG': 'G', 'G': 'G', 'GUA': 'G',
    '8OG': 'O', 'O8G': 'O', 'OG': 'O', '8OHG': 'O', 'OHG': 'O',
    'DU': 'U', 'U': 'U', 'URA': 'U',
}
NUM_FIELDS_STEP = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist', 'x_displacement', 'y_displacement', 'helical_rise', 'inclination', 'tip', 'helical_twist']
NUM_FIELDS_BP = ['shear', 'stretch', 'stagger', 'buckle', 'propeller', 'opening']


def tokenize_mmcif_line(line: str) -> list[str]:
    # mmCIF in these NDB loops uses simple quoted/unquoted tokens. shlex handles quotes.
    lx = line.strip()
    if not lx or lx.startswith('#'):
        return []
    try:
        return shlex.split(lx, posix=True)
    except ValueError:
        # Conservative fallback: split whitespace.
        return lx.split()


def parse_loop_category(path: Path, category: str) -> list[dict[str, str]]:
    """Parse a loop_ block for one category. Handles rows split over multiple lines."""
    lines = path.read_text(encoding='utf-8', errors='ignore').splitlines()
    rows: list[dict[str, str]] = []
    i = 0
    prefix = f'_{category}.'
    n = len(lines)
    while i < n:
        if lines[i].strip() != 'loop_':
            i += 1
            continue
        j = i + 1
        headers: list[str] = []
        # collect headers if they start with our category
        while j < n and lines[j].strip().startswith('_'):
            s = lines[j].strip()
            if s.startswith(prefix):
                headers.append(s[len(prefix):].split()[0])
                j += 1
            else:
                break
        if not headers:
            i += 1
            continue
        # If headers are mixed and broke before row start, skip this category block if no rows.
        values: list[str] = []
        # read data lines until #, new loop, or scalar data item.
        while j < n:
            s = lines[j].strip()
            if not s:
                j += 1
                continue
            if s == '#' or s == 'loop_' or (s.startswith('_') and not s.startswith(prefix)):
                break
            values.extend(tokenize_mmcif_line(s))
            while len(values) >= len(headers):
                row_values = values[:len(headers)]
                values = values[len(headers):]
                rows.append(dict(zip(headers, row_values)))
            j += 1
        i = j + 1
    return rows


def as_float(x: Any) -> float | None:
    if x is None:
        return None
    s = str(x).strip()
    if s in MISSING:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def clean_comp(comp: str | None) -> str:
    if not comp:
        return ''
    return str(comp).upper().strip().replace('"','').replace("'", '')


def base_letter(comp: str | None) -> str:
    c = clean_comp(comp)
    return BASE_MAP.get(c, c.replace('D','',1)[:1] if c.startswith('D') and len(c) > 1 else c[:1])


def contains_8og_comps(comps: list[str]) -> bool:
    return any(clean_comp(c) in EIGHT_OG_ALIASES or '8OG' in clean_comp(c) or 'O8G' in clean_comp(c) for c in comps)


def step_seq_from_row(row: dict[str, str]) -> str:
    a1 = base_letter(row.get('i_label_comp_id_1'))
    a2 = base_letter(row.get('i_label_comp_id_2'))
    if a1 and a2:
        return a1 + a2
    # fallback from step_name: AA_DC1DG2:DC23DG24_BB -> CD? not ideal
    return ''


def step_category(seq: str) -> str:
    if len(seq) < 2:
        return 'unknown'
    letters = [x for x in seq[:2]]
    if any(x not in {'A','T','G','C'} for x in letters):
        return 'modified_or_other'
    if all(x in {'A','T'} for x in letters):
        return 'AT-rich'
    if all(x in {'G','C'} for x in letters):
        return 'GC-rich'
    return 'mixed'


def is_terminal_step(row: dict[str, str]) -> bool:
    num = as_float(row.get('step_number'))
    if num is None:
        return False
    # Terminal detection can be refined by max per PDB later; initially first step marked terminal,
    # and the final step is marked after grouping.
    return int(num) == 1


def geometry_class(row: dict[str, Any]) -> str:
    vals = {k: as_float(row.get(k)) for k in ['rise','twist','roll','tilt','shift','slide']}
    if any(vals[k] is None or not math.isfinite(vals[k]) for k in vals):
        return 'missing'
    rise = vals['rise']; twist = vals['twist']; roll = vals['roll']; tilt = vals['tilt']; shift = vals['shift']; slide = vals['slide']
    assert rise is not None and twist is not None and roll is not None and tilt is not None and shift is not None and slide is not None
    # Hard extreme: frame discontinuity or impossible B-DNA-like geometry.
    if rise <= 0.5 or rise >= 6.5 or twist <= -5 or twist >= 80 or abs(roll) >= 45 or abs(tilt) >= 45 or abs(shift) >= 5 or abs(slide) >= 5:
        return 'extreme'
    # Plausible B-DNA-like range for baseline / primary analysis.
    if 2.4 <= rise <= 4.3 and 15 <= twist <= 55 and abs(roll) <= 25 and abs(tilt) <= 25 and abs(shift) <= 2.5 and abs(slide) <= 3.0:
        return 'plausible'
    return 'outside'


def baseline_subset(row: dict[str, str]) -> str:
    if not row.get('bucket','').startswith('normal'):
        return 'not_baseline'
    pr = row.get('priority','').lower()
    if 'screen' in pr:
        return 'screen_only'
    if 'sensitivity' in pr or 'tier 2' in pr:
        return 'sensitivity_only'
    if row.get('bucket') == 'normal_baseline_core':
        return 'core_original'
    return 'expanded_main'


def manifest_rows() -> dict[str, dict[str, str]]:
    with MANIFEST.open(newline='', encoding='utf-8-sig') as f:
        return {r['pdb_id']: r for r in csv.DictReader(f)}


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        keys = []
        for r in rows:
            for k in r:
                if k not in keys:
                    keys.append(k)
        fields = keys
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    manifest = manifest_rows()
    all_step_rows: list[dict[str, Any]] = []
    all_bp_rows: list[dict[str, Any]] = []
    parse_log: list[dict[str, Any]] = []

    for pdb_id, m in manifest.items():
        cif_path = MMCIF_DIR / f'{pdb_id}.cif'
        if not cif_path.exists():
            parse_log.append({'pdb_id': pdb_id, 'status': 'missing_file', 'base_pair_rows': 0, 'step_rows': 0, 'notes': ''})
            continue
        bp = parse_loop_category(cif_path, 'ndb_struct_na_base_pair')
        st = parse_loop_category(cif_path, 'ndb_struct_na_base_pair_step')
        parse_log.append({'pdb_id': pdb_id, 'status': 'parsed', 'base_pair_rows': len(bp), 'step_rows': len(st), 'notes': ''})

        for r in st:
            rr: dict[str, Any] = dict(r)
            rr['pdb_id'] = pdb_id
            rr['bucket'] = m.get('bucket','')
            rr['priority'] = m.get('priority','')
            rr['protein_bound_manifest'] = m.get('protein_bound','')
            rr['condition_risk'] = m.get('modification_ligand_risk','')
            rr['baseline_subset'] = baseline_subset(m)
            seq = step_seq_from_row(r)
            rr['step_seq'] = seq
            rr['step_category'] = step_category(seq)
            comps = [r.get('i_label_comp_id_1',''), r.get('j_label_comp_id_1',''), r.get('i_label_comp_id_2',''), r.get('j_label_comp_id_2','')]
            rr['contains_8OG'] = 'yes' if contains_8og_comps(comps) else 'no'
            for k in NUM_FIELDS_STEP:
                v = as_float(r.get(k))
                rr[k] = '' if v is None else f'{v:.6g}'
            tilt = as_float(rr.get('tilt'))
            roll = as_float(rr.get('roll'))
            rr['local_bend'] = '' if tilt is None or roll is None else f'{math.sqrt(tilt*tilt + roll*roll):.6g}'
            rr['geometry_class'] = geometry_class(rr)
            rr['terminal_step'] = 'yes' if is_terminal_step(r) else 'no'
            all_step_rows.append(rr)

        for r in bp:
            rr = dict(r)
            rr['pdb_id'] = pdb_id
            rr['bucket'] = m.get('bucket','')
            rr['priority'] = m.get('priority','')
            rr['protein_bound_manifest'] = m.get('protein_bound','')
            rr['baseline_subset'] = baseline_subset(m)
            comps = [r.get('i_label_comp_id',''), r.get('j_label_comp_id','')]
            rr['contains_8OG'] = 'yes' if contains_8og_comps(comps) else 'no'
            rr['base_pair_class'] = base_letter(r.get('i_label_comp_id')) + ':' + base_letter(r.get('j_label_comp_id'))
            # Opposite base to 8OG, if present
            opp = ''
            if rr['contains_8OG'] == 'yes':
                if contains_8og_comps([r.get('i_label_comp_id','')]): opp = base_letter(r.get('j_label_comp_id'))
                elif contains_8og_comps([r.get('j_label_comp_id','')]): opp = base_letter(r.get('i_label_comp_id'))
            rr['opposite_base_to_8OG'] = opp
            for k in NUM_FIELDS_BP:
                v = as_float(r.get(k))
                rr[k] = '' if v is None else f'{v:.6g}'
            all_bp_rows.append(rr)

    # mark final terminal step per pdb/model/chain group simply by maximum step_number per pdb_id
    max_step = defaultdict(int)
    for r in all_step_rows:
        num = as_float(r.get('step_number'))
        if num is not None:
            max_step[r['pdb_id']] = max(max_step[r['pdb_id']], int(num))
    for r in all_step_rows:
        num = as_float(r.get('step_number'))
        if num is not None and int(num) == max_step[r['pdb_id']]:
            r['terminal_step'] = 'yes'

    # Flanking step inference for lesion structures by step_number +/-1 next to contains_8OG.
    eight_steps_by_pdb = defaultdict(set)
    for r in all_step_rows:
        if r.get('contains_8OG') == 'yes':
            num = as_float(r.get('step_number'))
            if num is not None:
                eight_steps_by_pdb[r['pdb_id']].add(int(num))
    for r in all_step_rows:
        num = as_float(r.get('step_number'))
        flanking = False
        if num is not None:
            flanking = any(abs(int(num) - s) == 1 for s in eight_steps_by_pdb.get(r['pdb_id'], set()))
        r['flanks_8OG_sequence'] = 'yes' if flanking and r.get('contains_8OG') != 'yes' else 'no'
        r['near_8OG'] = 'yes' if r.get('contains_8OG') == 'yes' or r.get('flanks_8OG_sequence') == 'yes' else 'no'

    # Write all tables
    step_fields = [
        'pdb_id','bucket','priority','baseline_subset','protein_bound_manifest','condition_risk',
        'step_number','step_name','step_seq','step_category','contains_8OG','flanks_8OG_sequence','near_8OG','terminal_step','geometry_class',
        'shift','slide','rise','tilt','roll','twist','local_bend',
        'i_label_asym_id_1','i_label_comp_id_1','i_label_seq_id_1','j_label_asym_id_1','j_label_comp_id_1','j_label_seq_id_1',
        'i_label_asym_id_2','i_label_comp_id_2','i_label_seq_id_2','j_label_asym_id_2','j_label_comp_id_2','j_label_seq_id_2'
    ]
    bp_fields = [
        'pdb_id','bucket','priority','baseline_subset','protein_bound_manifest',
        'pair_number','pair_name','base_pair_class','contains_8OG','opposite_base_to_8OG',
        'shear','stretch','stagger','buckle','propeller','opening',
        'i_label_asym_id','i_label_comp_id','i_label_seq_id','j_label_asym_id','j_label_comp_id','j_label_seq_id'
    ]
    write_csv(OUT_DIR / 'base_pair_steps_all.csv', all_step_rows, step_fields)
    write_csv(OUT_DIR / 'base_pair_internal_all.csv', all_bp_rows, bp_fields)

    normal_steps = [r for r in all_step_rows if r['bucket'].startswith('normal')]
    lesion_steps = [r for r in all_step_rows if not r['bucket'].startswith('normal')]
    lesion_bp = [r for r in all_bp_rows if not r['bucket'].startswith('normal')]
    write_csv(OUT_DIR / 'normal_steps.csv', normal_steps, step_fields)
    write_csv(OUT_DIR / 'lesion_steps.csv', lesion_steps, step_fields)
    write_csv(OUT_DIR / 'base_pair_internal.csv', lesion_bp, bp_fields)

    # Baseline summary by subset/category and parameter
    subsets = {
        'core_original': lambda r: r['baseline_subset'] == 'core_original' and r['geometry_class'] == 'plausible',
        'expanded_main': lambda r: r['baseline_subset'] in {'core_original','expanded_main'} and r['geometry_class'] == 'plausible',
        'expanded_main_nonterminal': lambda r: r['baseline_subset'] in {'core_original','expanded_main'} and r['geometry_class'] == 'plausible' and r['terminal_step'] == 'no',
        'all_normal_sensitivity': lambda r: r['bucket'].startswith('normal') and r['geometry_class'] == 'plausible',
    }
    summary_rows: list[dict[str, Any]] = []
    params = ['shift','slide','rise','tilt','roll','twist','local_bend']
    for subset_name, filt in subsets.items():
        rows = [r for r in normal_steps if filt(r)]
        for cat in ['ALL','AT-rich','GC-rich','mixed']:
            cat_rows = rows if cat == 'ALL' else [r for r in rows if r['step_category'] == cat]
            for p in params:
                vals = [as_float(r.get(p)) for r in cat_rows]
                vals = [v for v in vals if v is not None]
                if vals:
                    sd = stdev(vals) if len(vals) > 1 else 0.0
                    summary_rows.append({
                        'baseline_subset': subset_name,
                        'step_category': cat,
                        'parameter': p,
                        'n': len(vals),
                        'mean': f'{mean(vals):.8g}',
                        'sd_sample': f'{sd:.8g}',
                        'min': f'{min(vals):.8g}',
                        'max': f'{max(vals):.8g}',
                    })
                else:
                    summary_rows.append({'baseline_subset': subset_name,'step_category':cat,'parameter':p,'n':0,'mean':'','sd_sample':'','min':'','max':''})
    write_csv(OUT_DIR / 'baseline_summary.csv', summary_rows, ['baseline_subset','step_category','parameter','n','mean','sd_sample','min','max'])

    # Add z-scores and D_step for lesion rows using expanded_main category baseline (fallback ALL).
    baseline_lookup: dict[tuple[str,str], tuple[float,float,int]] = {}
    for r in summary_rows:
        if r['baseline_subset'] == 'expanded_main' and int(r['n']) >= 3 and r['sd_sample'] not in {'','0','0.0'}:
            baseline_lookup[(r['step_category'], r['parameter'])] = (float(r['mean']), float(r['sd_sample']), int(r['n']))
    all_lookup = {(r['parameter']): (float(r['mean']), float(r['sd_sample']), int(r['n'])) for r in summary_rows if r['baseline_subset']=='expanded_main' and r['step_category']=='ALL' and r['sd_sample'] not in {'','0','0.0'}}
    scored_lesion = []
    for r in lesion_steps:
        rr = dict(r)
        zsq_sum = 0.0
        z_fields = []
        used_cat = rr.get('step_category','') if rr.get('step_category') in {'AT-rich','GC-rich','mixed'} else 'ALL'
        for p in ['shift','slide','rise','tilt','roll','twist']:
            x = as_float(rr.get(p))
            baseline = baseline_lookup.get((used_cat,p)) or all_lookup.get(p)
            z = None
            if x is not None and baseline and baseline[1] > 0:
                z = (x - baseline[0]) / baseline[1]
                zsq_sum += z*z
            rr[f'z_{p}'] = '' if z is None else f'{z:.6g}'
            z_fields.append(f'z_{p}')
        rr['D_step'] = f'{math.sqrt(zsq_sum):.6g}' if zsq_sum > 0 else ''
        rr['baseline_used'] = used_cat
        scored_lesion.append(rr)
    score_fields = step_fields + ['baseline_used','z_shift','z_slide','z_rise','z_tilt','z_roll','z_twist','D_step']
    write_csv(OUT_DIR / 'lesion_steps_scored.csv', scored_lesion, score_fields)

    # internal geometry summary and D_pair for 8OG-containing bp
    # Normal internal baselines from all normal plausible base pairs with canonical categories.
    def bp_class_norm(bpclass: str) -> str:
        b = bpclass.replace(':','')
        if set(b) <= {'A','T'} and len(b)==2: return 'AT_pair'
        if set(b) <= {'G','C'} and len(b)==2: return 'GC_pair'
        return 'other'
    normal_bp = [r for r in all_bp_rows if r['bucket'].startswith('normal')]
    bp_baseline: dict[tuple[str,str], tuple[float,float,int]] = {}
    for cls in ['AT_pair','GC_pair']:
        rows = [r for r in normal_bp if bp_class_norm(r['base_pair_class']) == cls]
        for p in ['stretch','opening','propeller']:
            vals = [as_float(r.get(p)) for r in rows]
            vals = [abs(v) for v in vals if v is not None]
            if len(vals) >= 3:
                bp_baseline[(cls,p)] = (mean(vals), stdev(vals), len(vals))
    scored_bp = []
    for r in lesion_bp:
        rr = dict(r)
        cls = 'AT_pair' if rr.get('opposite_base_to_8OG') == 'A' else 'GC_pair' if rr.get('opposite_base_to_8OG') == 'C' else bp_class_norm(rr.get('base_pair_class',''))
        zsq=0.0
        anyz=False
        for p in ['stretch','opening','propeller']:
            x = as_float(rr.get(p))
            baseline = bp_baseline.get((cls,p))
            z=None
            if x is not None and baseline and baseline[1] > 0:
                z = (abs(x) - baseline[0]) / baseline[1]
                zsq += z*z
                anyz=True
            rr[f'z_abs_{p}'] = '' if z is None else f'{z:.6g}'
        rr['D_pair'] = f'{math.sqrt(zsq):.6g}' if anyz else ''
        rr['pair_baseline_used'] = cls
        scored_bp.append(rr)
    bp_score_fields = bp_fields + ['pair_baseline_used','z_abs_stretch','z_abs_opening','z_abs_propeller','D_pair']
    write_csv(OUT_DIR / 'base_pair_internal_scored.csv', scored_bp, bp_score_fields)

    write_csv(OUT_DIR / 'week2_parse_log.csv', parse_log, ['pdb_id','status','base_pair_rows','step_rows','notes'])
    # Console summary
    print('Parsed files:', len(parse_log))
    print('Total step rows:', len(all_step_rows), Counter([r['geometry_class'] for r in all_step_rows]))
    print('Normal step rows:', len(normal_steps), Counter([r['geometry_class'] for r in normal_steps]))
    print('Lesion step rows:', len(lesion_steps), Counter([r['geometry_class'] for r in lesion_steps]))
    for subset in subsets:
        n = len([r for r in normal_steps if subsets[subset](r)])
        print(f'Baseline plausible {subset}: {n}')
    print('Wrote data_processed/*.csv')

if __name__ == '__main__':
    main()
