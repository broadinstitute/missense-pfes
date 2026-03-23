"""
pfes_batch_scorer.py  —  Batch PFES scorer
Usage:
    python pfes_batch_scorer.py <input.tsv> <output.tsv>

Input TSV columns : Gene, UniProt, ResID, RefAA, AltAA
Output TSV        : same + PFES, PFES_Physicochemical, PFES_Structure,
                    PFES_Domain, PFES_Function, PFES_Modification, PFES_PPI
"""

import argparse
import sys
import requests
import numpy as np
import pandas as pd
from io import StringIO

import g2papi


# ── Amino acid table ───────────────────────────────────────────────────────────
AA_TABLE = pd.DataFrame([
    {'one': 'A', 'prp': 'Aliphatic'},
    {'one': 'R', 'prp': 'Positively-Charged'},
    {'one': 'N', 'prp': 'Polar/Neutral'},
    {'one': 'D', 'prp': 'Negatively-Charged'},
    {'one': 'C', 'prp': 'Special'},
    {'one': 'E', 'prp': 'Negatively-Charged'},
    {'one': 'Q', 'prp': 'Polar/Neutral'},
    {'one': 'G', 'prp': 'Special'},
    {'one': 'H', 'prp': 'Positively-Charged'},
    {'one': 'I', 'prp': 'Aliphatic'},
    {'one': 'L', 'prp': 'Aliphatic'},
    {'one': 'K', 'prp': 'Positively-Charged'},
    {'one': 'M', 'prp': 'Aliphatic'},
    {'one': 'F', 'prp': 'Aromatic'},
    {'one': 'P', 'prp': 'Special'},
    {'one': 'S', 'prp': 'Polar/Neutral'},
    {'one': 'T', 'prp': 'Polar/Neutral'},
    {'one': 'W', 'prp': 'Aromatic'},
    {'one': 'Y', 'prp': 'Aromatic'},
    {'one': 'V', 'prp': 'Aliphatic'},
])

DISTANCE_MATRIX = {
    'A': {'A':0,'R':112,'N':111,'D':126,'C':195,'Q':91,'E':107,'G':60,'H':86,'I':94,'L':96,'K':106,'M':84,'F':113,'P':27,'S':99,'T':58,'W':148,'Y':112,'V':64},
    'R': {'A':112,'R':0,'N':86,'D':96,'C':180,'Q':43,'E':54,'G':125,'H':29,'I':97,'L':102,'K':26,'M':91,'F':97,'P':103,'S':110,'T':71,'W':101,'Y':77,'V':96},
    'N': {'A':111,'R':86,'N':0,'D':23,'C':139,'Q':46,'E':42,'G':80,'H':68,'I':149,'L':153,'K':94,'M':142,'F':158,'P':91,'S':46,'T':65,'W':174,'Y':143,'V':133},
    'D': {'A':126,'R':96,'N':23,'D':0,'C':154,'Q':61,'E':45,'G':94,'H':81,'I':168,'L':172,'K':101,'M':160,'F':177,'P':108,'S':65,'T':85,'W':181,'Y':160,'V':152},
    'C': {'A':195,'R':180,'N':139,'D':154,'C':0,'Q':154,'E':170,'G':159,'H':174,'I':198,'L':198,'K':202,'M':196,'F':205,'P':169,'S':112,'T':149,'W':215,'Y':194,'V':192},
    'Q': {'A':91,'R':43,'N':46,'D':61,'C':154,'Q':0,'E':29,'G':87,'H':24,'I':109,'L':113,'K':53,'M':101,'F':116,'P':76,'S':68,'T':42,'W':130,'Y':99,'V':96},
    'E': {'A':107,'R':54,'N':42,'D':45,'C':170,'Q':29,'E':0,'G':98,'H':40,'I':134,'L':138,'K':56,'M':126,'F':140,'P':93,'S':80,'T':65,'W':152,'Y':122,'V':121},
    'G': {'A':60,'R':125,'N':80,'D':94,'C':159,'Q':87,'E':98,'G':0,'H':98,'I':135,'L':138,'K':127,'M':127,'F':153,'P':42,'S':56,'T':59,'W':184,'Y':147,'V':109},
    'H': {'A':86,'R':29,'N':68,'D':81,'C':174,'Q':24,'E':40,'G':98,'H':0,'I':94,'L':99,'K':32,'M':87,'F':100,'P':77,'S':89,'T':47,'W':115,'Y':83,'V':84},
    'I': {'A':94,'R':97,'N':149,'D':168,'C':198,'Q':109,'E':134,'G':135,'H':94,'I':0,'L':5,'K':102,'M':10,'F':21,'P':95,'S':142,'T':89,'W':61,'Y':33,'V':29},
    'L': {'A':96,'R':102,'N':153,'D':172,'C':198,'Q':113,'E':138,'G':138,'H':99,'I':5,'L':0,'K':107,'M':15,'F':22,'P':98,'S':145,'T':92,'W':61,'Y':36,'V':32},
    'K': {'A':106,'R':26,'N':94,'D':101,'C':202,'Q':53,'E':56,'G':127,'H':32,'I':102,'L':107,'K':0,'M':95,'F':102,'P':103,'S':121,'T':78,'W':110,'Y':85,'V':97},
    'M': {'A':84,'R':91,'N':142,'D':160,'C':196,'Q':101,'E':126,'G':127,'H':87,'I':10,'L':15,'K':95,'M':0,'F':28,'P':87,'S':135,'T':81,'W':67,'Y':36,'V':21},
    'F': {'A':113,'R':97,'N':158,'D':177,'C':205,'Q':116,'E':140,'G':153,'H':100,'I':21,'L':22,'K':102,'M':28,'F':0,'P':114,'S':155,'T':103,'W':40,'Y':22,'V':50},
    'P': {'A':27,'R':103,'N':91,'D':108,'C':169,'Q':76,'E':93,'G':42,'H':77,'I':95,'L':98,'K':103,'M':87,'F':114,'P':0,'S':74,'T':38,'W':147,'Y':110,'V':68},
    'S': {'A':99,'R':110,'N':46,'D':65,'C':112,'Q':68,'E':80,'G':56,'H':89,'I':142,'L':145,'K':121,'M':135,'F':155,'P':74,'S':0,'T':58,'W':177,'Y':144,'V':124},
    'T': {'A':58,'R':71,'N':65,'D':85,'C':149,'Q':42,'E':65,'G':59,'H':47,'I':89,'L':92,'K':78,'M':81,'F':103,'P':38,'S':58,'T':0,'W':128,'Y':92,'V':69},
    'W': {'A':148,'R':101,'N':174,'D':181,'C':215,'Q':130,'E':152,'G':184,'H':115,'I':61,'L':61,'K':110,'M':67,'F':40,'P':147,'S':177,'T':128,'W':0,'Y':37,'V':88},
    'Y': {'A':112,'R':77,'N':143,'D':160,'C':194,'Q':99,'E':122,'G':147,'H':83,'I':33,'L':36,'K':85,'M':36,'F':22,'P':110,'S':144,'T':92,'W':37,'Y':0,'V':55},
    'V': {'A':64,'R':96,'N':133,'D':152,'C':192,'Q':96,'E':121,'G':109,'H':84,'I':29,'L':32,'K':97,'M':21,'F':50,'P':68,'S':124,'T':69,'W':88,'Y':55,'V':0},
}

MAX_ASA = {
    'A':129.0,'R':274.0,'N':195.0,'D':193.0,'C':167.0,'E':223.0,'Q':225.0,
    'G':104.0,'H':224.0,'I':197.0,'L':201.0,'K':236.0,'M':224.0,'F':240.0,
    'P':159.0,'S':155.0,'T':172.0,'W':285.0,'Y':263.0,'V':174.0,
}

SS9_LABELS    = ['B','E','G','H','I','P','C','S','T']
RSA_BINS      = [0, 0.05, 0.25, 0.50, 0.75, 1.01]
RSA_LABELS    = ['Core','Buried','Medium-buried','Medium-exposed','Exposed']
PLDDT_BINS    = [0, 50, 70, 90, 100.01]
PLDDT_LABELS  = ['Very low','Low','High','Very high']
PCHEM_CLASSES = ['Aliphatic','Aromatic','Polar/Neutral','Positively-Charged','Negatively-Charged','Special']

FUNCTION_FEATURES = ['Active site','Binding site','Site','DNA binding','Zinc finger']
DOMAIN_FEATURES   = [
    'Region/Disordered','Region/Interaction','Region/Others',
    'Motif','Coiled coil','Compositional bias','Repeat','Domain',
    'Topological domain','Transmembrane','Intramembrane',
    'Signal','Transit peptide','Propeptide','Peptide','Chain',
]
PTM_FEATURES = [
    'Acetylation','Methylation','Phosphorylation','SUMOylation',
    'Ubiquitination','O-GalNAc/GlcNAc','Lipidation','Glycosylation',
    'Cross-link','Modified residue',
]

FEATURE_COLS = (
    [f'RefAA:{c}'        for c in PCHEM_CLASSES] +
    [f'SS:{s}'           for s in SS9_LABELS] +
    [f'RSA:{b}'          for b in RSA_LABELS] +
    [f'pLDDT:{b}'        for b in PLDDT_LABELS] +
    ['PI:HB_intra','PI:SB_intra','PI:DS_intra','PI:NB_intra'] +
    [f'Function:{f}'     for f in FUNCTION_FEATURES] +
    [f'Domain:{d}'       for d in DOMAIN_FEATURES] +
    [f'Modification:{m}' for m in PTM_FEATURES] +
    ['PPI:HB_inter','PPI:SB_inter','PPI:DS_inter','PPI:NB_inter']
)

ATTR_PREFIXES = {
    'Physicochemical': ['RefAA:','AAchange:','Grantham:'],
    'Structure':       ['SS:','RSA:','pLDDT:','PI:'],
    'Domain':          ['Domain:'],
    'Function':        ['Function:'],
    'Modification':    ['Modification:'],
    'PPI':             ['PPI:'],
}

Q_THRESHOLD = 0.01


# ── Protein feature annotation ─────────────────────────────────────────────────

def _parse_ss9(val):
    if val == '-' or pd.isna(val): return None
    s = str(val).strip()
    return s[0] if s[0] in SS9_LABELS else None

def _is_annotated(val):
    return str(val).strip() != '-'

def _combine_glyco(row):
    a, b = str(row.get('O-GalNAc','-')).strip(), str(row.get('O-GlcNAc','-')).strip()
    if a == '-' and b == '-': return '-'
    if a == '-': return b
    if b == '-': return a
    return f'{a},{b}'

def _route_disulfide(val, ds_intra, ds_inter):
    if str(val).strip() == '-': return ds_intra, ds_inter
    parts = str(val).split(':', 1)
    if len(parts) == 2 and 'interchain' in parts[1].lower():
        ds_inter = val if ds_inter == '-' else f'{ds_inter},{val}'
    else:
        ds_intra = val if ds_intra == '-' else f'{ds_intra},{val}'
    return ds_intra, ds_inter

def _classify_region(val):
    matched = set()
    if str(val).strip() == '-': return matched
    for entry in str(val).split(';'):
        entry = entry.strip()
        if not entry: continue
        if entry == 'Disordered':           matched.add('Region/Disordered')
        elif 'interact' in entry.lower():   matched.add('Region/Interaction')
        else:                               matched.add('Region/Others')
    return matched

def annotate_protein_features(pf):
    df = pd.DataFrame(False, index=pf.index, columns=['ResID','RefAA'] + FEATURE_COLS)
    df['ResID'] = pf['residueId']
    df['RefAA']  = pf['AA']

    ss_col = 'Secondary structure (DSSP 9-state)*'
    if ss_col in pf.columns:
        parsed = pf[ss_col].apply(_parse_ss9)
        for s in SS9_LABELS:
            df[f'SS:{s}'] = parsed == s

    asa_col = 'Accessible surface area (Å²)*'
    if asa_col in pf.columns:
        asa  = pd.to_numeric(pf[asa_col], errors='coerce')
        rsa  = (asa / pf['AA'].map(MAX_ASA)).clip(0, 1)
        bins = pd.cut(rsa, bins=RSA_BINS, labels=RSA_LABELS, right=False)
        for b in RSA_LABELS:
            df[f'RSA:{b}'] = bins == b

    plddt_col = 'AlphaFold confidence (pLDDT)'
    if plddt_col in pf.columns:
        plddt = pd.to_numeric(pf[plddt_col], errors='coerce')
        bins  = pd.cut(plddt, bins=PLDDT_BINS, labels=PLDDT_LABELS, right=False)
        for b in PLDDT_LABELS:
            df[f'pLDDT:{b}'] = bins == b

    aa2prp  = AA_TABLE.set_index('one')['prp'].to_dict()
    ref_cls = pf['AA'].map(aa2prp)
    for c in PCHEM_CLASSES:
        df[f'RefAA:{c}'] = ref_cls == c

    intra_map = {
        'PI:HB_intra': ['Intra-chain Hydrogen bond (PDB)',          'Intra-chain Hydrogen bond (AlphaFold2)'],
        'PI:SB_intra': ['Intra-chain Salt bridge (PDB)',            'Intra-chain Salt bridge (AlphaFold2)'],
        'PI:DS_intra': ['Intra-chain Disulfide bond (PDB)',         'Intra-chain Disulfide bond (AlphaFold2)'],
        'PI:NB_intra': ['Intra-chain Non-bonded interaction (PDB)', 'Intra-chain Non-bonded interaction (AlphaFold2)'],
    }
    for feat, cols in intra_map.items():
        avail = [c for c in cols if c in pf.columns]
        if avail:
            df[feat] = pf[avail].apply(lambda row: any(_is_annotated(v) for v in row), axis=1)

    if 'Disulfide bond' in pf.columns:
        for i, val in enumerate(pf['Disulfide bond']):
            cur_intra = '-' if not df.at[i, 'PI:DS_intra']  else 'annotated'
            cur_inter = '-' if not df.at[i, 'PPI:DS_inter'] else 'annotated'
            new_intra, new_inter = _route_disulfide(val, cur_intra, cur_inter)
            df.at[i, 'PI:DS_intra']  = new_intra != '-'
            df.at[i, 'PPI:DS_inter'] = new_inter != '-'

    inter_map = {
        'PPI:HB_inter': 'Inter-chain Hydrogen bond (PDB)',
        'PPI:SB_inter': 'Inter-chain Salt bridge (PDB)',
        'PPI:DS_inter': 'Inter-chain Disulfide bond (PDB)',
        'PPI:NB_inter': 'Inter-chain Non-bonded interaction (PDB)',
    }
    for feat, col in inter_map.items():
        if col in pf.columns:
            df[feat] = df[feat] | pf[col].apply(_is_annotated)

    for f in FUNCTION_FEATURES:
        df[f'Function:{f}'] = pf[f].apply(_is_annotated)

    region_cats = pf['Region'].apply(_classify_region)
    for sub in ['Region/Disordered','Region/Interaction','Region/Others']:
        df[f'Domain:{sub}'] = region_cats.apply(lambda s: sub in s)
    for d in [d for d in DOMAIN_FEATURES if not d.startswith('Region/')]:
        df[f'Domain:{d}'] = pf[d].apply(_is_annotated)

    glyco = pf[['O-GalNAc','O-GlcNAc']].apply(_combine_glyco, axis=1)
    df['Modification:O-GalNAc/GlcNAc'] = glyco.apply(_is_annotated)
    for m in [m for m in PTM_FEATURES if m != 'O-GalNAc/GlcNAc']:
        df[f'Modification:{m}'] = pf[m].apply(_is_annotated)

    return df


# ── Variant-level features ─────────────────────────────────────────────────────

def add_variant_features(encoded_refaa, alt_aa):
    df = encoded_refaa.copy()
    aa2prp  = AA_TABLE.set_index('one')['prp'].to_dict()
    alt_cls = aa2prp.get(alt_aa)
    ref_cls = df['RefAA'].map(aa2prp)

    for ref_c in PCHEM_CLASSES:
        for alt_c in PCHEM_CLASSES:
            df[f'AAchange:{ref_c}>{alt_c}'] = (ref_cls == ref_c) & (alt_cls == alt_c)

    dist_bins   = [0, 50, 100, 150, np.inf]
    dist_labels = ['Mild','Moderate','Substantial','Severe']
    dists = df['RefAA'].apply(lambda r: DISTANCE_MATRIX.get(r, {}).get(alt_aa, np.nan))
    grantham_bin = pd.cut(dists, bins=dist_bins, labels=dist_labels, right=False)
    for lbl in dist_labels:
        df[f'Grantham:{lbl}'] = grantham_bin == lbl

    return df


# ── OR lookup ──────────────────────────────────────────────────────────────────

def get_or_lookup(odd_ratio_data):
    sig = odd_ratio_data[odd_ratio_data['q_value'] < Q_THRESHOLD].dropna(subset=['OR'])
    log_ors = {feat: np.log(row['OR']) for feat, row in sig.iterrows()}
    return log_ors


# ── PFES computation ───────────────────────────────────────────────────────────

def compute_pfes_row(feat_row, log_ors):
    scores = {attr: 0.0 for attr in ATTR_PREFIXES}
    for feat, log_or in log_ors.items():
        if feat not in feat_row.index or not feat_row[feat]:
            continue
        for attr, prefixes in ATTR_PREFIXES.items():
            if any(feat.startswith(p) for p in prefixes):
                scores[attr] += log_or
                break
    scores['PFES'] = sum(scores.values())
    return scores


# ── Remote data loaders ────────────────────────────────────────────────────────

_GITHUB_OR_URL = (
    'https://raw.githubusercontent.com/broadinstitute/missense-pfes'
    '/refs/heads/main/results/enrichment_OR_by_protein_class.csv'
)
_GCS_META_URL = (
    'https://storage.googleapis.com/g2p-portal'
    '/portal_data/2026_q1_data/uniprot_metadata.tsv'
)

def _load_resources():
    print("Loading enrichment table...")
    r = requests.get(_GITHUB_OR_URL); r.raise_for_status()
    enrichment_df = pd.read_csv(StringIO(r.text), header=[0, 1], index_col=0)

    print("Loading PANTHER protein class map...")
    r = requests.get(_GCS_META_URL); r.raise_for_status()
    meta = pd.read_csv(StringIO(r.text), sep='\t')
    panther_map = meta.set_index('UniprotKB_Entry')['PANTHER_protein_class'].to_dict()

    return enrichment_df, panther_map


def _resolve_protein_class(uniprot, panther_map, enrichment_df):
    """
    Look up PANTHER class for a UniProt ID.
    Falls back to 'unclassified' if: not found, marked 'not-available',
    or the class is absent from the enrichment table.
    """
    protein_class = panther_map.get(uniprot, 'unclassified')
    if protein_class == 'not-available':
        protein_class = 'unclassified'
    valid_classes = enrichment_df.columns.get_level_values(0)
    if protein_class not in valid_classes:
        print(f"  Warning: protein class '{protein_class}' not in enrichment table, falling back to 'unclassified'")
        protein_class = 'unclassified'
    return protein_class


# ── Per-UniProt cache ──────────────────────────────────────────────────────────

def _get_protein_data(gene, uniprot, enrichment_df, panther_map, cache):
    if uniprot in cache:
        return cache[uniprot]

    print(f"  Fetching protein features: {gene} ({uniprot})")
    protein_class = _resolve_protein_class(uniprot, panther_map, enrichment_df)

    pf = g2papi.get_protein_features(gene, uniprot)
    pf.fillna('-', inplace=True)
    pf.rename(columns=lambda c: c.replace(' (UniProt)', ''), inplace=True)

    encoded = annotate_protein_features(pf)

    odd_ratio_data = enrichment_df[protein_class][['OR','q_value']].copy()
    max_val = odd_ratio_data.loc[np.isfinite(odd_ratio_data['OR']), 'OR'].max()
    min_val = odd_ratio_data.loc[odd_ratio_data['OR'] > 0, 'OR'].min()
    odd_ratio_data['OR'] = np.where(odd_ratio_data['OR'] == np.inf, max_val, odd_ratio_data['OR'])
    odd_ratio_data['OR'] = np.where(odd_ratio_data['OR'] == 0,      min_val, odd_ratio_data['OR'])

    log_ors = get_or_lookup(odd_ratio_data)

    cache[uniprot] = {
        'encoded':        encoded,
        'log_ors':        log_ors,
        'protein_class':  protein_class,
    }
    return cache[uniprot]


# ── Main batch function ────────────────────────────────────────────────────────

SCORE_COLS = ['PFES','PFES_Physicochemical','PFES_Structure',
              'PFES_Domain','PFES_Function','PFES_Modification','PFES_PPI']

def run_pfes_batch(input_tsv, output_tsv):
    df = pd.read_csv(input_tsv, sep='\t')
    required = {'Gene','UniProt','ResID','RefAA','AltAA'}
    missing  = required - set(df.columns)
    if missing:
        sys.exit(f"Error: input TSV missing columns: {missing}")

    enrichment_df, panther_map = _load_resources()
    cache   = {}
    results = []

    for i, row in df.iterrows():
        gene, uniprot = row['Gene'], row['UniProt']
        res_id, ref_aa, alt_aa = int(row['ResID']), row['RefAA'], row['AltAA']
        label = f"{gene}:{ref_aa}{res_id}{alt_aa}"
        try:
            data      = _get_protein_data(gene, uniprot, enrichment_df, panther_map, cache)
            encoded   = data['encoded']
            log_ors   = data['log_ors']

            if not (encoded['ResID'] == res_id).any():
                raise ValueError(f"ResID {res_id} not found")

            enc       = add_variant_features(encoded, alt_aa)
            feat_cols = [c for c in enc.columns if c not in ('ResID','RefAA')]
            feat_row  = enc.loc[enc['ResID'] == res_id].iloc[0][feat_cols]
            scores    = compute_pfes_row(feat_row, log_ors)

            results.append({
                'PFES':                scores['PFES'],
                'PFES_Physicochemical': scores['Physicochemical'],
                'PFES_Structure':       scores['Structure'],
                'PFES_Domain':          scores['Domain'],
                'PFES_Function':        scores['Function'],
                'PFES_Modification':    scores['Modification'],
                'PFES_PPI':             scores['PPI'],
            })
            print(f"  [{i+1}/{len(df)}] {label}  PFES={scores['PFES']:.3f}")

        except Exception as e:
            print(f"  [{i+1}/{len(df)}] {label}  WARNING: {e}")
            results.append({col: np.nan for col in SCORE_COLS})

    out = pd.concat([df, pd.DataFrame(results)], axis=1)
    out.to_csv(output_tsv, sep='\t', index=False)
    print(f"\nDone. {len(out)} variants saved to: {output_tsv}")


# ── CLI ────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute PFES scores for a batch of missense variants.'
    )
    parser.add_argument('input',  help='Input TSV (Gene, UniProt, ResID, RefAA, AltAA)')
    parser.add_argument('output', help='Output TSV path')
    args = parser.parse_args()

    run_pfes_batch(args.input, args.output)