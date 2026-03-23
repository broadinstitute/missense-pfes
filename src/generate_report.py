"""
generate_report.py — Deterministic PFES variant report generator
Produces .md, .html, and a static .png of the landscape figure.

Depends on variables already in scope from compute_PFES.ipynb:
    gene, variant, protein_class, var_pos, ref_aa, alt_aa
    landscape, xi, yi, ATTR_ORDER
    sig_df       (significant features DataFrame)
    interp_df    (partitioning + p-values DataFrame)
    fig          (Plotly figure from run_pfes_landscape)
"""

import os
import numpy as np
import pandas as pd
import markdown
from jinja2 import Environment, FileSystemLoader
from pfes_scorer import DISTANCE_MATRIX

ATTR_ORDER = ['Physicochemical', 'Structure', 'Domain', 'Function', 'Modification', 'PPI']

# ── AA full names ──────────────────────────────────────────────────────────────
AA_FULLNAME = {
    'A':'Alanine','R':'Arginine','N':'Asparagine','D':'Aspartic acid',
    'C':'Cysteine','E':'Glutamic acid','Q':'Glutamine','G':'Glycine',
    'H':'Histidine','I':'Isoleucine','L':'Leucine','K':'Lysine',
    'M':'Methionine','F':'Phenylalanine','P':'Proline','S':'Serine',
    'T':'Threonine','W':'Tryptophan','Y':'Tyrosine','V':'Valine',
}

# ── Feature label and interpretation helpers ───────────────────────────────────
FEATURE_LABELS = {
    'SS:B': 'Secondary structure (B, beta-bridge)',
    'SS:E': 'Secondary structure (E, beta-strand)',
    'SS:G': 'Secondary structure (G, 3-10 helix)',
    'SS:H': 'Secondary structure (H, alpha-helix)',
    'SS:I': 'Secondary structure (I, pi-helix)',
    'SS:P': 'Secondary structure (P, polyproline helix)',
    'SS:C': 'Secondary structure (C, loop/coil)',
    'SS:S': 'Secondary structure (S, bend)',
    'SS:T': 'Secondary structure (T, turn)',
    'RSA:Core':           'Solvent accessibility (Core, RSA < 5%)',
    'RSA:Buried':         'Solvent accessibility (Buried, RSA 5–25%)',
    'RSA:Medium-buried':  'Solvent accessibility (Medium-buried, RSA 25–50%)',
    'RSA:Medium-exposed': 'Solvent accessibility (Medium-exposed, RSA 50–75%)',
    'RSA:Exposed':        'Solvent accessibility (Exposed, RSA > 75%)',
    'pLDDT:Very low':  'AlphaFold2 confidence (Very low, pLDDT < 50)',
    'pLDDT:Low':       'AlphaFold2 confidence (Low, pLDDT 50–70)',
    'pLDDT:High':      'AlphaFold2 confidence (High, pLDDT 70–90)',
    'pLDDT:Very high': 'AlphaFold2 confidence (Very high, pLDDT > 90)',
    'PI:HB_intra': 'Intra-protein hydrogen bond',
    'PI:SB_intra': 'Intra-protein salt bridge',
    'PI:DS_intra': 'Intra-protein disulfide bond',
    'PI:NB_intra': 'Intra-protein non-bonded interaction',
    'PPI:HB_inter': 'Inter-protein hydrogen bond',
    'PPI:SB_inter': 'Inter-protein salt bridge',
    'PPI:DS_inter': 'Inter-protein disulfide bond',
    'PPI:NB_inter': 'Inter-protein non-bonded interaction',
}

def _feature_label(feat, ref_aa, alt_aa):
    if feat in FEATURE_LABELS:
        return FEATURE_LABELS[feat]
    if feat.startswith('RefAA:'):
        return f'Reference residue class ({feat.replace("RefAA:", "")})'
    if feat.startswith('AAchange:'):
        return f'Residue class change ({feat.replace("AAchange:", "").replace(">", " → ")})'
    if feat.startswith('Grantham:'):
        sev  = feat.replace('Grantham:', '')
        dist = DISTANCE_MATRIX.get(ref_aa, {}).get(alt_aa, '?')
        return f'Grantham distance = {dist} ({sev})'
    if feat.startswith('Function:'):    return feat.replace('Function:', '')
    if feat.startswith('Domain:'):      return feat.replace('Domain:', '')
    if feat.startswith('Modification:'): return feat.replace('Modification:', '')
    return feat

def _feature_interpretation(feat, ref_aa, alt_aa, var_pos, OR):
    if feat.startswith('SS:'):     return f'{ref_aa}{var_pos} adopts {FEATURE_LABELS.get(feat, feat)} conformation'
    if feat.startswith('RSA:'):    return f'Residue is {feat.replace("RSA:", "").lower()} (solvent accessibility)'
    if feat.startswith('pLDDT:'): return f'{feat.replace("pLDDT:", "")} AlphaFold2 prediction confidence'
    if feat == 'PI:HB_intra':     return 'Residue participates in intra-protein hydrogen bonds'
    if feat == 'PI:SB_intra':     return 'Residue participates in intra-protein salt bridge'
    if feat == 'PI:DS_intra':     return 'Residue participates in intra-protein disulfide bond'
    if feat == 'PI:NB_intra':     return 'Residue participates in intra-protein non-bonded contacts'
    if feat == 'PPI:HB_inter':    return 'Residue participates in inter-protein hydrogen bond'
    if feat == 'PPI:SB_inter':    return 'Residue participates in inter-protein salt bridge'
    if feat == 'PPI:DS_inter':    return 'Residue participates in inter-protein disulfide bond'
    if feat == 'PPI:NB_inter':    return 'Residue participates in inter-protein non-bonded contacts'
    if feat.startswith('RefAA:'):
        return f'{ref_aa}{var_pos} is a {feat.replace("RefAA:", "").lower()} amino acid'
    if feat.startswith('AAchange:'):
        return f'Substitution class: {feat.replace("AAchange:", "").replace(">", " → ")}'
    if feat.startswith('Grantham:'):
        return f'{feat.replace("Grantham:", "").capitalize()} physicochemical shift'
    if feat.startswith('Function:'):
        return f'{ref_aa}{var_pos} is annotated as {feat.replace("Function:", "").lower()}'
    if feat.startswith('Domain:'):
        return f'{ref_aa}{var_pos} falls within a {feat.replace("Domain:", "").lower()} region'
    if feat.startswith('Modification:'):
        return f'{ref_aa}{var_pos} carries a {feat.replace("Modification:", "").lower()} annotation'
    return 'enriched among pathogenic variants' if OR > 1 else 'enriched among control variants'


# ── Build scores dict from interp_df ──────────────────────────────────────────

def _build_scores(interp_df):
    scores = {}
    for _, row in interp_df.iterrows():
        key = row['attribute'] if row['attribute'] != 'Overall' else 'PFES'
        scores[key] = {
            'score':      row['score'],
            'partition':  row['partition'],
            'p_enriched': None if pd.isna(row['p_enriched']) else row['p_enriched'],
            'p_depleted': None if pd.isna(row['p_depleted']) else row['p_depleted'],
        }
    return scores


# ── Build features dict grouped by attribute ──────────────────────────────────

def _build_features(sig_df, ref_aa, alt_aa, var_pos):
    features = {attr: [] for attr in ATTR_ORDER}
    attr_map = {
        'Physicochemical': ['RefAA:', 'AAchange:', 'Grantham:'],
        'Structure':       ['SS:', 'RSA:', 'pLDDT:', 'PI:'],
        'Domain':          ['Domain:'],
        'Function':        ['Function:'],
        'Modification':    ['Modification:'],
        'PPI':             ['PPI:'],
    }
    for _, row in sig_df.iterrows():
        feat = row['feature']
        for attr, prefixes in attr_map.items():
            if any(feat.startswith(p) for p in prefixes):
                features[attr].append({
                    'label':          _feature_label(feat, ref_aa, alt_aa),
                    'OR':             row['OR'],
                    'q_value':        row['q_value'],
                    'interpretation': _feature_interpretation(feat, ref_aa, alt_aa, var_pos, row['OR']),
                })
                break
    return features


# ── Build narrative paragraph ──────────────────────────────────────────────────

def _build_narrative(scores, gene, variant, ref_aa, alt_aa, var_pos):
    overall   = scores['PFES']
    partition = overall['partition']
    p_val     = overall['p_enriched'] if partition == 'PF-Enriched' else overall.get('p_depleted')
    p_str     = f"*p* = {p_val:.2e}" if p_val is not None else ''

    attr_scores  = {a: scores[a]['score'] for a in ATTR_ORDER if a in scores}
    enriched     = [a for a in ATTR_ORDER if scores.get(a, {}).get('partition') == 'PF-Enriched']
    depleted     = [a for a in ATTR_ORDER if scores.get(a, {}).get('partition') == 'PF-Depleted']

    if partition == 'PF-Enriched':
        dominant_pool, opposing = enriched, depleted
    elif partition == 'PF-Depleted':
        dominant_pool, opposing = depleted, enriched
    else:
        dominant_pool, opposing = [], []

    dominant = (max(dominant_pool, key=lambda a: abs(attr_scores.get(a, 0)))
                if dominant_pool else None)

    if partition == 'PF-Neutral':
        return (f"The variant {variant} in {gene} is classified as **PF-Neutral**, "
                f"indicating that its protein feature profile is statistically consistent "
                f"with both pathogenic and control variant distributions. "
                f"No single attribute shows significant enrichment in either direction.")

    direction = ('enriched among known pathogenic variants' if partition == 'PF-Enriched'
                 else 'depleted of features associated with pathogenic variants')

    s1 = (f"The variant {variant} in {gene} is classified as **{partition}** ({p_str}), "
          f"indicating that its protein feature profile is {direction}.")
    s2 = f"This classification is primarily driven by **{dominant.lower()}** features." if dominant else ''
    s3 = ''
    if opposing:
        opp_str = ' and '.join(f'**{o.lower()}**' for o in opposing)
        s3 = (f"A modest opposing contribution from {opp_str} features "
              f"is present but does not change the overall classification.")

    return ' '.join(s for s in [s1, s2, s3] if s)


# ── Main render function ───────────────────────────────────────────────────────

def generate_report(gene, variant, protein_class, var_pos, ref_aa, alt_aa,
                    interp_df, sig_df, fig, out_dir, template_path='.',
                    landscape_png_name=None):
    os.makedirs(out_dir, exist_ok=True)

    scores    = _build_scores(interp_df)
    features  = _build_features(sig_df, ref_aa, alt_aa, var_pos)
    narrative = _build_narrative(scores, gene, variant, ref_aa, alt_aa, var_pos)

    # Use shared PNG name if provided, otherwise save variant-specific PNG
    png_name = landscape_png_name or f'{gene}_{variant}_landscape.png'
    png_path = os.path.join(out_dir, png_name)
    if not os.path.exists(png_path):
        fig.write_image(png_path, width=1200, height=1000, scale=2)

    context = {
        'gene':          gene,
        'variant':       variant,
        'protein_class': protein_class,
        'var_pos':       var_pos,
        'ref_aa':        ref_aa,
        'alt_aa':        alt_aa,
        'ref_aa_full':   AA_FULLNAME.get(ref_aa, ref_aa),
        'alt_aa_full':   AA_FULLNAME.get(alt_aa, alt_aa),
        'scores':        scores,
        'features':      features,
        'narrative':     narrative,
        'attr_order':    ATTR_ORDER,
        'landscape_html': f'{gene}_{variant}_landscape.html',
        'landscape_png':  png_name,
    }

    env = Environment(loader=FileSystemLoader(template_path))
    env.filters['format_or'] = lambda v: f'{v:.2f}' if v < 1 else f'{v:.1f}'
    template = env.get_template('report_template.md.j2')
    md_text  = template.render(**context)

    # Save markdown
    md_path = os.path.join(out_dir, f'{gene}_{variant}_report.md')
    with open(md_path, 'w') as f:
        f.write(md_text)

    # Convert to HTML
    CSS = """
        body { font-family: Arial, sans-serif; max-width: 900px; margin: 40px auto; font-size: 13px; }
        table { border-collapse: collapse; width: 100%; margin: 16px 0; }
        th, td { border: 1px solid #ccc; padding: 6px 10px; text-align: left; }
        th { background: #f2f2f2; }
        h1 { font-size: 18px; } h2 { font-size: 15px; }
        img { max-width: 100%; }
    """
    html_body = markdown.markdown(md_text, extensions=['tables'])
    html_full = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><style>{CSS}</style></head>
<body>{html_body}</body></html>"""

    html_path = os.path.join(out_dir, f'{gene}_{variant}_report.html')
    with open(html_path, 'w') as f:
        f.write(html_full)

    print(f"Report saved:\n  {md_path}\n  {html_path}\n  {png_path}")
    return md_text


if __name__ == '__main__':
    pass