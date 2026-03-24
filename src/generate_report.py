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

ATTR_ORDER = ['Physicochemical', 'Structure', 'Domain', 'Function', 'Modification', 'PPI']

# ── AA full names ──────────────────────────────────────────────────────────────
AA_FULLNAME = {
    'A':'Alanine','R':'Arginine','N':'Asparagine','D':'Aspartic acid',
    'C':'Cysteine','E':'Glutamic acid','Q':'Glutamine','G':'Glycine',
    'H':'Histidine','I':'Isoleucine','L':'Leucine','K':'Lysine',
    'M':'Methionine','F':'Phenylalanine','P':'Proline','S':'Serine',
    'T':'Threonine','W':'Tryptophan','Y':'Tyrosine','V':'Valine',
}

# ── Grantham distance matrix (partial, for demonstration) ──────────────────────
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
    'PI:HB_intra': 'Intra-protein hydrogen bond',
    'PI:SB_intra': 'Intra-protein salt bridge',
    'PI:DS_intra': 'Intra-protein disulfide bond',
    'PI:NB_intra': 'Intra-protein non-bonded interaction',
    'PPI:HB_inter': 'Inter-protein hydrogen bond',
    'PPI:SB_inter': 'Inter-protein salt bridge',
    'PPI:DS_inter': 'Inter-protein disulfide bond',
    'PPI:NB_inter': 'Inter-protein non-bonded interaction',
}

def _feature_label(feat, raw_df, ref_aa, alt_aa):
    if feat in FEATURE_LABELS:
        return FEATURE_LABELS[feat]
    if feat.startswith('RefAA:'):
        return f'Reference amino acid class: {feat.replace("RefAA:", "")}'
    if feat.startswith('AAchange:'):
        return f'Change in amino acid class: {feat.replace("AAchange:", "").replace(">", " → ")}'
    if feat.startswith('Grantham:'):
        sev  = feat.replace('Grantham:', '')
        dist = raw_df.loc[raw_df['feature']=='raw_Grantham','value'].iloc[0]
        return f"Grantham's distance: {sev} (D = {int(dist)})"
    if feat.startswith('RSA:'):
        rsa = raw_df.loc[raw_df['feature']=='raw_RSA','value'].iloc[0]
        return f'Solvent accessibility: {feat.replace("RSA:", "")} (RSA = {rsa:.2f})'
    if feat.startswith('pLDDT:'):
        plddt = raw_df.loc[raw_df['feature']=='raw_pLDDT','value'].iloc[0]
        return f'AlphaFold2 confidence: {feat.replace("pLDDT:", "")} (pLDDT = {plddt:.2f})'
    if feat.startswith('Function:'):    return feat.replace('Function:', '')
    if feat.startswith('Domain:'):      return feat.replace('Domain:', '')
    if feat.startswith('Modification:'): return feat.replace('Modification:', '')
    return feat

FEAT_INTERPRETATION_RANGE = {
    'RSA:Core'              :       "Residue is at the protein core (RSA ≤ 0.05)",
    'RSA:Buried'            :       "Residue is buried (0.05 < RSA ≤ 0.15)",
    'RSA:Medium-buried'     :  "Residue is medium-buried (0.15 < RSA ≤ 0.25)",
    'RSA:Medium-exposed'    :      "Residue is medium-exposed (0.25 < RSA ≤ 0.50)",
    'RSA:Exposed'           :      "Residue is exposed (RSA > 0.50)",
    'pLDDT:Very low'        :   "Very low AlphaFold2 confidence score (pLDDT < 50)",
    'pLDDT:Low'             :        "Low AlphaFold2 confidence score (pLDDT 50–70)",
    'pLDDT:High'            :       "High AlphaFold2 confidence score (pLDDT 70–90)",
    'pLDDT:Very high'       :  "Very high AlphaFold2 confidence score (pLDDT > 90)",
    'Grantham:Mild'         : "Mild physicochemical property shift (Grantham distance 0–50)",
    'Grantham:Moderate'     : "Moderate physicochemical property shift (Grantham distance 51–100)",
    'Grantham:Substantial'  : "Substantial physicochemical property shift (Grantham distance 101–150)",
    'Grantham:Severe'       : "Severe physicochemical property shift (Grantham distance >150)"
}
def _feature_interpretation(feat, raw_df, ref_aa, alt_aa, var_pos, OR):
    if feat in FEAT_INTERPRETATION_RANGE:
        return FEAT_INTERPRETATION_RANGE[feat]
    if feat.startswith('SS:'):     return f'{ref_aa}{var_pos} adopts {FEATURE_LABELS.get(feat, feat)} conformation according to the 9-class DSSP algorithm'
    if feat == 'PI:HB_intra':     return f'{ref_aa}{var_pos} participates in intra-protein hydrogen bonds'
    if feat == 'PI:SB_intra':     return f'{ref_aa}{var_pos} participates in intra-protein salt bridge'
    if feat == 'PI:DS_intra':     return f'{ref_aa}{var_pos} participates in intra-protein disulfide bond'
    if feat == 'PI:NB_intra':     return f'{ref_aa}{var_pos} participates in intra-protein non-bonded contacts'
    if feat == 'PPI:HB_inter':    return f'{ref_aa}{var_pos} participates in inter-protein hydrogen bond'
    if feat == 'PPI:SB_inter':    return f'{ref_aa}{var_pos} participates in inter-protein salt bridge'
    if feat == 'PPI:DS_inter':    return f'{ref_aa}{var_pos} participates in inter-protein disulfide bond'
    if feat == 'PPI:NB_inter':    return f'{ref_aa}{var_pos} participates in inter-protein non-bonded contacts'
    if feat.startswith('RefAA:'):
        if ref_aa == 'G':
            return f'{ref_aa}{var_pos} is a {feat.replace("RefAA:", "").lower()} amino acid (smallest amino acid, introducing high flexibility in structures)'
        elif ref_aa == 'P':
            return f'{ref_aa}{var_pos} is a {feat.replace("RefAA:", "").lower()} amino acid (no backbone hydrogen, often introduces rigid kinks in structures)'
        elif ref_aa == 'C':
            return f'{ref_aa}{var_pos} is a {feat.replace("RefAA:", "").lower()} amino acid (a very reactive sulfhydryl group can form disulfide bonds)'
        else:
            return f'{ref_aa}{var_pos} is a {feat.replace("RefAA:", "").lower()} amino acid'
    if feat.startswith('AAchange:'):
        return f'{feat.split(">")[0].replace("AAchange:", "")} amino acid ({ref_aa}{var_pos}) is substituted to {feat.split(">")[1].lower()} amino acid ({alt_aa}{var_pos})'
    if feat.startswith('Function:'):
        return f'{ref_aa}{var_pos} is annotated as {feat.replace("Function:", "").lower()} (see UniProtKB entry for details)'
    if feat.startswith('Domain:'):
        return f'{ref_aa}{var_pos} is annotated as {feat.replace("Domain:", "").lower()} (see UniProtKB entry for details)'
    if feat.startswith('Modification:'):
        if feat in ['Modification:Lipidation', 'Modification:Glycosylation', 'Modification:Cross-link', 'Modification:Modified residue']:
            return f'{ref_aa}{var_pos} is annotated as {feat.replace("Modification:", "").lower()} site (see UniProtKB entry for details)'
        else:
            return f'{ref_aa}{var_pos} is annotated as {feat.replace("Modification:", "").lower()} site (see PhosphoSitePlus entry for details)'
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

def _build_features(sig_df, raw_df, ref_aa, alt_aa, var_pos):
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
                    'label':          _feature_label(feat, raw_df, ref_aa, alt_aa),
                    'OR':             row['OR'],
                    'q_value':        row['q_value'],
                    'interpretation': _feature_interpretation(feat, raw_df, ref_aa, alt_aa, var_pos, row['OR']),
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
                    interp_df, sig_df, raw_df, fig, out_dir, template_path='.',
                    landscape_png_name=None):
    os.makedirs(out_dir, exist_ok=True)

    scores    = _build_scores(interp_df)
    features  = _build_features(sig_df, raw_df,  ref_aa, alt_aa, var_pos)
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