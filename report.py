from __future__ import annotations
from pathlib import Path
import pandas as pd

TARGET_LABELS={
 'human_tyrosinase':'Human tyrosinase',
 'mushroom_tyrosinase':'Mushroom tyrosinase',
 'collagenase':'Collagenase',
 'elastase':'Elastase',
}

def make_report(ranked_csv, output_html, validation_csv=None):
    df=pd.read_csv(ranked_csv)
    out=Path(output_html); out.parent.mkdir(parents=True,exist_ok=True)
    cols=['rank','name','canonical_smiles','MW','logP','TPSA','topical_suitability',
          'human_tyrosinase_activity','mushroom_tyrosinase_activity','collagenase_activity',
          'elastase_activity','activity_mean','overall_score']
    cols=[c for c in cols if c in df.columns]
    table=df[cols].head(25).to_html(index=False,float_format=lambda x:f'{x:.3f}')
    val_html=''
    if validation_csv and Path(validation_csv).exists():
        v=pd.read_csv(validation_csv)
        val_html='<h2>Validation summary</h2>'+v.to_html(index=False,float_format=lambda x:f'{x:.3f}')
    n=len(df)
    html=f'''<!doctype html><html><head><meta charset="utf-8"><title>Skinnesh report</title>
<style>body{{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Arial,sans-serif;margin:36px;line-height:1.45}}
h1{{margin-bottom:4px}} .note{{background:#f4f4f4;padding:14px;border-radius:10px}}
table{{border-collapse:collapse;width:100%;font-size:11px}}th,td{{border:1px solid #ddd;padding:6px;text-align:left}}th{{background:#eee}}
.small{{color:#666;font-size:13px}} .pill{{display:inline-block;background:#eee;padding:5px 9px;border-radius:999px;margin-right:6px}}</style></head>
<body><h1>Skinnesh — ranked plant compounds</h1>
<p>Computational prioritisation for skin-relevant enzyme activity and topical suitability.</p>
<p><span class="pill">{n} input compounds</span><span class="pill">4 target scores</span><span class="pill">Transparent ranking</span></p>
<div class="note"><b>Interpretation:</b> scores are computational prioritisation signals. They are not proof of enzyme inhibition,
topical delivery, safety, efficacy, or clinical/cosmetic benefit. Top-ranked molecules are candidates for experimental follow-up.</div>
<h2>Top candidates</h2>{table}
{val_html}
<h2>Biological context</h2><p>Tyrosinase-related scores are intended to support prioritisation around melanin-associated enzyme activity.
Collagenase and elastase scores are intended to support prioritisation around extracellular-matrix and elastin-related enzyme activity.
The computational ranking should be interpreted together with the experimental literature and target-specific assay context.</p>
<p class="small">Skinnesh keeps human and mushroom tyrosinase as separate target channels. Similarity scores are based on Morgan fingerprints and Tanimoto similarity to the curated reference-active set. Topical suitability is a transparent heuristic based on physicochemical descriptors; it is not a skin-permeation prediction.</p>
</body></html>'''
    out.write_text(html,encoding='utf-8')
