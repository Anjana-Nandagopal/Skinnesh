#!/usr/bin/env python
"""Convert an IMPPAT chemical-information CSV into Skinnesh input.

The exact column names can vary by IMPPAT release. Pass --name-column and
--smiles-column when necessary. The script does not redistribute IMPPAT data;
it operates on the file you downloaded under IMPPAT's terms.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from skinnesh.core import standardize_smiles


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',required=True)
    ap.add_argument('--output',default='data/plant_compounds.csv')
    ap.add_argument('--name-column',default='name')
    ap.add_argument('--smiles-column',default='smiles')
    args=ap.parse_args()
    df=pd.read_csv(args.input)
    if args.smiles_column not in df.columns:
        raise ValueError(f"SMILES column '{args.smiles_column}' not found. Available: {list(df.columns)}")
    name=df[args.name_column] if args.name_column in df.columns else pd.Series([f'plant_compound_{i+1}' for i in range(len(df))])
    out=pd.DataFrame({'name':name.astype(str),'smiles':df[args.smiles_column].astype(str)})
    out[['canonical_smiles','status']]=out['smiles'].apply(lambda s: pd.Series(standardize_smiles(s)))
    out=out[out.status=='ok'].copy().drop_duplicates('canonical_smiles')
    out[['name','canonical_smiles']].rename(columns={'canonical_smiles':'smiles'}).to_csv(args.output,index=False)
    print(f'Wrote {len(out)} standardized plant compounds to {args.output}')

if __name__=='__main__': main()
