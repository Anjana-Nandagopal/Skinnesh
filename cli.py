from __future__ import annotations
import argparse
from pathlib import Path
import yaml
from .core import rank_compounds
from .validation import score_validation_set
from .report import make_report


def load_cfg(path):
    with open(path,'r',encoding='utf-8') as f: return yaml.safe_load(f)

def main():
    parser=argparse.ArgumentParser(prog='skinnesh',description='Screen and rank plant compounds for skin-relevant activity and topical suitability.')
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('rank',help='Standardise, score and rank plant compounds')
    p.add_argument('--input',required=True); p.add_argument('--references',required=True); p.add_argument('--output',required=True); p.add_argument('--config',default='skinnesh/config.yaml')
    p=sub.add_parser('validate',help='Evaluate similarity scoring against held-out actives/decoys')
    p.add_argument('--validation',required=True); p.add_argument('--references',required=True); p.add_argument('--output',required=True)
    p=sub.add_parser('report',help='Create HTML report from ranked CSV')
    p.add_argument('--ranked',required=True); p.add_argument('--output',required=True); p.add_argument('--validation-metrics',default=None)
    args=parser.parse_args()
    if args.command=='rank':
        result=rank_compounds(args.input,args.references,args.output,load_cfg(args.config))
        print(result[['rank','name','overall_score']].head(10).to_string(index=False)); print(f'\nSaved: {Path(args.output)/"ranked_compounds.csv"}')
    elif args.command=='validate':
        metrics=score_validation_set(args.validation,args.references,args.output)
        print(metrics.to_string(index=False) if not metrics.empty else 'No target had both classes.')
        print(f'\nSaved validation outputs to: {args.output}')
    else:
        make_report(args.ranked,args.output,args.validation_metrics); print(f'Saved: {args.output}')

if __name__=='__main__': main()
