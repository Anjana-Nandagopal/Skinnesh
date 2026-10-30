#!/usr/bin/env python
"""One-command final pipeline after datasets have been prepared."""
from __future__ import annotations
import subprocess,sys
from pathlib import Path


def run(cmd):
    print('\n$',' '.join(cmd)); subprocess.run(cmd,check=True)

run([sys.executable,'-m','skinnesh.cli','validate','--validation','data/validation.csv','--references','data/reference_actives.csv','--output','results/validation'])
run([sys.executable,'-m','skinnesh.cli','rank','--input','data/plant_compounds.csv','--references','data/reference_actives.csv','--output','results/final'])
run([sys.executable,'-m','skinnesh.cli','report','--ranked','results/final/ranked_compounds.csv','--output','results/final/report.html'])
print('\nFINAL OUTPUTS: results/validation/ and results/final/')
