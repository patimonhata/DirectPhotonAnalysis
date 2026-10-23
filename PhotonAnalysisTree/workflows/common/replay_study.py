#!/usr/bin/env python3
"""Replay a recorded historical study without overwriting its saved results."""
import argparse
import json
import os
from pathlib import Path
import subprocess

module = Path(__file__).resolve().parents[2]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('workflow', choices=('conditional_efficiency', 'conditional_efficiency_reweighted', 'pair_matching_efficiency'))
p.add_argument('configuration')
p.add_argument('--check', action='store_true', help='validate input and recipe without running ROOT')
p.add_argument('--output-base', type=Path, help='required new output base for a replay')
a = p.parse_args()
recipes = json.loads((module / 'workflows' / 'studies' / a.workflow / 'reproduction.json').read_text())
if a.configuration not in recipes:
    p.error('Unknown configuration. Available: ' + ', '.join(sorted(recipes)))
r = recipes[a.configuration]
input_path = module / r['input']
macro = module / r['macro']
if not input_path.is_file() or not macro.is_file():
    p.error('Recorded input or macro is missing: ' + str(input_path))
if a.check:
    print(json.dumps({'configuration': a.configuration, 'input': str(input_path), 'function': r['function'], 'arguments': r['arguments']}, indent=2))
    raise SystemExit(0)
if a.output_base is None:
    p.error('--output-base is required; saved results are never overwritten')
out = a.output_base.resolve()
if out.parent.exists() and any(out.parent.glob(out.name + '*')):
    p.error('Output base already has files: ' + str(out))
out.parent.mkdir(parents=True, exist_ok=True)
args = [json.dumps(str(input_path)), json.dumps(str(out))] + [json.dumps(v) for v in r['arguments']]
expression = 'gSystem->Exit(' + r['function'] + '(' + ','.join(args) + '));'
env = dict(os.environ)
env.pop('BASH_ENV', None)
subprocess.run(['bash', '--noprofile', '--norc', '-c',
    'set -e; source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565; exec "$@"', 'replay-study',
    'root', '-l', '-b', '-q', '-e', '.L ' + str(macro), '-e', expression], check=True, env=env, cwd=module.parent)
