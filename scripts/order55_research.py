"""Record independently verified incumbents and exact incomplete screens."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sympy as sp
from maxdet.order55 import (exact55, cyclotomic_norms, correlations, fold,
                            affine_canonical, sector_moments, ryser_global,
                            balanced_squares, fourth_bound, PRIMES)

ROOT = Path(__file__).resolve().parents[1]


def record(path):
    raw = json.loads(path.read_text(encoding='utf-8-sig'))
    word = raw.get('bits', raw.get('word'))
    k = word.count('1')
    det, residues = exact55(word)
    matrix = sp.Matrix(55, 55, lambda i, j: int(word[(j-i) % 55]))
    independent = int(matrix.det(method='domain-ge'))
    norms = cyclotomic_norms(word)
    assert det == independent == k*norms[0]*norms[1]*norms[2]
    complement = ''.join('1' if c == '0' else '0' for c in word)
    comp_det, comp_residues = exact55(complement)
    assert comp_det*k == det*(55-k)
    commit = subprocess.check_output(['git', '-c', f'safe.directory={ROOT.as_posix()}',
                                      'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    result = dict(word=word, weight=k, determinant=str(det), complement=complement,
                  complement_determinant=str(comp_det), primes=PRIMES,
                  modular_residues=residues, complement_residues=comp_residues,
                  cyclotomic_norms=list(map(str, norms)),
                  affine_canonical=affine_canonical(word),
                  correlations=correlations(word), mod5=fold(word, 5), mod11=fold(word, 11),
                  sector_moments={m: list(map(str, v)) for m, v in sector_moments(word).items()},
                  search_log=str(path.relative_to(ROOT)), search_log_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  search_method=raw.get('method', 'ALETHEIA generic annealing and local search; see command log'),
                  seed=raw.get('seed'), timestamp=datetime.now(timezone.utc).isoformat(), commit=commit,
                  verification=['finite-field Fourier signed CRT', 'SymPy domain-ge', 'cyclotomic resultants'],
                  status='incumbent only; no global maximality claim')
    with (ROOT/'results/order55_incumbents.jsonl').open('a', encoding='utf8') as f:
        f.write(json.dumps(result)+'\n')
    print(json.dumps(result), flush=True)


def screen():
    records = [json.loads(line) for line in (ROOT/'results/order55_incumbents.jsonl').read_text().splitlines()]
    incumbent = max(max(int(r['determinant']), int(r['complement_determinant'])) for r in records)
    weights = []
    for k in range(1, 28):
        baseline = balanced_squares(k*(k-1)//2, 27)
        ryser = ryser_global(k)
        bound = fourth_bound(k, baseline)
        row = dict(weight=k, complement_weight=55-k, ryser_floor=int(ryser),
                   minimum_half_correlation_squares=baseline, fourth_floor=int(bound),
                   excluded=int(min(ryser, bound)) < incumbent)
        if not row['excluded']:
            threshold = baseline
            while int(fourth_bound(k, threshold+1)) >= incumbent:
                threshold += 1
            row['maximum_half_correlation_squares_not_excluded'] = threshold
            row['excess_square_budget'] = threshold-baseline
            row['first_excluded_floor'] = int(fourth_bound(k, threshold+1))
        weights.append(row)
    result = dict(status='INCOMPLETE: weight screen only', incumbent=str(incumbent), weights=weights)
    (ROOT/'results/order55_weight_screen.json').write_text(json.dumps(result, indent=2), encoding='utf8')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('paths', nargs='*', type=Path)
    parser.add_argument('--screen', action='store_true')
    args = parser.parse_args()
    for path in args.paths:
        record(path.resolve())
    if args.screen:
        screen()
