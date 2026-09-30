"""Regenerate the intermediate pickles / arrays written by the scripts in 02, 03 and 04.

Usage (from the repository root):
    python tools/generate_intermediates.py [--timeout SECONDS] [--only 02|03|04]

Each chain runs its steps in order (later steps read earlier outputs); the chains, one per
directory, run in parallel. Logs go to intermediates/logs/, and intermediates/MANIFEST.txt
records the command, wall time, exit status, output files with SHA-256, and the Python/library
versions.
"""
import argparse, hashlib, os, platform, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
D2, D3, D4 = '02_axis_direct_method', '03_scalefree_direct_method', '04_transition_2d_direct_method'

CHAINS = {
    '02': [(D2, s) for s in ['n4_pert.py', 'n4_elim2.py']],
    '03': [(D3, s) for s in [
        'sf_classify.py', 'sf_classify2.py', 'sf_classify3.py', 'sf_struct.py', 'sf_master2.py',
        'sf_w2c.py', 'sf_mu.py 1',
        'sf_general.py 6 4 1', 'deg6_d4_solve.py',
        'sf_first.py 6 4 1', 'sf_first.py 6 2 1', 'sf_first.py 6 0 1',
        'sf_scan.py 4 1', 'deg6_close.py',
        # flagged in 03/README as long; 'sf_fast.py 6 0 1' exceeded 3 h in the v0.1.1 run
        'sf_fast.py 6 2 1', 'sf_fast.py 6 0 1',
    ]],
    '04': [(D4, s) for s in [
        'fo2d.py 0 2', 'fo2d.py 2 2', 'fo2d.py 4 4', 'fo2d.py 6 6', 'fo2d_m6.py',
        'fo2d2.py hernquist 6 2 2', 'fo2d2.py hernquist 6 6 6',
        'so2d_seq.py 2', 'so2d_seq.py 0', 'so2d_step1.py',
        'so2d_step2.py 2', 'so2d_step2.py 0', 'so2d_step2.py 4',
    ]],
}

OUT, LOGS = ROOT / 'intermediates', ROOT / 'intermediates' / 'logs'


def outputs(d):
    return {p.name: p.stat().st_mtime for p in (ROOT / d).iterdir() if p.suffix in ('.pkl', '.npy')}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run_chain(name, steps, timeout):
    records = []
    for d, cmd in steps:
        before = outputs(d)
        log = LOGS / f"{d[:2]}_{cmd.replace(' ', '_').replace('.py', '')}.log"
        t0 = time.time()
        with open(log, 'w') as fh:
            try:
                rc = subprocess.run([sys.executable, '-u', *cmd.split()], cwd=ROOT / d, stdout=fh,
                                    stderr=subprocess.STDOUT, timeout=timeout).returncode
            except subprocess.TimeoutExpired:
                rc = 'timeout'
        dt = time.time() - t0
        new = [f for f, m in outputs(d).items() if before.get(f) != m]
        records.append((d, cmd, rc, dt, new))
        print(f"[{name}] {d[:2]} {cmd}: rc={rc} {dt:.0f}s -> {', '.join(new) or '-'}", flush=True)
    return records


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--timeout', type=float, default=3 * 3600)
    ap.add_argument('--only', choices=list(CHAINS))
    a = ap.parse_args()
    LOGS.mkdir(parents=True, exist_ok=True)
    chains = {a.only: CHAINS[a.only]} if a.only else CHAINS
    with ThreadPoolExecutor(len(chains)) as ex:
        results = dict(zip(chains, ex.map(lambda k: run_chain(k, chains[k], a.timeout), chains)))
    import sympy, numpy, scipy
    lines = [f"python {platform.python_version()}  sympy {sympy.__version__}  numpy {numpy.__version__}"
             f"  scipy {scipy.__version__}  platform {platform.platform()}", '']
    for name, recs in results.items():
        for d, cmd, rc, dt, new in recs:
            lines.append(f"{d}/{cmd}   rc={rc}   {dt:.0f}s")
            lines += [f"    {sha(ROOT / d / f)}  {d}/{f}" for f in sorted(new)]
    (OUT / (f'MANIFEST_{a.only}.txt' if a.only else 'MANIFEST.txt')).write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
