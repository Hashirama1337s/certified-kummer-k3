#!/usr/bin/env python
"""run_chain.py -- re-run the certificate chain in the documented order and compare with the shipped results.

Order (each step reads the previous step's result file):
  1. WP2a  core constants          code/wp2/wp2a.py   --main result_wp2a.json    (about 4 minutes)
  2. WP2b  curvature inputs        code/wp2/wp2b.py   --main result_wp2b.json    (about 10 minutes)
  3. WP3   outer residual, Ricci   code/wp3/wp3_a1.py --main result_wp3.json     (about 22 minutes; the AMEND-WP3-1
                                   runner imports the sealed wp3.py unchanged and sets two refinement tolerances)
  4. WP5   assembly + core         code/wp5/wp5.py    --main result_wp5.json     (about 1 minute)

The sealed scripts are run unchanged (see PATCHES.md for the path lines). This runner only prepares a work folder
with the layout the scripts expect:
  <work>/wp2/  wp2a.py, wp2b.py, GO-WP2a, GO-WP2b, VERDICT-WP2a.md, VERDICT-WP2b.md
  <work>/wp3/  wp3.py, wp3_a1.py, GO-WP3, VERDICT-WP3.md
  <work>/wp5/  wp5.py, GO-WP5
A GO file is the lab's launch gate (a step refuses to start without it); here the runner writes it.
Each --main invocation stops after a time budget with exit code 3 and resumes from its state folder when run
again; the runner repeats the same command until the exit code is 0. wp5.py reads the WP2a, WP2b and WP3 result
files and the first line (seal ids) of their verdict files from K3_ROOT (set to <work> by this runner).

Usage:
  python run_chain.py                      # all four steps in a fresh folder ./run, then compare
  python run_chain.py --from-shipped       # step 4 only, from the shipped WP2a/WP2b/WP3 results (about 1 minute)
  python run_chain.py --work DIR --python "py -3.11"
Compare rule: every JSON field must be identical except timing fields (names containing 'seconds'), recorded file
paths ('source', 'verdict_file'), the recorded sha256 of wp5.py itself ('sha256_wp5') and the input-file sha256
entries under 'inputs' (they change with the timing fields of freshly produced inputs).
Exit code 0 = every produced result identical to the shipped one under that rule.
"""
import argparse, json, os, shlex, shutil, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))     # repository root
CODE, RES = os.path.join(ROOT, 'code'), os.path.join(ROOT, 'results')
IGNORE = {'source', 'verdict_file', 'sha256_wp5'}

STEPS = [  # name, folder, script, output
    ('WP2a', 'wp2', 'wp2a.py', 'result_wp2a.json'),
    ('WP2b', 'wp2', 'wp2b.py', 'result_wp2b.json'),
    ('WP3', 'wp3', 'wp3_a1.py', 'result_wp3.json'),
    ('WP5', 'wp5', 'wp5.py', 'result_wp5.json'),
]


def prepare(work):
    for sub, files, gos, verdicts in [
        ('wp2', ['wp2a.py', 'wp2b.py'], ['GO-WP2a', 'GO-WP2b'], ['VERDICT-WP2a.md', 'VERDICT-WP2b.md']),
        ('wp3', ['wp3.py', 'wp3_a1.py'], ['GO-WP3'], ['VERDICT-WP3.md']),
        ('wp5', ['wp5.py'], ['GO-WP5'], []),
    ]:
        d = os.path.join(work, sub)
        os.makedirs(d, exist_ok=True)
        for f in files:
            shutil.copy2(os.path.join(CODE, sub, f), os.path.join(d, f))
        for g in gos:
            with open(os.path.join(d, g), 'w') as fh:
                fh.write('GO written by run_chain.py (release re-run)\n')
        for v in verdicts:
            shutil.copy2(os.path.join(RES, sub, v), os.path.join(d, v))


def run_step(work, py, sub, script, out):
    d = os.path.join(work, sub)
    env = dict(os.environ, K3_ROOT=work)
    for attempt in range(1, 200):
        cmd = py + [script, '--main', out]
        print('[run_chain] %s (invocation %d): %s' % (sub, attempt, ' '.join(cmd)), flush=True)
        rc = subprocess.call(cmd, cwd=d, env=env)
        if rc == 0:
            return
        if rc != 3:
            sys.exit('[run_chain] %s/%s stopped with exit code %d' % (sub, script, rc))
    sys.exit('[run_chain] %s/%s did not finish' % (sub, script))


def diff(a, b, path=''):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k in IGNORE or 'seconds' in k or (k == 'sha256' and path.startswith('/inputs/')):
                continue
            if k not in a or k not in b:
                out.append(path + '/' + k + ' (missing on one side)')
            else:
                out += diff(a[k], b[k], path + '/' + k)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(path + ' (length %d vs %d)' % (len(a), len(b)))
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                out += diff(x, y, '%s[%d]' % (path, i))
    elif a != b:
        out.append('%s: %r vs %r' % (path, a, b))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--work', default=os.path.join(ROOT, 'run'))
    ap.add_argument('--python', default=sys.executable, help='interpreter command, e.g. "py -3.11"')
    ap.add_argument('--from-shipped', action='store_true', help='run WP5 only, from the shipped WP2a/WP2b/WP3 results')
    ap.add_argument('--compare-only', action='store_true', help='compare an existing work folder, run nothing')
    a = ap.parse_args()
    work = os.path.abspath(a.work)
    py = [a.python] if os.path.exists(a.python) else shlex.split(a.python, posix=(os.name != 'nt'))
    steps = STEPS
    if a.from_shipped:
        steps = STEPS[3:]
    if not a.compare_only:
        prepare(work)
        if a.from_shipped:
            for sub, f in [('wp2', 'result_wp2a.json'), ('wp2', 'result_wp2b.json'), ('wp3', 'result_wp3.json')]:
                shutil.copy2(os.path.join(RES, sub, f), os.path.join(work, sub, f))
        for name, sub, script, out in steps:
            run_step(work, py, sub, script, out)
    ok = True
    for name, sub, script, out in steps:
        mine = json.load(open(os.path.join(work, sub, out)))
        ref = json.load(open(os.path.join(RES, sub, out)))
        d = diff(mine, ref)
        print('[run_chain] %s: %s' % (name, 'IDENTICAL' if not d else 'DIFFERENT (%d fields)' % len(d)))
        for line in d[:10]:
            print('    ' + line)
        ok = ok and not d
    print('[run_chain] ALL IDENTICAL' if ok else '[run_chain] DIFFERENCES FOUND')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
