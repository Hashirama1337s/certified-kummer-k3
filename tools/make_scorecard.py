#!/usr/bin/env python
"""make_scorecard.py -- build scorecard.svg and SCORECARD.md from the verdict files in results/.

Every sealed test of this project is listed: the per-prediction verdicts are read from the scoring table of each
VERDICT file (rows whose first cell starts with a prediction id such as P1, W3, V5), and the test-level outcome
(HOLD / FAILED / VOID) is read from the lab ledger export results/LEDGER-K3.txt. FAILED and VOID rows are never
dropped or moved. Layout: the lab's house scorecard (white card, section bands, pill chips), then the neon palette
of the sibling repositories. Run from the repository root:  python tools/make_scorecard.py
"""
import os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, 'results')

# (section, test id, seal, verdict file or None, one-line description from the sealed claim)
TESTS = [
    ('PROBES THAT SHAPED THE APPROXIMANT (SAMPLED, NOT CERTIFIED)', [
        ('K1', '08b7d89f85', 'probes/VERDICT-K1.md', 'Second-order approximant: residual of order a^8'),
        ('K1b', '8a8380e6ef', None, 'Re-sealed before any run (seal fields added)'),
        ('K1b', 'ad344efc9b', None, 'Re-sealed before any run (seal format only)'),
        ('K1b', 'ff0c52846b', 'probes/VERDICT-K1b.md', 'Corrected approximant: bounded a^8 residual'),
        ('K1c', '214e4b3168', 'probes/VERDICT-K1c.md', 'Inner cut-off (0.1a): Ricci and curvature bars'),
        ('K1d', '4143022d85', 'probes/VERDICT-K1d.md', 'Cut-off outside the core: smooth approximant'),
    ]),
    ('CERTIFICATE RUNS (INTERVAL ARITHMETIC, EVERY a IN (0, 0.03])', [
        ('WP2a', '984fff60a3', None, 'Superseded before any code (path field fixed)'),
        ('WP2a', '2d19816000', 'wp2/VERDICT-WP2a.md', 'Core profile constants c6, P*, R6*'),
        ('WP2b', 'cdf666b973', 'wp2/VERDICT-WP2b.md', 'Curvature inputs e, dR, B+, C0'),
        ('WP3', '4158aba4f8', 'wp3/VERDICT-WP3.md', 'Residual and Ricci outside the cores'),
        ('WP5', '055009634d', 'wp5/VERDICT-WP5.md', 'Full core + assembly: d2(a), a0 = 0.03'),
    ]),
    ('INDEPENDENT CHECKS', [
        ('REP', '700278e11f', 'rep/VERDICT-REP.md', 'Blind replication: assembly and bulk peak'),
        ('REP3', '93b705fabc', 'rep/VERDICT-REP3.md', 'Blind replication: full core (WP5 Stage A)'),
    ]),
    ('STAGE 2: SHARPENED BOUND (NOT THE MAIN THEOREM)', [
        ('S2', 'a4c39ac007', 'stage2/VERDICT-S2.md', 'L5\' + spectral constant; float comparator'),
        ('S2b', '1e719251da', 'stage2/VERDICT-S2b.md', 'Re-test, cancellation-free comparator'),
    ]),
]
PENDING = []
# Verdicts corrected in the verdict file after the ledger entry was written: the corrected verdict is used, and the
# ledger's earlier outcome is reported next to it.
CORRECTED = {
    'a4c39ac007': ('VOID', [('HOLD', 'S1-S3 HOLD'), ('FAILED', 'S4 FAILED'), ('VOID', 'VOID')],
                   'S1 HOLD; S2 HOLD; S3 HOLD; S4 FAILED as sealed (in-run comparator cancellation), so the test is VOID '
                   'by its bar rule ("Verdict as sealed (CORRECTED ...)" in the verdict file; the ledger export still '
                   'shows the pre-correction HOLD). Certified values not in question.'),
}

ID_RE = re.compile(r'^\s*([A-Z]\d+)\b')


def ledger():
    out = {}
    for line in open(os.path.join(R, 'LEDGER-K3.txt'), encoding='utf-8'):
        m = re.match(r'([0-9a-f]{10})\s.*->\s+(HOLD|FAILED|VOID|open)\b', line)
        if m:
            out[m.group(1)] = m.group(2)
    return out


def predictions(path):
    rows = []
    for line in open(os.path.join(R, path), encoding='utf-8'):
        if not line.startswith('|'):
            continue
        cells = [c.strip() for c in line.strip().strip('|').split('|')]
        m = ID_RE.match(cells[0])
        if not m:
            continue
        last = [c for c in cells if c][-1]
        v = re.match(r'(HOLD|FAILED|PASS|VOID)', last)
        if v:
            rows.append((m.group(1), v.group(1), last))
    return rows


def compress(ids):
    """['P1','P3','P4','P5'] -> 'P1 P3-P5'"""
    out, i = [], 0
    while i < len(ids):
        j = i
        while (j + 1 < len(ids) and ids[j + 1][0] == ids[i][0]
               and int(ids[j + 1][1:]) == int(ids[j][1:]) + 1):
            j += 1
        out.append(ids[i] if j - i < 2 else ids[i] + '-' + ids[j])
        if j - i == 1:
            out.append(ids[j])
        i = j + 1
    return ' '.join(out)


# house palette (light) and the neon mapping used by the sibling repositories
FONT = 'Segoe UI, Helvetica, Arial, sans-serif'
INK, MUTED, ACCENT, BAND, LINE, BORDER = '#18192b', '#555b74', '#3a3fa8', '#eef0fa', '#eceef5', '#d9dbe8'
CHIP = {'HOLD': ('#dff3e8', '#1f7a4d'), 'PASS': ('#dff3e8', '#1f7a4d'), 'FAILED': ('#fbe3e1', '#b3322a'),
        'VOID': ('#fbf0d4', '#8a6100')}
NEON = {'#ffffff': '#1c2026', '#d9dbe8': '#2f353d', '#18192b': '#e8edf2', '#555b74': '#9aa4ae', '#eef0fa': '#262b31',
        '#3a3fa8': '#39ff14', '#eceef5': '#22272e', '#dff3e8': '#12301c', '#1f7a4d': '#39ff14', '#fbe3e1': '#3a1a1a',
        '#b3322a': '#ff6b6b', '#fbf0d4': '#33290f', '#8a6100': '#ffc94d'}


def tw(s, per=7.4):
    return round(len(s) * per)


def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def build():
    led = ledger()
    sections, md_rows, counts = [], [], {'HOLD': 0, 'FAILED': 0, 'VOID': 0}
    for head, tests in TESTS:
        rows = []
        for tid, seal, vf, desc in tests:
            outcome = led.get(seal)
            if outcome not in counts:
                sys.exit('ledger outcome missing for %s %s' % (tid, seal))
            chips, detail = [], ''
            if seal in CORRECTED:
                txt = open(os.path.join(R, vf), encoding='utf-8').read()
                if 'Verdict as sealed (CORRECTED' not in txt or 'VOID' not in txt:
                    sys.exit('correction not found in ' + vf)
                outcome, chips, detail = CORRECTED[seal]
            counts[outcome] += 1
            if seal in CORRECTED:
                pass
            elif vf:
                preds = predictions(vf)
                if not preds:
                    sys.exit('no scoring rows in ' + vf)
                for key in ('HOLD', 'PASS', 'FAILED', 'VOID'):
                    ids = [p for p, v, _ in preds if v == key]
                    if ids:
                        chips.append((key, compress(ids) + ' ' + key))
                detail = '; '.join('%s %s' % (p, v) for p, v, _ in preds)
                derived = 'FAILED' if any(v == 'FAILED' for _, v, _ in preds) else 'HOLD'
                if derived != outcome:
                    sys.exit('verdict file and ledger disagree for %s %s' % (tid, seal))
            else:
                chips.append(('VOID', 'VOID'))
                detail = 'no run (superseded before data)'
            rows.append((tid + ' ' + seal, desc, chips))
            md_rows.append((head, tid, seal, desc, outcome, detail, vf or 'results/LEDGER-K3.txt'))
        sections.append((head, rows))
    title = 'Scorecard: every settled sealed test (%d HOLD, %d FAILED, %d VOID)' % (
        counts['HOLD'], counts['FAILED'], counts['VOID'])
    sub = 'Sealed = prediction and pass bar hashed before the run. Failures are shown as failures. Details: SCORECARD.md.'
    W, y, body = 860, 66, []
    for head, rows in sections:
        body.append(f'<rect x="10" y="{y}" width="{W - 20}" height="22" fill="{BAND}"/>')
        body.append(f'<text x="20" y="{y + 16}" font-size="11.5" font-weight="600" fill="{ACCENT}" letter-spacing="0.5">{esc(head)}</text>')
        y += 22
        for tid, what, chips in rows:
            body.append(f'<text x="20" y="{y + 22}" font-size="12.5" font-family="Consolas, Menlo, monospace" fill="{INK}">{esc(tid)}</text>')
            body.append(f'<text x="160" y="{y + 22}" font-size="12.5" fill="{INK}">{esc(what)}</text>')
            x = 520.0
            for key, label in chips:
                fill, fg = CHIP[key]
                t = tw(label)
                body.append(f'<rect x="{x}" y="{y + 8}" width="{t + 20}" height="20" rx="10" fill="{fill}"/>')
                body.append(f'<text x="{x + 10}" y="{y + 22}" font-size="11.5" font-weight="600" fill="{fg}" '
                            f'textLength="{t}" lengthAdjust="spacingAndGlyphs">{esc(label)}</text>')
                x += t + 26
            y += 30
            body.append(f'<line x1="10" y1="{y + 1}" x2="{W - 10}" y2="{y + 1}" stroke="{LINE}"/>')
        y += 8
    H = y + 6
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="{FONT}" '
           f'role="img" aria-label="{esc(title)}">\n<title>{esc(title)}</title>\n'
           f'<rect width="{W}" height="{H}" rx="10" fill="#ffffff" stroke="{BORDER}"/>\n'
           f'<text x="20" y="34" font-size="18" font-weight="600" fill="{INK}">{esc(title)}</text>\n'
           f'<text x="20" y="54" font-size="12" fill="{MUTED}">{esc(sub)}</text>\n' + '\n'.join(body) + '\n</svg>\n')
    pat = re.compile('|'.join(re.escape(k) for k in NEON), re.I)
    svg = pat.sub(lambda m: NEON[m.group(0).lower()], svg)
    with open(os.path.join(ROOT, 'scorecard.svg'), 'w', encoding='utf-8', newline='\n') as f:
        f.write(svg)
    chip = {'HOLD': '![HOLD](brand/chips/hold.svg)', 'FAILED': '![FAILED](brand/chips/failed.svg)',
            'VOID': '![VOID](brand/chips/void.svg)'}
    md = ['# Scorecard', '', title + '.', '',
          'Generated by `tools/make_scorecard.py` from the verdict files in `results/` and the ledger export',
          '`results/LEDGER-K3.txt`. A sealed test is one whose predictions and pass bars were hashed before its run.', '',
          '| Section | Test | Seal | What was predicted | Verdict | Per prediction (as scored) | Source |',
          '|---|---|---|---|---|---|---|']
    for head, tid, seal, desc, outcome, detail, src in md_rows:
        sec = head.split(' (')[0].capitalize()
        srcl = src if src.startswith('results/') else 'results/' + src
        md.append(f'| {sec} | {tid} | `{seal}` | {desc} | {chip[outcome]} | {detail} | [{srcl}]({srcl}) |')
    md += ['', 'Notes:',
           '- WP2a is FAILED as sealed: its predictions W1 and W2 used rounded constants as containment targets, and the',
           '  exact enclosures do not contain the rounded values. The enclosures themselves are valid and are what the',
           '  proof uses (results/wp2/VERDICT-WP2a.md).',
           '- WP5 prediction V5 is marked POST HOC in its verdict: it passed only after the comparator was re-read.',
           '- REP, REP3: Q1-Q5 are controls of the sealed bars (any failure would have voided the test). REP3 results were',
           '  produced before its instrument was pinned (disclosed in results/rep/VERDICT-REP3.md).',
           '- The probes K1-K1d sampled the approximant in double precision; they shaped the construction and certify nothing.',
           '- All three load-bearing pieces (assembly arithmetic, bulk residual peak, full core) are blind-replicated.']
    with open(os.path.join(ROOT, 'SCORECARD.md'), 'w', encoding='utf-8', newline='\n') as f:
        f.write('\n'.join(md) + '\n')
    print(title)


if __name__ == '__main__':
    build()
