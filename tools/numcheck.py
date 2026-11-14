#!/usr/bin/env python
"""numcheck.py -- trace every number printed in the write-ups to a result JSON, with the rounding direction checked.

  python tools/numcheck.py table            print the 30-row certified table (Markdown and LaTeX), rounded up from JSON
  python tools/numcheck.py check FILE ...   scan FILEs; every numeric token must be a traced claim or an allowed token

Rules. An upper bound must be the JSON value rounded UP at the printed precision; a lower bound rounded DOWN;
'exact' values must equal the JSON value rounded to nearest at the printed precision (or be an exact constant).
A token is allowed without a JSON source only if it is a year, an arXiv id or DOI fragment, a reference locator
(section, theorem, equation, page), a small count (<= 64), or listed in ALLOW with its reason.
Exit code 0 = every claim passes and no untraced token remains. CLAIMS.md is written by `check`.
"""
import json, math, os, re, sys
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, 'results')


def J(path):
    return json.load(open(os.path.join(R, path), encoding='utf-8'))


W5, W2A, W2B, W3 = J('wp5/result_wp5.json'), J('wp2/result_wp2a.json'), J('wp2/result_wp2b.json'), J('wp3/result_wp3.json')
REP1, REP2, REP3 = J('rep/rep1_result.json'), J('rep/rep2_result.json'), J('rep/rep3_result.json')
X4 = J('wp5/controls_wp5.json')['X4']['rows']['0.03']
S2B = J('stage2/result_s2b.json')
B = W5['stageB']


def iv(k):
    return B[k]


def maxrow(field):
    return max(r[field] for r in W2B['rows'])


C6_UP = Decimal(W2A['results']['c6']['upper'])
P_UP = Decimal(W2A['results']['P_cinf']['upper'])
R1_UP = Decimal(repr(W2A['results']['R_cinf']['upper_c6=1']))

# (printed token, direction, value, source description)
CLAIMS = [
    # theorem, section 0 of the chain
    ('2275.00', 'up', W5['leading']['(0,1e-4]']['d2_over_a3'], 'wp5 leading/(0,1e-4]/d2_over_a3'),
    ('2326.6', 'up', iv('[0, 0.001]')['d2_over_a3'], 'wp5 stageB[0,0.001]/d2_over_a3'),
    ('2441.3', 'up', iv('[0.002, 0.003]')['d2_over_a3'], 'wp5 stageB[0.002,0.003]/d2_over_a3'),
    ('2851.0', 'up', iv('[0.009, 0.01]')['d2_over_a3'], 'wp5 stageB[0.009,0.01]/d2_over_a3'),
    ('4517.6', 'up', iv('[0.029, 0.03]')['d2_over_a3'], 'wp5 stageB[0.029,0.03]/d2_over_a3'),
    ('2.327e-6', 'up', iv('[0, 0.001]')['d2'], 'wp5 stageB[0,0.001]/d2'),
    ('6.592e-5', 'up', iv('[0.002, 0.003]')['d2'], 'wp5 stageB[0.002,0.003]/d2'),
    ('2.851e-3', 'up', iv('[0.009, 0.01]')['d2'], 'wp5 stageB[0.009,0.01]/d2'),
    ('3.876e-3', 'up', iv('[0.01, 0.011]')['d2'], 'wp5 stageB[0.01,0.011]/d2'),
    ('0.121974', 'up', iv('[0.029, 0.03]')['d2'], 'wp5 stageB[0.029,0.03]/d2'),
    ('0.122', 'up', iv('[0.029, 0.03]')['d2'], 'wp5 stageB[0.029,0.03]/d2 (credit paragraph)'),
    ('1.000003', 'up', iv('[0, 0.001]')['smallness'], 'wp5 stageB[0,0.001]/smallness'),
    ('1.000081', 'up', iv('[0.002, 0.003]')['smallness'], 'wp5 stageB[0.002,0.003]/smallness'),
    ('1.003501', 'up', iv('[0.009, 0.01]')['smallness'], 'wp5 stageB[0.009,0.01]/smallness'),
    ('1.1694', 'up', W5['summary']['smallness_max'], 'wp5 summary/smallness_max'),
    ('1.169325', 'up', W5['summary']['smallness_max'], 'wp5 summary/smallness_max'),
    ('2269.27', 'up', W5['summary']['leading_coefficient_upper'], 'wp5 summary/leading_coefficient_upper (a = 0 value)'),
    ('1.015342', 'up', W5['leading']['a=0']['D'], 'wp5 leading/a=0/D'),
    ('17226.21', 'up', W3['results']['C_out_upper'], 'wp3 results/C_out_upper (= C_q)'),
    ('4.530841', 'up', W5['constants']['Gh0_upper'], 'wp5 constants/Gh0_upper'),
    ('14416.38', 'up', W5['summary']['C_core_F_max'], 'wp5 summary/C_core_F_max'),
    ('1,731,732.5', 'up', W5['summary']['C_core_R_max'], 'wp5 summary/C_core_R_max'),
    ('1695936.40', 'up', iv('[0, 0.001]')['input_sources']['CR_used'], 'wp5 stageB[0,0.001]/input_sources/CR_used'),
    ('2.0000353', 'up', maxrow('a2_C0_upper'), 'wp2b max rows/a2_C0_upper'),
    ('1.99999862', 'down', W2B['lower_bound_a003']['a2_C0_lower'], 'wp2b lower_bound_a003/a2_C0_lower'),
    ('8.34e-5', 'up', maxrow('e'), 'wp2b max rows/e'),
    ('3.27e-3', 'up', maxrow('a2_dR'), 'wp2b max rows/a2_dR'),
    ('1.37e-4', 'up', iv('[0.029, 0.03]')['input_sources']['e_used'], 'wp5 stageB/input_sources/e_used'),
    ('1.77518', 'up', iv('[0.029, 0.03]')['input_sources']['D'], 'wp5 stageB[0.029,0.03]/input_sources/D'),
    ('65064.3', 'up', W3['results']['R_out_upper'], 'wp3 results/R_out_upper'),
    ('16918.28', 'down', W3['results']['F8bulk_lower']['abs_lower'], 'wp3 results/F8bulk_lower/abs_lower'),
    ('14123.85', 'up', 3 * C6_UP / 16 * P_UP, 'wp2a: (3 c6/16) P* from c6 upper and P_cinf upper'),
    ('1695905.922', 'up', R1_UP * C6_UP, 'wp2a: R_cinf upper_c6=1 x c6 upper'),
    ('15476.93654158358624852', 'up', C6_UP, 'wp2a results/c6/upper (centre of the printed ball; radius 4e-19 covers the interval)'),
    ('4.86705830003957', 'down', W2A['results']['P_cinf']['lower'], 'wp2a results/P_cinf/lower'),
    ('4.86705830003983', 'up', W2A['results']['P_cinf']['upper'], 'wp2a results/P_cinf/upper'),
    ('0.002851', 'up', iv('[0.009, 0.01]')['d2'], 'wp5 stageB[0.009,0.01]/d2 (plain-words form)'),
    ('2.0000352762', 'up', maxrow('a2_C0_upper'), 'wp2b rows[0.025,0.03]/a2_C0_upper (zone Z1 sets it)'),
    ('1.98843008', 'up', '1.9884300760979716', 'wp2b zone Z2 value on [0.025,0.03] (literal in result_wp2b.json)'),
    ('0.70784651', 'up', '0.7078465004722245', 'wp2b zone NB value on [0.025,0.03] (literal in result_wp2b.json)'),
    ('16918.37305985292092614861', 'exact', REP2['enclosure_z_quarter']['intersection_hi'].lstrip('-'),
     'rep2 enclosure_z_quarter/intersection (absolute value of the centre)'),
    ('16918.37306', 'exact', REP2['sampled_sup']['min_F'].lstrip('-'), 'rep2 sampled_sup/min_F (absolute value)'),
    # Stage 2 (sharpened bound; result_s2b.json, certified values bit-identical to result_s2.json)
    ('765.75', 'up', S2B['summary']['sup_d2_over_a3_0_1em4'], 'stage2 result_s2b summary/sup_d2_over_a3_0_1em4'),
    ('763.14', 'up', S2B['summary']['limsup_d2_over_a3'], 'stage2 result_s2b summary/limsup_d2_over_a3 (a = 0 value)'),
    ('8.7221e-4', 'up', S2B['summary']['d2_0p01'], 'stage2 result_s2b summary/d2_0p01'),
    ('2.85599e-2', 'up', S2B['summary']['d2_0p03'], 'stage2 result_s2b summary/d2_0p03'),
    ('0.0285599', 'up', S2B['summary']['d2_0p03'], 'stage2 result_s2b summary/d2_0p03 (decimal form)'),
    # blind replication
    ('1.695996527e6', 'exact', REP3['a=0.03']['refined_max_Ric_over_a6'], 'rep3 a=0.03/refined_max_Ric_over_a6'),
    ('14125.12', 'up', REP3['a=0.03']['refined_max_absF_over_a8'], 'rep3 a=0.03/refined_max_absF_over_a8 (1.412511755e4)'),
    ('1.746877e6', 'exact', REP3['rigorous_window']['Ric_over_a6_upper_bound'], 'rep3 rigorous_window/Ric_over_a6_upper_bound'),
    ('1,695,994', 'exact', X4['sampled_max_Ric_over_a6'], 'wp5 controls X4 a=0.03 sampled_max_Ric_over_a6'),
    ('14125.11', 'exact', X4['sampled_max_F_over_a8'], 'wp5 controls X4 a=0.03 sampled_max_F_over_a8'),
    ('2269.2912', 'exact', float(REP1['limit_d2_over_a3_upper']['LemmaG']), 'rep1 limit_d2_over_a3_upper/LemmaG'),
    ('0.12197295', 'exact', float(REP1['d2_upper']['a0.03_G_rule']), 'rep1 d2_upper/a0.03_G_rule'),
    ('1.1612989', 'exact', float(REP1['smallness_a0.03']['G_rule']['M2(2,2kappa)e^eta_upper']), 'rep1 smallness_a0.03/G_rule'),
]
REP2_P3 = None
for k, v in REP2.items():
    if isinstance(v, str) and '16918.37305985' in v:
        REP2_P3 = (k, v)
if REP2_P3 is None:
    for k, v in REP2.items():
        if isinstance(v, dict):
            for k2, v2 in v.items():
                if isinstance(v2, str) and '16918.37305985' in v2:
                    REP2_P3 = (k + '/' + k2, v2)
if REP2_P3:
    CLAIMS.append(('16918.37306', 'exact', 16918.37305985292, 'rep2 %s = %s' % REP2_P3))

# tokens allowed without a JSON source, with the reason
ALLOW = {
    '0.03': 'a0 = 0.03 and the top of the a-range (wp5 a0 field, wp3 A0)',
    '0.001': 'a-interval width (wp5 config/width)',
    '1e-4': 'the interval (0, 1e-4] of the stated bound (wp5 leading key)',
    '0.01': 'a value of a',
    '0.011': 'a-interval endpoint',
    '0.009': 'a-interval endpoint',
    '0.002': 'a-interval endpoint',
    '0.003': 'a-interval endpoint',
    '0.029': 'a-interval endpoint',
    '1e-6': 'a value of a (REFEREE-WP5 item 2, quoted in chain section 0)',
    '2269.2701': 'exact bound at a = 1e-6 quoted from REFEREE-WP5 item 2 via chain section 0 (shows limsup is not a bound at a > 0)',
    '0.5': 'Lambda0 = 2 / bootstrap constant 1/2 (definition)',
    '1/2': 'definition',
    '2^-2260': 'interval endpoint of Fu-Liu-Xiao-Yau, quoted in the credit paragraph (GATE-FINAL, OPENED)',
    '4.5871': '(27e/16) closed form, exact arithmetic: 27*e/16 = 4.58710...',
    '4.5308': 'double-precision CHECK value of chain 3.4 (not a bound)',
    '128': 'arb precision bits (wp2a arb_prec_bits)',
    '160': 'arb precision bits (wp5 ctx.prec, code line 27)',
    '1152': '|W(F4)|, STANDARD, enumeration CHECK',
    '192': '|G_J|, enumeration CHECK',
    '1/4': 'window radius a/4 (definition)',
    'a/4': 'window (definition)',
    '1/32': 'covolume of L (definition)',
    '62,310': 'wp2b bulk/boxes_final',
    '25,477': 'wp3 bulk/n_eval',
    '22,601': 'wp3 bulk/n_leaves',
    '2001': 'wp5 s-grid size (chain 3.4)',
    '4000': 'wp5 config/ghN',
    '256': 'wp5 config/nT',
    '1.0000000023': 'min ratio, VERDICT-WP5 V5 and REFEREE-WP5 t2chk (chain section 9)',
    '13': 'Lean theorem count, K3Chain.lean',
    '220': 'clean-room WP2a runtime in seconds (CLEANROOM-REPORT)',
    '59': 'clean-room WP5 runtime in seconds (CLEANROOM-REPORT)',
    '0.0702': 'flat-torus sharp value, chain 3.6 CHECK',
    '2.0': 'licence version / Lambda0',
    '4.0': 'CC BY 4.0',
    '1.0': 'version 1.0',
    '1.2.0': 'cff-version',
    '3.11': 'Python version',
    '3.11.9': 'Python version tested',
    '0.9.0': 'python-flint version',
    '1.3.0': 'mpmath version',
    '2.4.4': 'numpy version',
    '4.36.0-rc2': 'Lean toolchain version',
    '1e-8': 'order of the comparator shortfall in BAR-S2 (VERDICT-S2: 6.8e-9 and 1.1e-8 relative for d2, 1.8e-12 for U)',
    '2e-12': 'order of the U comparator shortfall in BAR-S2 (VERDICT-S2: 1.8e-12 relative)',
    '70': 'upper end of the observed full-chain run time in minutes (release re-run 67, clean room 36)',
    '67': 'release re-run wall time in minutes (00:11:43 to 01:18:30 log stamps)',
    '0.3': 'percent form of d2(0.01) <= 2.851e-3 (0.2851 percent, stated as "under 0.3 percent")',
    '4e-19': 'radius of the printed c6 ball; covers the certified c6 interval (wp2a results/c6)',
    '4.88e-21': 'radius of the rep2 enclosure (rep2 enclosure_z_quarter/intersection)',
    '17/16': 'constant in the diameter lemma L5 (definition)',
    '2.0036': 'cruder global-e, global-dR value of a^2 C0 on [0.025, 0.03] (2.0035992), referee report on Lemma A7 item 5; not used',
    '2.851030e-3': 'independent float re-implementation by a referee (chain section 10)',
    '0.1219727': 'independent float re-implementation by a referee (chain section 10)',
    '2269.266': 'independent float re-implementation by a referee (chain section 10)',
    '696': 'Lemma E constant kappa(1,4), exact (referee report on the lemmas)',
    '8400': 'Lemma E constant kappa(1,5), exact (referee report on the lemmas)',
    '1e-9': 'relative slack of libm powers in WP3 (Definitions section 7 item 6)',
    '0.047': 'an earlier headline value of a0, quoted only as withdrawn (chain section 0, historical note)',
}
LOCATOR = re.compile(r'(Section|section|Sec\.|sec\.|§|Thm|Theorem|Prop\.?|Proposition|Lemma|Lemmas|eq\.|eqs\.|Eq\.|equation|'
                     r'Remark|Cor\.?|item|items|p\.|pp\.|version|Version|v|Table|Appendix|chain|DEFINITIONS|Definitions|'
                     r'step|Step|Stage|stage|ASPM|Math\.|Acta|AGAG|Geom\.|vol\.|L|S|A|I|GAP|Fix|W|P|Q|U|V|B|D|K|R|X|E|F|Z)'
                     r'[\s~]*\(?$')
NUM = re.compile(r'(?<![\w.^/-])(\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+(?:\.\d+)*(?:e-?\d+)?(?:-rc\d)?|\d+/\d+|2\^-\d+|a/4)(?![\w])')


def dec(value):
    return Decimal(value) if isinstance(value, (str, Decimal)) else Decimal(repr(float(value)))


def rounded(value, printed, mode):
    s = printed.replace(',', '')
    d = Decimal(s)
    exp = d.as_tuple().exponent
    q = Decimal(1).scaleb(exp)
    v = dec(value)
    rnd = {'up': ROUND_CEILING, 'down': ROUND_FLOOR, 'exact': ROUND_HALF_EVEN}[mode]
    if 'e' in s.lower():
        mant, e = s.lower().split('e')
        md = Decimal(mant)
        q = Decimal(1).scaleb(md.as_tuple().exponent + int(e))
    return v.quantize(q, rounding=rnd) == d.quantize(q) or (mode == 'up' and d >= v and d - v < q)


def check_claims():
    rows, bad = [], 0
    for tok, mode, val, src in CLAIMS:
        ok = rounded(val, tok, mode)
        d = Decimal(tok.replace(',', ''))
        if mode == 'up':
            ok = ok and d >= dec(val)
        if mode == 'down':
            ok = ok and d <= dec(val)
        bad += not ok
        rows.append((tok, mode, repr(val), src, 'OK' if ok else 'WRONG'))
    return rows, bad


def scan(files):
    known = {c[0] for c in CLAIMS}
    unknown = []
    for f in files:
        txt = open(f, encoding='utf-8').read()
        if f.endswith('.tex'):                                        # LaTeX number notation -> plain
            txt = re.sub(r'(\d+(?:\.\d+)?)\\times10\^\{(-?\d+)\}', r'\1e\2', txt)
            txt = re.sub(r'(?<![\d.])10\^\{(-?\d+)\}', r'1e\1', txt)
            txt = txt.replace('{,}', ',').replace('--', '-').replace('~', ' ').replace('\\ ', ' ').replace('$', ' ')
            txt = re.sub(r'\\tfrac12|\\frac12', '1/2', txt)
            txt = re.sub(r'\\(le|ge)(?=\d)', ' <= ', txt)
            txt = re.sub(r'\((\d\.\d)\)', r'(Section \1)', txt)          # "(1.3)" section references
        if f.endswith('.json'):
            txt = txt.replace('&le;', '<=')
        txt = re.sub(r'https?://\S+', ' ', txt)                     # URLs
        txt = re.sub(r'10\.\d{4,5}/[^\s)\]]+', ' ', txt)               # DOIs
        txt = re.sub(r'\b\d{4}\.\d{5}(v\d)?\b', ' ', txt)              # arXiv ids
        txt = re.sub(r'\b(?=[0-9a-f]*[a-f])[0-9a-f]{10}\b', ' ', txt)  # seal ids (always contain a letter)
        txt = re.sub(r'`?[0-9a-f]{64}`?', ' ', txt)                   # hashes
        for ln, line in enumerate(txt.split('\n'), 1):
            for m in NUM.finditer(line):
                tok = m.group(1)
                pre = line[max(0, m.start() - 14):m.start()]
                if tok in known or tok in ALLOW:
                    continue
                if re.fullmatch(r'\s*#+\s*', line[:m.start()]) or re.search(r'(Chain|Definitions|Paper|\|)\s*$', pre) and re.fullmatch(r'\d+(\.\d+)*', tok) and len(tok) <= 4:
                    continue                                              # headings and section cross-references
                if re.fullmatch(r'\d\.\d', tok) and re.search(r'(,|-|definitions|interfaces|Sec|\|)\s*$', pre):
                    continue                                              # section lists such as "2, 3.1" or "1.1-1.6"
                if re.search(r'(arXiv:|\(\d{4}\),)\s*$', pre):
                    continue                                              # old-style arXiv id, page range after a year
                if re.fullmatch(r'0\.0[0-3]\d?', tok) and re.search(r'[\[(,]\s*$', pre):
                    continue                                              # a-interval endpoints inside [ , ]
                if re.fullmatch(r'(19|20)\d\d', tok):
                    continue
                if re.fullmatch(r'\d+', tok) and int(tok) <= 64:
                    continue
                if LOCATOR.search(pre.rstrip()):
                    continue
                if re.fullmatch(r'\d+(\.\d+){1,3}', tok) and re.search(r'(Section|section|§|item|Thm|Lemma|L\d|S\d)\s*[\w.]*\s*(,|and|-)?\s*$', pre):
                    continue
                unknown.append((os.path.relpath(f, ROOT), ln, tok, line.strip()[:110]))
    return unknown


def table():
    keys = list(B.keys())
    md = ['| a-interval | U (continuity step) <= | d2 <= | d2/a^3 <= |', '|---|---|---|---|']
    tex = []
    out = []
    for k in keys:
        r = B[k]
        U = Decimal(repr(r['smallness'])).quantize(Decimal('1e-6'), rounding=ROUND_CEILING)
        d2 = Decimal(repr(r['d2']))
        e = d2.adjusted()
        d2u = d2.quantize(Decimal(1).scaleb(e - 3), rounding=ROUND_CEILING)
        d2s = ('%.3fe%d' % (float(d2u.scaleb(-e)), e)) if e < -1 else str(d2u)
        da = Decimal(repr(r['d2_over_a3'])).quantize(Decimal('0.1'), rounding=ROUND_CEILING)
        md.append('| %s | %s | %s | %s |' % (k, U, d2s, da))
        tex.append('$%s$ & %s & %s & %s \\\\' % (k, U, d2s.replace('e', '\\times10^{') + '}' if 'e' in d2s else d2s, da))
        out.append((d2s, 'up', r['d2']))
        out.append((str(U), 'up', r['smallness']))
        out.append((str(da), 'up', r['d2_over_a3']))
    return md, tex, out


TABLE_CLAIMS = table()[2]
for tok, mode, val in TABLE_CLAIMS:
    if tok not in {c[0] for c in CLAIMS}:
        CLAIMS.append((tok, mode, val, 'wp5 stageB (30-row table)'))


def main():
    if len(sys.argv) >= 2 and sys.argv[1] == 'table':
        md, tex, _ = table()
        print('\n'.join(md)); print(); print('\n'.join(tex)); return 0
    if len(sys.argv) < 3 or sys.argv[1] != 'check':
        sys.exit(__doc__)
    rows, bad = check_claims()
    unknown = scan(sys.argv[2:])
    L = ['# CLAIMS: every number in the write-ups, traced to a result file', '',
         'Generated by `tools/numcheck.py check`. Direction: `up` = printed upper bound is the JSON value rounded up at the',
         'printed precision; `down` = lower bound rounded down; `exact` = rounded to nearest (replication values, quoted).', '',
         '| Printed | Direction | JSON value | Source (results/...) | Check |', '|---|---|---|---|---|']
    for tok, mode, val, src, st in rows:
        L.append(f'| {tok} | {mode} | {val} | {src} | {st} |')
    L += ['', 'Tokens allowed without a JSON source (reason):', '']
    for k, v in ALLOW.items():
        L.append(f'- `{k}`: {v}')
    open(os.path.join(ROOT, 'CLAIMS.md'), 'w', encoding='utf-8', newline='\n').write('\n'.join(L) + '\n')
    print('claims: %d checked, %d wrong' % (len(rows), bad))
    for r in rows:
        if r[4] != 'OK':
            print('  WRONG', r)
    print('untraced tokens: %d' % len(unknown))
    for u in unknown:
        print('  %s:%d  %s   | %s' % u)
    return 1 if (bad or unknown) else 0


if __name__ == '__main__':
    sys.exit(main())
