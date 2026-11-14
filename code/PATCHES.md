# PATCHES: shipped files that differ from the lab originals

Only absolute lab paths were changed: `<lab>` stands for the lab project folder (the layout of `code/` and
`results/` mirrors it), `<lab-tools>` for the lab tool folder, `<scratch>` for a session scratch folder.
Four lines in code are functional path lines (marked FUNCTIONAL); every other change is in a comment, a docstring,
a recorded file path inside a result JSON, or prose. No number, formula or other logic was changed.

Seal ids in the verdict files refer to the lab copies (lab sha256 below); the shipped sha256 is that of the file
in this repository. The certified run of WP5 records the sha256 of the lab `wp5.py`; a re-run with the shipped file
records the shipped hash instead (the runner `code/run_chain.py` ignores that field when comparing).

## code/wp2/wp2a.py (lab: wp2/wp2a.py)
- lab sha256: `22b2bd300e0a561849809e7636f07b0a7ee52caebb5a1089e0d3479caa55a75e`
- shipped sha256: `375283b863faebfcd9f8c61ca0978272bcc0dfff8da0ba461f5648af71e095c1`
- line 16: path wording only; now `--main OUT.json       the certified run; refuses unless <lab>/wp2/GO-WP2a exists; chunked and`
- line 17: path wording only; now `resumable (state in <lab>/wp2/state/); rerun until status == COMPLETE`
- line 18: path wording only; now `Run:  python <lab-tools>/memcap.py 8 python <lab>/wp2/wp2a.py --controls <lab>/wp2/controls.json`
- line 890 FUNCTIONAL: was `REF_JSON = '<lab>/wp1/v2/v2_profile.json'`
  now `REF_JSON = os.path.join(os.path.dirname(HERE), 'wp1', 'v2', 'v2_profile.json')`

## code/wp2/wp2b.py (lab: wp2/wp2b.py)
- lab sha256: `9cd117312fedb699cd6c67e7b1e7ad3f52427fbcc8851a8028be96c95362e85f`
- shipped sha256: `200f66aeeeb563a01a76ad50990907eaa1d3415f1e2a01f72676eed33af59273`
- line 20: path wording only; now `--main OUT.json       refuses unless <lab>/wp2/GO-WP2b exists; resumable (state_wp2b/), exit 3 = resume.`
- line 21: path wording only; now `Run:  python <lab-tools>/memcap.py 8 python <lab>/wp2/wp2b.py --controls <lab>/wp2/controls_wp2b.json`
- line 1091 FUNCTIONAL: was `sys.path.insert(0, '<lab>/probe1')`
  now `sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'probe1'))`
- line 1191 FUNCTIONAL: was `sys.path.insert(0, '<lab>/probe1'); import k1c`
  now `sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'probe1')); import k1c`

## code/wp3/wp3.py (lab: wp3/wp3.py)
- lab sha256: `80d5492f5a5a19b06ef1208c239e735f61c3e7897d80a6640d0392303a8f365f`
- shipped sha256: `cfda2e3e8e15fd351324459169ac862a85e5e5f29f6e1d94b7e4b11eb49a92c0`
- line 22: path wording only; now `--main OUT.json       certified run; refuses unless <lab>/wp3/GO-WP3 exists; chunked (exit 3 = resume)`
- line 23: path wording only; now `Run:  python <lab-tools>/memcap.py 8 python <lab>/wp3/wp3.py --controls <lab>/wp3/controls_wp3.json`
- line 2568: path wording only; now `REF = 14569.136077584635          # <lab>/wp1/check_num.json 'bulk_cell'[0] (full precision)`

## code/wp5/wp5.py (lab: wp5/wp5.py)
- lab sha256: `b3919b17955766b501d0b1e7f0d1113a3c3e0cb074f1d5e38546bc416067287a`
- shipped sha256: `e0182e9ae54f8f73de5a4fc01974110329c7dd4d8c65bff23991537a5142ff38`
- line 15: path wording only; now `--main OUT.json       refuses unless <lab>/wp5/GO-WP5 exists; resumable (state_wp5/), exit 3 = rerun.`
- line 16: path wording only; now `Run:  python <lab-tools>/memcap.py 8 py -3.11 <lab>/wp5/wp5.py --controls <lab>/wp5/controls_wp5.json`
- line 23 FUNCTIONAL: was `K3 = '<lab>'`
  now `K3 = os.environ.get('K3_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))`

## code/rep/rep2/rep2.py (lab: rep/rep2/rep2.py)
- lab sha256: `71714dffb36220e7938ab54c9f2894f765b3e8cea60e85654a8a46e72b5b2142`
- shipped sha256: `fc6df51c95a7b22230f29bb0b8bccbb64e5d35fa40e94fe2904629254d0e452a`
- line 19: path wording only; now `Run:  python <lab-tools>/memcap.py 8 py -3.11 rep2.py`

## results/wp5/result_wp5.json (lab: wp5/result_wp5.json)
- lab sha256: `bcfe3a7581cad50daabf177d30901bd8f9d9eb87720b1c204281c0e39c447012`
- shipped sha256: `c4e4de3ea759a8a97a88bf8ab271304ca210eb31d9019b94ed5abc3bf31248f3`
- line 13: path wording only; now `"source": "<lab>/wp2/result_wp2a.json",`
- line 15: path wording only; now `"verdict_file": "<lab>/wp2/VERDICT-WP2a.md",`
- line 27: path wording only; now `"source": "<lab>/wp2/result_wp2b.json",`
- line 29: path wording only; now `"verdict_file": "<lab>/wp2/VERDICT-WP2b.md",`
- line 41: path wording only; now `"source": "<lab>/wp3/result_wp3.json",`
- line 43: path wording only; now `"verdict_file": "<lab>/wp3/VERDICT-WP3.md",`

## results/wp5/controls_wp5.json (lab: wp5/controls_wp5.json)
- lab sha256: `cb31f2e833482a86aaa1632eab81afee17730e378b5843ee1c8b384d284eb53e`
- shipped sha256: `d406a5e0a58589d3945b92477ba769d731a2e20458c6770d10242613513194b5`
- line 12: path wording only; now `"source": "<lab>/wp2/result_wp2a.json",`
- line 14: path wording only; now `"verdict_file": "<lab>/wp2/VERDICT-WP2a.md",`
- line 36: path wording only; now `"source": "<lab>/wp2/result_wp2b.json",`
- line 38: path wording only; now `"verdict_file": "<lab>/wp2/VERDICT-WP2b.md",`
- line 106: path wording only; now `"source": "<lab>/wp3/result_wp3.json",`
- line 108: path wording only; now `"verdict_file": "<lab>/wp3/VERDICT-WP3.md",`

## results/wp5/VERDICT-WP5.md (lab: wp5/VERDICT-WP5.md)
- lab sha256: `1b8926d668503fc95ce56c5172147e938d00ce5287ace3ca95e58a7632b2c0c7`
- shipped sha256: `e1dfa7fc248cb2605b537acd3b683453bb87bae883a0143356e4a92e92470050`
- line 24: path wording only; now `V5 note (disclosed): the lead's FIRST check (v5check.py) used wp4/chain.py's double-precision Gamma grid (Gamh_up rounds K D^2 up to the nex...`

## results/wp5/VERDICT-WP5.pre-erratum.md (lab: wp5/VERDICT-WP5.pre-erratum.md)
- lab sha256: `f234c8461de0af7ced0784f061ea932041e0a14c0e6a5e23a35b9a9b5a0d6438`
- shipped sha256: `1777894f80f1ac7728db76c13e4f418d7c4cac294e5ef8ade908aae595b7655f`
- line 22: path wording only; now `V5 note (disclosed): the lead's FIRST check (v5check.py) used wp4/chain.py's double-precision Gamma grid (Gamh_up rounds K D^2 up to the nex...`

## results/cleanroom/CLEANROOM-REPORT.md (lab: cleanroom/CLEANROOM-REPORT.md)
- lab sha256: `b9feeee94b13b0b79f54273a947dd9f47c947a4192c8eb76313b14c4c6c84c47`
- shipped sha256: `3fe6881254ce247c3a75b4d07eeed65622d2d4d1b2bd7502c974d5f0b9f6df9a`
- line 11: path wording only; now `Limitation (disclosed): wp5.py reads its input result files from the fixed path <lab>/wp2, wp3 (not the clean-room copies); since the clean-...`

## results/probes/BAR-K1.md (lab: probe1/BAR-K1.md)
- lab sha256: `a6396f080a0c3e7c0ba346761e29523c43c3841e3729821fca40faa59398cf4e`
- shipped sha256: `c71860dce3bbea3c1ba996bf4b09d0c706150a548c4bc13ce581d35ed0848931`
- line 36: path wording only; now `- command: python <lab-tools>/memcap.py 8 python k1.py result_k1.json ; then k1_mp.py result_k1.json c4_k1.json`

## code/stage2/s2.py (lab: stage2/s2.py)
- lab sha256: `b4b6c96c35dc93b5a3e418493ad461ce9e358a0faae5b4e1a0828a5f48728c26`
- shipped sha256: `4d70d13835473618fca023ee2d6b7eb1df04a9d87244f645f948d739db11cd52`
- line 29: path wording only; now `Run:  python <lab-tools>/memcap.py 8 py -3.11 <lab>/stage2/s2.py --controls <lab>/stage2/controls_s2.json`
- line 36 FUNCTIONAL: was `K3 = '<lab>'`
  now `K3 = os.environ.get('K3_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))`

## results/stage2/result_s2.json (lab: stage2/result_s2.json)
- lab sha256: `47cfc2cf0ff6b2c676c10cd3996b4ac69c3ab55a402f11421b77524b750938be`
- shipped sha256: `abfdbd8daeda088f1a2e310f071972b6e861b03aea40c4519869a202426c0212`
- line 16: path wording only; now `"source": "<lab>/wp2/result_wp2a.json",`
- line 18: path wording only; now `"verdict_file": "<lab>/wp2/VERDICT-WP2a.md",`
- line 30: path wording only; now `"source": "<lab>/wp2/result_wp2b.json",`
- line 32: path wording only; now `"verdict_file": "<lab>/wp2/VERDICT-WP2b.md",`
- line 44: path wording only; now `"source": "<lab>/wp3/result_wp3.json",`
- line 46: path wording only; now `"verdict_file": "<lab>/wp3/VERDICT-WP3.md",`
- line 60: path wording only; now `"source": "<lab>/wp5/result_wp5.json",`

## results/stage2/result_s2b.json (lab: stage2/result_s2b.json)
- lab sha256: `bb61e64fd21c8b30357927a3eae266f526bf25dc8259bc4c058277256f3a1736`
- shipped sha256: `fc68f1256dafb1ed2d6e5a48139d486d940f85e14153780461ee14c5d91c2f93`
- line 16: path wording only; now `"source": "<lab>/wp2/result_wp2a.json",`
- line 18: path wording only; now `"verdict_file": "<lab>/wp2/VERDICT-WP2a.md",`
- line 30: path wording only; now `"source": "<lab>/wp2/result_wp2b.json",`
- line 32: path wording only; now `"verdict_file": "<lab>/wp2/VERDICT-WP2b.md",`
- line 44: path wording only; now `"source": "<lab>/wp3/result_wp3.json",`
- line 46: path wording only; now `"verdict_file": "<lab>/wp3/VERDICT-WP3.md",`
- line 60: path wording only; now `"source": "<lab>/wp5/result_wp5.json",`

## docs/ANALYSIS-CHAIN-v2.2.md (lab: ANALYSIS-CHAIN-v2.2.md)
- lab sha256: `0027a64f7c8351a6b00c26b5bcf0ee1e22a6717f11bf79b2129bb2d3137fe43f`
- shipped sha256: `36bef341b82baab4f80e77a59afc84e68545acd2a548f90614e059e6b1d2f60b`
- line 7: path wording only; now `Scripts of v1 in <lab>/wp4/ (green_const.py, torus_ctrl.py, referee_gamma.py, chain.py, explicit.py, each with a .log); Li-Yau page images w...`
- line 157: path wording only; now `- Independent float re-implementation of the chain with the WP5 inputs (REFEREE-CHAIN, scratch folder refchain.py): d2(0.01) = 2.851030e-3, ...`

## docs/DEFINITIONS-AND-LEMMAS-v4.2.md (lab: DEFINITIONS-AND-LEMMAS-v4.2.md)
- lab sha256: `d2f02579268f3a03f481583d19c6a34334d8c27467705777bf3ba319a4cccd1b`
- shipped sha256: `e20d2c1e4ff8241346a4166aa4979fdeebfa765c72dbda2ecc85e90440630c6f`
- line 25: path wording only; now `Scripts for every CHECK line: `<lab>/wp1/check_sym.py` (exact, sympy), `check_num.py`, `check_peak.py`, `check_misc.py` (double precision); ...`
- line 518: path wording only; now `Proof. (i)-(iii) are identities of rational functions in zeta_1, zeta_2, zetabar_1, zetabar_2 (treated as independent symbols), Pi, Pi', Pi'...`

## lean/K3Chain.lean (lab: LEAN/k3lean/K3Chain.lean)
- lab sha256: `993512abdf90f728842847550c99b215d0dec32bce8aab5a7ade6af16d636261`
- shipped sha256: `61381a682009d44ee9e15064e093deab12b66b53447ea9d326af0368351ca3e7`
- line 6: path wording only; now `Source: <lab>/ANALYSIS-CHAIN.md, section 4 (Proposition C0, step S1) and`
- line 7: path wording only; now `section 6 (trace to eigenvalues, Route 2), and <lab>/seats/REFEREE-WP4.md item 4.`
