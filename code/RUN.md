# RUN: re-run the certificate chain

Requirements: Python 3.11 (tested with 3.11.9), python-flint 0.9.0, numpy (tested with 2.4.4), mpmath (tested with
1.3.0). The lab ran every step under an 8 GB memory cap; the largest step (WP3) peaked near 3.3 GB.

## One command

From the repository root:

```
python code/run_chain.py --from-shipped     # WP5 only, from the shipped WP2a/WP2b/WP3 results: about 1 minute
python code/run_chain.py                    # the whole chain from empty state: about 35 to 70 minutes on one core
```

Expected last line: `[run_chain] ALL IDENTICAL`. The runner works in `./run` (or `--work DIR`); pass another
interpreter with `--python "py -3.11"`.

## The order, and what the runner does

| Step | Command (inside the work folder) | Reads | Writes |
|---|---|---|---|
| 1. WP2a core constants | `wp2/wp2a.py --main result_wp2a.json` | nothing | `wp2/result_wp2a.json` |
| 2. WP2b curvature inputs | `wp2/wp2b.py --main result_wp2b.json` | imports `wp2a.py` | `wp2/result_wp2b.json` |
| 3. WP3 outer residual and Ricci | `wp3/wp3_a1.py --main result_wp3.json` | imports the sealed `wp3.py` | `wp3/result_wp3.json` |
| 4. WP5 full core and assembly | `wp5/wp5.py --main result_wp5.json` | the three results above and the first line of their verdict files, from `K3_ROOT` | `wp5/result_wp5.json` |

- **GO files.** Each `--main` refuses to start unless its launch file (`GO-WP2a`, `GO-WP2b`, `GO-WP3`, `GO-WP5`)
  exists next to the script. In the lab the lead wrote it after reading the controls; the runner writes it.
- **Exit code 3 = resume.** Each invocation stops after a time budget (about 9 minutes) with exit code 3 and keeps its
  state in a `state*` folder; running the same command again continues. The runner repeats until exit code 0.
- **WP3 through the amendment runner.** `wp3_a1.py` imports the sealed `wp3.py` unchanged and sets only the two
  refinement-stopping tolerances of AMEND-WP3-1 (`results/wp3/AMEND-WP3-1.md`), as in the certified run.
- **Input paths.** `wp5.py` reads its inputs from `K3_ROOT` (an environment variable; default: the parent of its own
  folder). The runner sets `K3_ROOT` to the work folder, so no lab path is needed. This is one of four path lines
  changed from the lab copies; all are listed in [PATCHES.md](PATCHES.md).
- **Compare.** Every produced JSON is compared field by field with `results/`; ignored are timing fields, recorded
  file paths (`source`, `verdict_file`), the recorded hash of `wp5.py` itself (the shipped file differs from the
  sealed one in its path lines) and the input-file hashes WP5 records (freshly produced inputs differ from the shipped
  ones only in their timing fields; every other field of those inputs is compared directly).

## Controls

The control modes (`--controls OUT.json`) of WP2a, WP2b and WP5 need lab reference files that are not in this
repository (the double-precision probe scripts and two earlier reference outputs); their recorded outputs are in
`results/` (`controls.json`, `controls_wp2b_sealed.json`, `controls_wp3.json`, `controls_wp5.json`), each with
`ALL_CONTROLS_PASS` true.

## Other scripts

- `wp2/r6sym.py`: exact symbolic check of the R6 angular reduction (needs sympy; log in `results/wp2/r6sym.log`).
- `rep/rep1/rep1.py`, `rep/rep2/rep2.py`, `rep/rep3/rep3.py`: the three blind replication seats (rep2 also needs
  scipy); results in `results/rep/`.
- `stage2/s2.py`, `stage2/s2b.py`: the Stage 2 sharpened bound (not the main theorem). They read the WP2a, WP2b, WP3
  and WP5 results from `K3_ROOT` (default: the parent of their folder), like wp5.py; results and verdicts in
  `results/stage2/`.

## Figures

The two README figures are drawn from shipped files only: `python figures/make_figures.py .` from the repository root
(Python 3.11, matplotlib 3.11; seconds; reads only shipped files). It rewrites `figures/*.svg`; two runs give
byte-identical SVGs.
