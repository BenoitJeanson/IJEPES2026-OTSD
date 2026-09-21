# OTSD campaign replay

An interactive replay of the component-ablation campaign behind the IJEPES optimal
transmission switching paper: 162 runs, 480 outer local-search phases, and **88,149
individual Benders candidates**, every one of them a topology you can look at.

**→ [Open the viewer](https://benoitjeanson.github.io/IJEPES2026-OTSD/)**

No install, no server, no build step: one HTML file and a folder of JSON. It works on a
phone — the panels become swipe-over sheets behind ☰ and ⓘ, one finger pans and two
pinch-zoom, and the whole corpus index is 4 KB gzipped (the largest single run, 17,743
candidates, is 300 KB).

---

## What you are looking at

The campaign has three nested levels, and the viewer gives each of them a control:

| level | what it is | count | where it is in the UI |
|---|---|---|---|
| run | a `(system, H, d_viol, d_sol, config, seed)` configuration | 162 | left-hand list |
| phase | one outer local-search iteration | 480 | the `i1 i2 i3` buttons |
| candidate | one Benders iteration — an open-branch set the master proposed | 88,149 | the scrubber |

### The scrub paints; it does not animate

The obvious thing to build would be an animation of the topology evolving. The data
will not support it, and it is worth being explicit about why: **89,758 of the 91,486
logged candidates differ from the one before them**, and consecutive candidates differ
by two to six branches — one to three swaps at a time. This is a branch-and-cut node
sequence, not a descent. Played back frame by frame it is a strobe of unrelated
topologies.

So the scrubber *accumulates* instead. At scrub position *t* every branch is shaded by
how often it has been open across candidates 1…*t* — the warm "heat" — with the
current candidate drawn bold white on top. Drag right and the neighbourhood the search
actually lives in paints itself in; drag left and it un-paints. That view survives the
one phase with 17,743 candidates as comfortably as the median one with 22.

### The objective was not in the logs

The callback logs record elapsed time, the open-branch set, the violated contingencies
and the cut count — but no objective. It is recoverable: `manifest.parameters` gives
`ov_lostload_coef = 1`, `ov_n1violation_coef = null`, so the campaign objective *is*
`eval_risk`, which `tools/coordedit/pf.py` already implements and which was validated
against the Julia original to 8.9e-16.

`tools/prep.py` therefore re-evaluates every distinct topology in the corpus (44,223 of
them, ~1 minute on 8 cores) and attaches a true objective to every single candidate.

**Cross-check against ground truth:** for all **372** phases whose `result.json`
recorded an opening set, the recomputed `eval_risk` of that set equals the reported
objective exactly, to the 3-decimal rounding the campaign used.

---

## Things the viewer tells you honestly

These are all real properties of the data, surfaced in the UI rather than smoothed over.

- **The answer is usually not the last frame.** The phase's reported solution is the
  final logged candidate in only **82 of 480** phases. It appears *somewhere* in the
  trace in 372, at a median of 74 % of the way through — the solver keeps exploring
  after it has already seen the winner. The caption under the objective plot says which
  case you are in, every time.
- **108 phases recorded `openings: []`.** For those, `result.json` reports an objective
  but no topology, so it cannot be placed among the candidates. The viewer omits the
  reference line and says so instead of drawing something plausible. This is a logging
  gap in the campaign, not a modelling difference — every phase that *did* record a set
  reconciles exactly.
- **Wall time is confounded** and carries a ⚠ wherever it is shown. Configurations ran
  in a fixed order within each block, so timing picks up machine-state drift. LP-solve
  and Benders-iteration counts are the deterministic quantities; sort by those.
- **Objectives are per unit on a 100 MVA base.** The viewer leads with MW and prints
  the p.u. value underneath, never a bare number.
- **TLF differs by system** — 1.5 for `ieee118`, 1.0 for `ieee57` — and is applied per
  run from the manifest.
- **The outer search is heuristic**, so two configurations ending at different local
  optima is a finding, not a bug. That is what the ghost overlay is for.

## Controls

| | |
|---|---|
| scrubber, `←` `→`, `shift`+arrow | move one / twenty-five candidates |
| `space`, **play** | run the phase, 0.4× to 400× — **1× is 2.5 candidates per second** |
| **ghost a run…** then pick one | overlay a second run's whole residency in violet |
| **heat** / **violations** / **bus labels** | layers on and off |
| drag, wheel — or one finger / two-finger pinch | pan, zoom |
| ☰ and ⓘ (phone only) | the run list and the run details, as sheets |
| **open in coordedit →** | send this exact candidate to the local tool for real flows |

`&play=1` starts a link playing, and `&speed=5` picks the rate it plays at.

The URL is a deep link and updates as you go —
`?run=118_H4_d2_hop1_NO-EMBED&p=0&f=end&ghost=118_H4_d2_hop1_REF` is a specific
candidate of a specific phase with a specific ghost, ready to paste into an email.

### What it deliberately does not draw

No flows, no overloads, no N−1 panel. Those need the DC solver, which is Python, which
means a server — and this had to be something you open from a link. `tools/coordedit`
in the main repo already does all of that and was validated to 1e-14, so rather than
build a second, subtly different copy, every frame carries a link into it.

## Data provenance

Compiled by `tools/prep.py` from the raw campaign output (`tmp/IJEPES/ablation/`, not
redistributed here: 3.6 GB of it is serialized Julia solver state). Each run's
`manifest.json` — `git_sha`, Julia and Gurobi versions, CPU, seed, tolerances, and every
algorithm knob — is carried into `docs/data/runs/<run>.json`.

Two parsing hazards worth recording, both handled:

- The callback logs print Julia `Set` literals via `show`, whose iteration order is not
  stable. Everything is compared as a set of branch indices, never as a string.
- **44 phases have two callback logs in the same directory** — those runs were executed
  twice and both campaigns' logs coexist. `result.json` describes the later one, so the
  parser pins the log by `manifest.started_at` rather than taking whichever sorts first.

## Rebuilding the data

Needs the raw campaign directory and the `coordedit` venv from the main repo:

```bash
python tools/prep.py docs/data
```

---

Objective values, network dumps and the DC model all come from
[`TNROpt`](https://github.com/BenoitJeanson) and its `coordedit` bench.
