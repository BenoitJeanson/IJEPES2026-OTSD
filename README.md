# OTSD campaign replay

**OTSD — Optimal Transmission Switching with De-energization.** Choose transmission
lines to open so the remaining network carries the least load at risk while staying
N−1 secure. Paper:
[preprint on SSRN](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7352475).

An interactive replay of the component-ablation campaign behind that paper: 162 runs, 480 outer local-search phases, and **88,149
individual Benders candidates**, every one of them a topology you can look at.

**→ [Open the viewer](https://benoitjeanson.github.io/IJEPES2026-OTSD/)** ·
[concepts & glossary](https://benoitjeanson.github.io/IJEPES2026-OTSD/?help=1)

Every term below is also in the app itself, behind the **? what am I looking at** button
at the top of the right panel (the **?** in the bar on a phone): what OTSD is, what a
candidate, phase, SBS, `H`, `d_viol` and `d_sol` are, what each `NO-…` configuration
switches off, and the two caveats worth reading before drawing conclusions.

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

### The SBS is drawn, and it had to be dug out

The switchable branch set is the surface the master is allowed to act on; everything
outside it is fixed for that phase. So it is drawn as the surface: branches in the SBS
sit forward, branches outside it recede to near-invisible, and the ones that phase
*added* carry a thicker teal underlay. In phase 1 the whole SBS is "added" — that is
the phase which brings it into being from the heuristic seed. Toggle it with **SBS**.

Membership is per phase, not permanent: in `118_H4_d2_hop1_REF`, 17 branches that sit
outside the SBS in phase 1 are inside it by phase 3, so a branch goes faint → teal →
ordinary as the search widens its surface. Across all 480 phases a median of **63 %**
of the network is outside the SBS.

No text log records SBS *membership* — only `|SBS|` as a number. The tempting shortcut,
taking the union of branches ever opened across a phase's candidates, is badly wrong:
measured over all 480 phases it covers a **median 39 %** of the real set. That would
show the search's footprint and pass it off as its permitted surface.

So `tools/sbs.jl` reads the authoritative set out of the serialized solver state:
`i<N>_state.jls` holds `accumulated_sbs` after phase *N* grew it. Phase 1's set is never
serialized, so it is recomputed by calling TNROpt's own `sa_induced_followed` — the same
call `ablation.jl` makes — rather than reimplementing it.

**Cross-checks, both exact:** every recomputed and every deserialized set matches the
logged `|SBS|` for all 162 runs and all 480 phases; and all 88,149 candidates lie inside
the SBS their phase searched, which they must, since the master cannot open anything else.

> **A labelling trap.** `result.json`'s `sbs_size` for phase *N* is the SBS **after**
> phase *N* grew it — that is, the set phase *N+1* searches. The set phase *N* actually
> searched is the previous one (for phase 1, `init_sbs`). The viewer always shows the
> set the phase searched, and the `+n` beside it is what that phase added.

### The gantt: every branch against every candidate

Under the scrubber, one row per branch the phase ever opened, one column per
candidate, sharing the scrubber's x-axis so the two line up. Rows are ordered by how
often the branch is open, which separates two things at a glance: the solid block at
the top is the stable core the search never gives up, and the scatter below is what it
churns through. Columns where the subproblem found **no violated contingency** — the
feasible proposals — are washed green with a solid foot, and the phase's reported
result is the blue column — and those same feasible candidates are the green dots in
the objective plot, the only ones the incumbent can step to.

Click or drag anywhere on it to seek. Drag its top edge to resize; past about 8 px a
row it grows a gutter with the branch names. **gantt** toggles it — on by default on
desktop, off on a phone, where it also gets capped so it cannot starve the network.

### Comparing two runs

**pin this run as reference** freezes the run you are looking at as a violet halo —
where it spent its candidates, across all its phases — and then lets you pick another
run to view against it. The violet appears the moment you click, so the button's effect
is never in doubt, and the run you pick is drawn on top. The difference between the two
is the divergence between the searches, which for a heuristic outer loop is a finding
rather than a bug.

### The green line is the incumbent, not a running minimum

The objective plot draws every candidate as a dot, but the line through them is the
best **feasible** objective so far — feasible meaning the subproblem found no violated
contingency. A running minimum over all candidates would track proposals the subproblem
went on to reject, which bounds nothing. What you get instead is a step function that
moves only when a feasible candidate beats the incumbent, and that lands on the phase
result. Feasible candidates are the larger green dots; they are a median **30 %** of a
phase's candidates, and every phase has at least one.

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
| **pin this run as reference**, then pick one from the list | the run you were viewing turns violet and stays as the reference; the run you pick comes to the front, drawn over it |
| **heat** / **violations** / **SBS** / **bus labels** | layers on and off |
| drag, wheel — or one finger / two-finger pinch | pan, zoom |
| click or drag the gantt | seek to that candidate |
| drag the gantt's top edge, or a panel's inner edge | resize (remembered per browser) |
| ☰ and ⓘ (phone only) | the run list and the run details, as sheets |

`&play=1` starts a link playing, `&speed=5` picks the rate, and `&gantt=1` forces the
matrix on.

The URL is a deep link and updates as you go —
`?run=118_H4_d2_hop1_NO-EMBED&p=0&f=end&ghost=118_H4_d2_hop1_REF` is a specific
candidate of a specific phase with a specific ghost, ready to paste into an email.

### What it deliberately does not draw

No flows, no overloads, no N−1 panel. Those need the DC solver, which is Python, which
means a server — and this had to be something you open from a link. `tools/coordedit`
in the main repo already does all of that and was validated to 1e-14, so rather than
build a second, subtly different copy, this one stays out of that business.

## Data provenance

Compiled by `tools/prep.py` from the raw campaign output (`tmp/IJEPES/ablation/`, not
redistributed here: 3.6 GB of it is serialized Julia solver state). Each run's
`manifest.json` — `git_sha`, Julia and Gurobi versions, CPU, seed, tolerances, and every
algorithm knob — is carried into `docs/data/runs/<run>.json`.

Two parsing hazards worth recording, both handled:

- The callback logs print Julia `Set` literals via `show`, whose iteration order is not
  stable. Everything is compared as a set of branch indices, never as a string.
- **44 phases have two callback logs in the same directory** — those runs were executed
  twice and both campaigns' logs coexist. `result.json` describes the later one, so both
  tools pin the log by the stamp *nearest* `manifest.started_at`. Nearest, not equal: the
  filename stamp is read from the clock just before the manifest is written, and in
  `118_H4_d2_hop1_REF_s3` the two differ by one second.

## Rebuilding the data

Needs the raw campaign directory and the `coordedit` venv from the main repo:

```bash
julia --project=/path/to/tnr tools/sbs.jl     # SBS membership out of the .jls
python tools/prep.py docs/data                # everything else, merging the above
```

---

Objective values, network dumps and the DC model all come from
[`TNROpt`](https://github.com/BenoitJeanson) and its `coordedit` bench.
