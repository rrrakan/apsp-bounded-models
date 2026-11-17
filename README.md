# APSP bounded models 

**Rakan Al Qawasmi (20210161), Sameer Amjad Abboud (20220778)**

Project report: [REPORT.md](REPORT.md) (PDF: [report.pdf](report.pdf)).

Implementation of the four bounded models of van Os, Geilen, Hendriks, Basten,
*Optimal Resource Allocation and Periodic Scheduling*, RTAS 2026
(DOI 10.1109/RTAS68450.2026.00015), derived from the authors' APSP-toolset
release 1.0.2 (https://github.com/TUE-EE-ES/APSP-toolset) with the corrections
listed below. The paper is the authority; every deviation from release 1.0.2 is
listed in this file with its reason. `EQUATIONS.md` maps every equation to file
and line.

Models:
1. **BM** — Bounded Monolithic, eqs. (12)–(18)
2. **BB** — Bounded Benders (Magnanti–Wong–Papadakos), eqs. (19)–(37), Algorithm 1
3. **BQM** — Bounded Quinton Monolithic (Sec. VIII-A, model 5)
4. **BQB** — Bounded Quinton Benders (Sec. VIII-A, model 6)

## Setup (this environment — already done; recorded for reproduction)

- Windows, PowerShell. Python **3.10.11** in `C:\apsp_project\.venv`
  (`.venv\Scripts\Activate.ps1`; every command runs inside it).
- **cplex 22.2.0.0** (pip) + **docplex 2.32.264**, plus pandas, interruptingcow,
  pytest. CPLEX Studio 22.2 is installed at
  `C:\Program Files\IBM\ILOG\CPLEX_Studio222`.
- From CPLEX 22.1.2 on, the Python API no longer ships with CPLEX Studio. The pip
  `cplex` package was installed and then **unlocked to the full edition** with:
  `docplex config --upgrade "C:\Program Files\IBM\ILOG\CPLEX_Studio222"`.
- **Never `pip install` / `--upgrade` / reinstall `cplex`.** That replaces the
  unlocked full edition with the size-limited Community Edition. If it happens,
  rerun the `docplex config --upgrade` command above.
- Version note: the paper used CPLEX 22.1.1 and docplex 2.29.241; this
  environment uses CPLEX 22.2 and docplex 2.32.264. Results on the example
  instance match the paper's µ; iteration counts and times differ with solver
  version (and with hash seed, see below).

## How to run

From the repository root, inside the venv:

```
python solve_instance.py --model all --instance example/graph__10-1000-0.10-02_0.xml --time-limit 3600
```

- `--model {BM,BB,BQM,BQB,all}` (default `all`)
- `--instance PATH` (default: the example graph; expected µ = 6494 for all four)
- `--time-limit N` seconds per model (default 3600; must be positive — the
  Benders models always need a time limit)
- `--legacy` sets `paper_faithful=False`: the authors' release-1.0.2 core-point
  logic and ε = 10⁻⁹ in Bounded Benders (corrections 3–4 stay on in both modes).
  Default is the paper-faithful mode.
- The Quinton models maximize τ = 1/µ; the script reports the period as
  `1/get_objective()`.
- **Repeatable Benders iteration counts need `PYTHONHASHSEED=0`**
  (PowerShell: `$env:PYTHONHASHSEED='0'`). Resource sets iterate in hash order,
  which reorders CPLEX rows and changes the solver's path run-to-run.
  Objectives do not depend on it; the script prints a note when the seed is not
  pinned.

Tests: `python -m pytest tests` from the repository root (log of the delivered
run: `evidence/final_pytest_log.txt`; the delivered example run is
`evidence/final_example_log.txt`).

## Model-to-section map

| Part | Paper | Code |
|---|---|---|
| Dependencies (1), no-overlap (2), lowest-depth rule | Sec. III–V | `sdf3_python_utilities/`, `verify_instance.py` |
| Theorem 1 bound B = \|T\|+\|R\| | Sec. V | each model's `setup_variables` (see EQUATIONS.md) |
| Bounded Monolithic (12)–(18) | Sec. VI | `models/bounded/bounded_monolithic/` |
| Master MP (19)–(25) | Sec. VII-B | `models/bounded/bounded_Benders/bounded_master.py` |
| SP (26)–(31) / DSP (32)–(35) | Sec. VII-C | `bounded_primal_sub.py` / `bounded_dual_sub.py` |
| ASP (36)–(37), core point, cut | Sec. VII-D | `bounded_auxilary_sub.py`, `bounded_magnanti_benders.py` |
| Algorithm 1 | Sec. VII-D | `combinatorial_optimization_tools/magnanti_wong_papadakos_benders_decomposition.py` |
| Quinton models, −B ≤ K ≤ B | Sec. VIII-A | `models/bounded/bounded_Quinton_*/` |
| Tolerances (10⁻⁶ optimality; 10⁻⁹ subproblems) | Sec. VIII-B | model setup / engine / driver |

## Changes to the models (each with its reason)

| # | Change | Where | Reason |
|---|---|---|---|
| C1 | Core point now averages over iterations (both Benders models): v̄′₀ = ½v̄₀; v̄′ₖ = ½v̄′ₖ₋₁ + ½v̄ₖ (paper Sec. VII-D) | `bounded_magnanti_benders.py`, `bounded_quinton_magnanti_benders.py` (`_calculate_core_point_paper`) | The paper averages the current point with the previous one; the original (1.0.2) code skipped the average and just used the current point (and its µ/z component has a `0.5 +` typo). The authors' logic is preserved verbatim (`_calculate_core_point_legacy`) behind the switch. |
| C2 | Stopping tolerance ε = 10⁻⁶ in Bounded Benders | `bounded_magnanti_benders.py:28-29` | The value the paper states (Sec. VIII-B: the 10⁻⁶ optimality tolerance "is reflected in the value of epsilon in Algorithm 1"); the original code used a stricter 10⁻⁹ (the engine default, never overridden). |
| ε-BQB | Stopping tolerance **kept** at 10⁻⁹ in Bounded Quinton Benders (both modes) | `bounded_quinton_magnanti_benders.py` | This model works with very small numbers (its objective is τ = 1/µ ≈ 0.0001), so the paper's 10⁻⁶ would be too loose: on the example the period could be off by about 42 time units. Recorded project decision. |
| C3 | Constraint (17)/(30) applies to **every** task; π^nd consistently over all tasks in DSP, ASP and the cut | `bounded_monolithic.py`, `bounded_primal_sub.py`, `bounded_dual_sub.py`, `bounded_auxilary_sub.py`, `bounded_magnanti_benders.py` | The original code skipped it for any task that takes zero time on SOME resource, even when the task runs on a resource where it takes longer. That allowed invalid schedules (e.g. a period of 3 for a task that takes 5). Applying it to every task is always correct and matches the paper (∀t ∈ T); the master's (25) never had the skip. Active in both modes (no effect on the benchmark, which has no zero execution times). |
| C4 | Duplicate dependencies keep the **lowest** depth | `sdf3_data_conversion.py:20-24` | The paper (Thm 1 proof) says the lowest depth counts; the original code kept whichever duplicate came last in the file. Active in both modes. |
| S | `paper_faithful` switch (default = the paper's behaviour); `--legacy` CLI flag restores the 1.0.2 core point and ε | both Benders drivers, `solve_instance.py` | So the authors' original results (Figure 7 was produced with the 1.0.2 behaviour) can still be reproduced in the next phase. `--legacy` reproduces the seed-0 baseline exactly (µ and iteration counts). |
| (33) | **No code change.** | `bounded_dual_sub.py`, `bounded_auxilary_sub.py` | The paper prints a wrong sign in the dual constraint (33); the original code already had the correct one (−Σkπ^nl **+** Σkπ^nr). See "Paper typos" below. |

## Fixes outside the models (no model logic changed)

| # | Change | Where | Reason |
|---|---|---|---|
| R1 | Crash fix: `solve(time_out=None)` means "no limit" | engine, 4 guarded checks | 1.0.2 crashes with `TypeError: None - float`, then crashes again inside its own error handler. |
| R2 | Crash fix: `write_error_log` guarded when `best_master_problem_extract is None` | engine `write_error_log` | 1.0.2's error handler crashes (`AttributeError`) when the error occurs before the first completed iteration. |
| V1 | Verifier operator-precedence fix: `(y or z) and value==1` | `verify_instance.py`, both parse branches | In 1.0.2 the value check binds only to `z`, so an exported epsilon binary (observed `y = 4.17e-11`) is parsed as an ordering decision → false violation on a valid schedule. |
| V2 | Verifier default fills consolidated before `verify()` | `verify_instance.py` | In 1.0.2 the order/start defaults only run if some y/z/K (or s/u) variable was exported; with none (no shared resource), verification crashes with `KeyError`. |
| CLI | Run script `solve_instance.py` with `--model/--instance/--time-limit/--legacy` | `solve_instance.py` | Single-instance runs with machine-readable `RESULT` lines; keeps the authors' solve / write_solution / verify flow. |
| T | Trim: plotting removed entirely (figure methods, `Show`, matplotlib/fontTools imports, the engine's abstract `build_figure`), experiment scripts and benchmark data removed (one example graph kept), unbounded Quinton baselines removed; `models/` stays a namespace package | throughout | Only what this phase needs, and no plotting dependency. No formulation change — verified by identical µ and identical seed-0 iteration counts before/after the trim. |

Documented but **not** changed (from 1.0.2): implicit imports through
`import *`; the Quinton master's `z ≤ 1e20 − 1e10` safety bound (replacing it
was tested and does not fix the numerical fragility).

## Paper typos found (code keeps the mathematically correct form)

- **(33) sign:** the printed µ-row of the DSP reads `− Σ (kπ^nl + kπ^nr)`
  (one Σ, the minus covering both terms). Deriving the dual of SP (26)–(31)
  gives −Σkπ^nl **+** Σkπ^nr (µ enters (28) and (29) with opposite senses), and
  the printed (35) is consistent only with the + reading. The code (1.0.2 and
  ours) implements the correct dual.
- **(28) index:** the paper quantifies (28) with "∀m ∈ I" while the constraint
  body uses k. The code uses k throughout.

## Known issues and behaviour

- **Bounded Benders converges slowly by design.** The big-M cuts are weak; the
  LB can stay at 0 for many iterations and the master MILP takes 95–98% of the
  run time. The authors report 52–1421 iterations on hard 10-task instances.
- **The Quinton models are numerically fragile** (τ ≈ 10⁻⁴ with master
  coefficients spanning ~10⁷). Solvers can return a false "optimal" master or
  a suboptimal µ; the paper reports this itself. CPLEX reporting BQM as
  "integer optimal, tolerance" is normal (optimal within CPLEX's default gap).
- **Hash seed:** Benders iteration counts vary run-to-run unless
  `PYTHONHASHSEED=0` (set iteration order → CPLEX row order → solver path).
  With seed 0 on the example: default mode BB 29 / BQB 55 iterations; `--legacy`
  BB 55 / BQB 92 (identical to the authors' 1.0.2 under the same seed). The
  default-vs-legacy difference is an observation on this one instance with
  seed 0, not a general convergence result.
- **MIP-gap note:** the paper's 10⁻⁶ "through solver settings" is implemented
  (as in 1.0.2) as the LP simplex optimality tolerance; the CPLEX MIP gap stays
  at its default 10⁻⁴. The paper does not name the parameter. Documented only.
- **Strict-< note:** the engine stops when UB − LB **<** ε (strict), the paper
  writes ≤. Practically irrelevant; documented only.
- **Bound after a timeout:** after a timeout, the Benders lower bound comes from
  the master's best solution so far, not a proven bound (inherited from 1.0.2;
  documented only).
- **Finding B1 — Bounded Benders master can be infeasible on tiny instances.**
  The authors chose the **minimum** M the paper allows ("at least
  max e·\|T\|·B"). In the master's predefined cuts (22)–(25) the pair cut at
  k = B leaves only M − e − B·µ̂ of big-M slack, and when execution times are
  near-uniform this slack can be smaller than a dependency-forced w-spread,
  making the master integer infeasible (observed on the hand-made paper
  Example 1: slack 3 vs. forced spread 6; conflict-refiner proof in
  `tests/test_ground_truth.py`). **This cannot occur on the authors' 120
  benchmark graphs:** there, M − B·µ̂ is 11 928 / 38 250 / 101 898
  (min/median/max) while the largest depth-0 dependency-chain spread is
  3 162 / 6 550 / 12 357, and no instance comes within a margin of 3 237 of
  the condition. Property of the instances, computed without solving.
- **Finding B2 — Bounded Benders can abort on a subproblem-infeasible master
  proposal.** The master (22)–(25) has no bound on w analogous to (18)/(31),
  so it may propose an ordering ȳ (e.g. y at k = B) that no schedule with
  s ≤ Bµ can realize; Algorithm 1 has no feasibility cuts, so the run aborts
  with "dual sub problem unbounded" and reports the incumbent (observed µ = 35
  vs. optimum 32 on a hand-made tiny instance; release 1.0.2 behaves
  identically). Observed on hand-made tiny instances only; never observed on
  the benchmark. Both findings are pinned by tests in
  `tests/test_ground_truth.py` and were left unfixed deliberately: the code
  matches the paper, and changing M or adding w-bounds/feasibility cuts would
  change the model beyond what the paper calls for.
