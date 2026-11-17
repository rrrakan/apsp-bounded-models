# EQUATIONS.md — paper equations mapped to the code

Paper: van Os, Geilen, Hendriks, Basten, *Optimal Resource Allocation and Periodic
Scheduling*, RTAS 2026 (DOI 10.1109/RTAS68450.2026.00015). All paths are relative to
the repository root; line numbers refer to the delivered files.

## Shared definitions

| Item | Paper | File:line |
|---|---|---|
| Theorem 1 bound B = \|T\|+\|R\| | Sec. V | set in each model's `setup_variables`: `bounded_monolithic.py:16`, `bounded_master.py:16`, `bounded_primal_sub.py:14`, `bounded_dual_sub.py:16`, `bounded_auxilary_sub.py:16`, `bounded_quinton_monolithic.py:16`, `bounded_quinton_master.py:17` |
| I = {0, …, B} | Sec. VI | `range(self.data.BOUND + 1)` wherever repetitions are enumerated (e.g. `bounded_monolithic.py:22,62`) |
| M ≥ max e · \|T\| · B (the paper's stated minimum is used) | Sec. VI (p. 8) | `bounded_monolithic.py:38`, `bounded_master.py:39`, `bounded_primal_sub.py:34` |
| µ̂ = Σ_t max_r e(t,r) | Sec. VII-B (p. 9) | `bounded_master.py:37` (`mu_ub`) |
| Lowest-depth duplicate dependency dominates | Sec. V (Thm 1 proof) | `sdf3_python_utilities/sdf3_data_conversion.py:20-24` (correction 4) |

## Bounded Monolithic (12)–(18) — `models/bounded/bounded_monolithic/bounded_monolithic.py`

| Eq. | Content | Line |
|---|---|---|
| (12) | min µ | 73 |
| (13) | Σ_r x_{t,r} = 1 ∀t | 47 |
| (14) | dependencies: s_{t2} ≥ s_{t1} − nµ + Σ x·e | 41–44 |
| (15) | no-overlap left, M(3 − x − x − y) | 63–66 |
| (16) | no-overlap right, M(2 − x − x + y) | 67–70 |
| (17) | s_t + Σ x·e ≤ s_t + µ, **∀t ∈ T** (correction 3) | 49–53 |
| (18) | s_t ≤ Bµ | 55 |

Tolerance 10⁻⁶ via solver settings (Sec. VIII-B): line 13.

## Benders master MP (19)–(25) — `models/bounded/bounded_Benders/bounded_master.py`

| Eq. | Content | Line |
|---|---|---|
| (19) | min z (variable named `mu` in code) | 22, 68 |
| (20) | Σ_r x_{t,r} = 1 | 42 |
| (21) | z ≥ c ∀c ∈ C (Benders cuts) | added by the driver: `bounded_magnanti_benders.py:101-144` (`__apply_cuts`) |
| (22) | predefined dependency cut with w, µ̂ | 46–48 |
| (23) | predefined no-overlap left with w, µ̂ | 54–57 |
| (24) | predefined no-overlap right with w, µ̂ | 58–61 |
| (25) | w_t + Σ x·e ≤ w_t + µ̂, ∀t ∈ T | 64–66 |

## Primal subproblem SP (26)–(31) — `models/bounded/bounded_Benders/bounded_primal_sub.py`

| Eq. | Content | Line |
|---|---|---|
| (26) | min µ | 59 |
| (27) | dependencies with x̄ fixed | 37–39 |
| (28) | no-overlap left with x̄, ȳ fixed | 50–53 |
| (29) | no-overlap right with x̄, ȳ fixed | 54–57 |
| (30) | s_t + Σ x̄·e ≤ s_t + µ, **∀t ∈ T** (correction 3) | 41–44 |
| (31) | s_t ≤ Bµ | 45 |

Subproblem tolerance 10⁻⁹ (Sec. VIII-B): line 12.

## Dual subproblem DSP (32)–(35) — `models/bounded/bounded_Benders/bounded_dual_sub.py`

| Eq. | Content | Line |
|---|---|---|
| π variables (π^p, π^nl, π^nr, π^nd, π^s ≥ 0; π^nd over **all** tasks, correction 3) | | 19–35 |
| (33) | µ-row: Σnπ^p − Σkπ^nl **+ Σkπ^nr** + Σπ^nd + BΣπ^s ≤ 1 (the printed (33) has a sign typo; the code is the correct dual — see README) | 55–70 |
| (34) | s_t-row ≤ 0 | 72–96 |
| (32)+(35) | max ob(π, x̄, ȳ) | 98–112 |

## Magnanti–Wong ASP (36)–(37) — `models/bounded/bounded_Benders/bounded_auxilary_sub.py`

| Eq. | Content | Line |
|---|---|---|
| π variables (π^nd over all tasks) | | 19–35 |
| (33) repeated | | 59–74 |
| (34) repeated | | 76–100 |
| (37) | ob(π, x̄, ȳ) = ob̄ (equality) | 102–122 |
| (36) | max ob(π, x̄′, ȳ′) (core point in objective) | 124–140 |

## Cut, core point, Algorithm 1

| Item | Paper | File:line |
|---|---|---|
| Cut z ≥ ob(π̄, x, y) via (35) with open x, y (π^nd over all tasks, correction 3) | (21)+(35) | `bounded_magnanti_benders.py:101-144` |
| Papadakos core point v̄′₀ = ½v̄₀; v̄′ₖ = ½v̄′ₖ₋₁ + ½v̄ₖ (correction 1) | Sec. VII-D (p. 10) | `bounded_magnanti_benders.py:44-76` (`_calculate_core_point_paper`; the authors' version is preserved verbatim as `_calculate_core_point_legacy`, line 78, for `paper_faithful=False`); Quinton analogue `bounded_quinton_magnanti_benders.py:41-73` (legacy path 75-95) |
| Algorithm 1 loop (MP → DSP → core point → ASP → cut) | Sec. VII-D | `combinatorial_optimization_tools/magnanti_wong_papadakos_benders_decomposition.py:104-254` (`__solve`) |
| Stop when UB − LB ≤ ε (engine uses strict <; documented note) | Alg. 1 | engine line 261 |
| ε = 10⁻⁶ in Bounded Benders (correction 2); BQB keeps 10⁻⁹ (recorded decision) | Sec. VIII-B | `bounded_magnanti_benders.py:28-29`; engine default line 11 |

## Quinton models (Section VIII-A, models 5–6)

The Quinton formulations implement eqs. (6)–(11) of Quinton et al. (2020), as the
paper describes; the in-code comments carry that numbering
(`bounded_quinton_monolithic.py`, `bounded_Quinton_Benders/*.py`).

| Item | Paper | File:line |
|---|---|---|
| Bound −B ≤ K_{t1,t2} ≤ B (the paper's addition to Quinton's model) | Sec. VIII-A (p. 11) | `bounded_quinton_monolithic.py:21-22`, `bounded_quinton_master.py:20-21` (`lb=-self.data.BOUND, ub=self.data.BOUND`) |
| Objective max τ, period µ = 1/τ | Sec. VIII-A | `bounded_quinton_monolithic.py:73`, `bounded_quinton_primal_sub.py:59`; `solve_instance.py` reports `1/get_objective()` |
