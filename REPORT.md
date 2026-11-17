# Phase 1 Report

Rakan Al Qawasmi (20210161), Sameer Amjad Abboud (20220778)  
*Optimal Resource Allocation and Periodic Scheduling*, paper by van Os, Geilen, Hendriks and Basten, in RTAS 2026. DOI [10.1109/RTAS68450.2026.00015](https://doi.org/10.1109/RTAS68450.2026.00015)

## Objective
We implemented the four models presented in the paper, using CPLEX and following the paper: Bounded Monolithic (Section VI, Eqs. 12-18), Bounded Benders (Section VII, Eqs. 19-37 and Algorithm 1), and Bounded Quinton Monolithic and Bounded Quinton Benders, which are Quinton's models with −B ≤ K ≤ B (Section VIII-A).

## Implementation
The code is adapted from the authors' public APSP-toolset (version 1.0.2, MIT License) and trimmed to the four models. We used the full edition of CPLEX 22.2, docplex 2.32 and Python 3.10. Running `python solve_instance.py --model all --instance example/graph__10-1000-0.10-02_0.xml` computes the solution of the benchmark example for all four models. `EQUATIONS.md` links each equation to the code that implements it.

## Changes to the models
| Change | Reason |
|---|---|
| Core point now averages over iterations (both Benders models) | The paper takes the average of the previous and current points. The original code simply used the current point. |
| Stopping tolerance 10⁻⁶ in Bounded Benders | This is the tolerance the paper gives. The original code was using a much tighter 10⁻⁹. |
| Stopping tolerance kept at 10⁻⁹ in Bounded Quinton Benders | This model works with very small numbers (around 0.0001). With 10⁻⁶ the tolerance would be too loose, and the period could be off by approximately 42. |
| Constraint (17)/(30) applies to every task | The original code ignored this constraint for any task with zero time on some resource, even when the task was performed on a slower one. This enabled invalid schedules, for example a period of 3 for a task that takes 5. |
| Duplicate dependencies keep the lowest depth | According to the paper, the lowest depth is the one that counts. The original code kept the last one in the file. |
| `paper_faithful` switch, set to the paper's behaviour by default | The authors' original behaviour remains reproducible for the following phase through `--legacy`. |

We also corrected two crashes and two bugs in the authors' solution checker, added a command-line run script and deleted the plotting code. The details are in `README.md`.

## Tests
We checked every equation, Algorithm 1 and the core-point formula against the code, and a separate review checked them again. For reference, the authors' original code returns µ = 6494 for all four models on our CPLEX. Each change is accompanied by a test that does not pass on the original code but does pass after the change. We also created five small instances whose optimal period we proved by hand (5, 9, 7, 7 and 32). These check that the models agree, that all schedules are valid and that start times remain within Theorem 1's bound. All four models achieve µ = 6494 and pass verification on the benchmark example. Overall, 50 tests pass and 4 are skipped because of the two Bounded Benders restrictions below.

## Findings
- There are some errors in the paper: a wrong sign in Eq. (33), a wrong index in (28) and a misprint in the starting case of the core point. The code already had the correct sign for (33).
- We found two bugs in the authors' solution checker. One marked valid schedules as invalid and the other crashed. We fixed both.
- Bounded Benders has two restrictions on tiny instances. Its master problem can be infeasible, because the authors use the smallest big-M that the paper allows. The master can also suggest an order that the subproblem is unable to schedule, and Algorithm 1 has no way to recover from this. The first restriction cannot occur on the paper's benchmark, since the safety margin is at least 11,928 across all 120 graphs. The second did not occur on the benchmark example.
- The number of Benders iterations depends on the run and changes unless `PYTHONHASHSEED=0` is set. The optimal period remains unchanged.
- As the paper points out, the Quinton models are numerically delicate, since they are based on τ = 1/µ ≈ 0.0001.
- On the benchmark example with seed 0, the corrected core point reduced the number of iterations from 55 to 29 for Bounded Benders and from 92 to 55 for Bounded Quinton Benders. This is just one case, and we do not assume it is a general result.

## Conclusion
All four models are solved with CPLEX and obtain the verified optimum on the benchmark example. Each change from the authors' code is listed above with its reason and has a test.
