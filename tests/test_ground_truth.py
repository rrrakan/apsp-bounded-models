"""Stage 4 ground-truth and property tests.

Four tiny hand-made SDF3 instances, each with the optimal period mu* proven by
hand below. For every instance, all four models are solved, every solution is
exported and checked with the authors' verifier, all models must agree on mu
(Quinton models: mu = 1/tau), and for BM/BB every start time must satisfy the
Theorem 1 bound s <= B*mu, eqs. (18)/(31). The Quinton models bound K, not the
start times, so no start-time bound is asserted for them.

HAND PROOFS (each optimum is integer by construction; the verifier rounds
start times to integers):

example1.xml (paper Example 1, allocation forced: x,z only on r1=proc_0, y only
on r2=proc_1; e(x)=3, e(y)=3, e(z)=2; deps x->y depth 0, y->x depth 3, y->z
depth 0). LOWER BOUND: proc_0 must execute x and z once per period without
overlap, so by (7): mu >= e(x)+e(z) = 3+2 = 5. ATTAINED by the paper's schedule
sigma(x)=0, sigma(y)=7, sigma(z)=13, mu=5: x->y: 7 >= 0+3; y->x (depth 3):
0 >= 7+3-3*5 = -5; y->z: 13 >= 7+3; proc_0 occupancy mod 5: x in [0,3),
z in [3,5) - no overlap. Hence mu* = 5.

single_resource.xml (a0,a1,a2 only on proc_0, e = 2,3,4, no dependencies).
LOWER BOUND: one resource executes all three tasks once per period without
overlap, so by (7): mu >= 2+3+4 = 9. ATTAINED by sigma = 0,2,5 (back to back).
Hence mu* = sum(e) = 9.

cycle_bound.xml (a0,a1,a2 each alone on its own resource, e = 4,4,6; cycle
a0->a1 (depth 0), a1->a2 (depth 0), a2->a0 (depth 2)). LOWER BOUND from the
cycle, chaining constraint (1): sigma(a1) >= sigma(a0)+4, sigma(a2) >=
sigma(a1)+4, sigma(a0) >= sigma(a2)+6-2*mu => 0 >= 14 - 2*mu => mu >= 7.
Resource loads are only 4,4,6 < 7, so the CYCLE determines the optimum.
ATTAINED by sigma = 0,4,8 with mu=7: a2->a0: 0 >= 8+6-14 = 0; no resource is
shared. Hence mu* = (4+4+6)/2 = 7.

balanced_load.xml (a0..a3 independent, each runnable on proc_0 or proc_1 with
e = 3,3,4,4). LOWER BOUND: every task executes once per period on its allocated
resource; total load 3+3+4+4 = 14 over 2 resources, so some resource carries
>= 14/2 = 7, and by (7) mu >= 7. ATTAINED by the balanced split {a0,a2} on
proc_0 and {a1,a3} on proc_1 (loads 3+4 = 7 each), sigma = 0,0,3,3.
Hence mu* = 7.

example1_scaled.xml (Example 1's structure with e(x) = 30 instead of 3).
LOWER BOUND: proc_0 runs x and z once per period, so mu >= 30+2 = 32; the
dependency cycle x->y->x gives only mu >= (30+3)/3 = 11, so the load binds.
ATTAINED by sigma(x)=0, sigma(y)=30, sigma(z)=62, mu=32: x->y: 30 >= 0+30;
y->x (depth 3): 0 >= 30+3-3*32 = -63; y->z: 62 >= 30+3; proc_0 occupancy
mod 32: x in [0,30), z in [62,64) -> [30,32) - no overlap. Hence mu* = 32.

KNOWN PAPER LIMITATIONS of Bounded Benders (documented and pinned, not
"fixed"; the authors' original 1.0.2 behaves identically, and the code matches
the paper's eqs. (19)-(25) exactly per AUDIT.md):

1. example1.xml: the MASTER is integer infeasible, although BM, BQM and BQB
   all find mu* = 5. Proof from the paper's own constants: predefined cut (22)
   chains w_a1 >= w_a0 + 3 and w_a2 >= w_a1 + 3, so w_a2 >= w_a0 + 6; the
   pair cut (23) for (a2,a0) at k = B = 5 gives, even with y = 0,
   w_a2 + e_a2 + B*mu_hat <= w_a0 + M, i.e. w_a2 <= w_a0 + (45-2-40) =
   w_a0 + 3. Contradiction. The paper's "sufficiently large" M (at least
   max_e*|T|*B = 45) is too small relative to B*mu_hat = 40 when execution
   times are near-uniform; on the authors' benchmark M - B*mu_hat is huge and
   the cuts never bind this way. Pinned by
   test_bb_example1_master_infeasible (infeasibility is model property, so
   this is solver-path independent).

2. example1_scaled.xml: the master may propose an ordering whose SUBPROBLEM is
   infeasible, and Algorithm 1 has no feasibility-cut path, so the engine
   aborts with status "dual sub problem unbounded" and a suboptimal incumbent
   (observed: mu = 35 after iteration 1, while mu* = 32). Mechanism
   (verified by dumping the trial point): the master sets y_(a0,a2,k) = 1 for
   ALL k = 0..5 - master-feasible because w is unbounded above in (19)-(25) -
   but y at k = B forces s_a0 + e_a0 + B*mu <= s_a2 <= B*mu (eq. (31)), i.e.
   s_a0 <= -e_a0 < 0: SP infeasible, DSP unbounded. The paper's claim that
   the predefined cuts make every master proposal schedulable does not hold
   here. Which master optimum CPLEX returns is solver-path dependent, so the
   pinning test accepts either the documented abort or correct convergence to
   32 (test_bb_example1_scaled_abort_or_optimum).
"""
import pytest

from conftest import DATA_DIR

from models.utilities import APSPModelData
import models.bounded.bounded_monolithic as bounded_monolithic
import models.bounded.bounded_Benders as bounded_Benders
import models.bounded.bounded_Quinton_monolithic as bounded_Quinton_monolithic
import models.bounded.bounded_Quinton_Benders as bounded_Quinton_Benders
from verify_instance import verify_solution

INSTANCES = {
    "example1": 5,
    "single_resource": 9,
    "cycle_bound": 7,
    "balanced_load": 7,
    "example1_scaled": 32,
}
# Known Bounded Benders limitations of the paper's formulation on these
# instances (see module docstring); pinned by their own tests below instead of
# solved in the fixture.
BB_KNOWN_LIMITATIONS = {"example1", "example1_scaled"}
MODELS = ["BM", "BB", "BQM", "BQB"]
IDENTIFIERS = {
    "BM": "Bounded_Monolithic",
    "BB": "Bounded_Benders",
    "BQM": "Bounded_Quinton_Monolithic",
    "BQB": "Bounded_Quinton_Benders",
}
TIME_LIMIT = 600


@pytest.fixture(scope="module")
def solutions(tmp_path_factory):
    """Solve all instances with all four models once; verify each solution.

    Returns {(instance, model): {"mu": float, "verify": str,
                                 "starts": {task: s} or None, "B": int}}.
    """
    results = {}
    results_root = str(tmp_path_factory.mktemp("results")) + "\\"

    for graph, _ in INSTANCES.items():
        data = APSPModelData(DATA_DIR, graph + ".xml")
        B = data.NR_TASKS + data.NR_RESOURCES
        instance_dir = results_root + graph + "\\"
        import os
        os.makedirs(instance_dir, exist_ok=True)

        # Bounded Monolithic
        model = bounded_monolithic.BoundedMonolithic("BM_" + graph, data)
        model.setup_model("")
        model.solve(False, time_limit=TIME_LIMIT)
        model.write_solution(instance_dir, graph + "_Bounded_Monolithic")
        v = verify_solution("Bounded_Monolithic", DATA_DIR, graph, results_root)
        starts = {t: var.solution_value for t, var in model.decision_vars.s.items()}
        results[(graph, "BM")] = {"mu": model.get_objective(), "verify": v,
                                  "starts": starts, "B": B}

        # Bounded Benders (skipped on the documented paper limitations,
        # pinned by their own tests below)
        if graph not in BB_KNOWN_LIMITATIONS:
            driver = bounded_Benders.BoundedMagnantiWongBenders("BB_" + graph, data)
            driver.set_iterations(100000)
            driver.solve(False, log_verbosity=0, time_out=TIME_LIMIT)
            driver.write_solution(instance_dir, graph + "_Bounded_Benders")
            v = verify_solution("Bounded_Benders", DATA_DIR, graph, results_root)
            # write_solution solved the primal subproblem for the best trial values.
            psub = driver.primal_sub_problem
            starts = {t: var.solution_value for t, var in psub.decision_vars.s.items()}
            results[(graph, "BB")] = {"mu": driver.get_objective(), "verify": v,
                                      "starts": starts, "B": B,
                                      "primal_mu": psub.get_objective()}

        # Bounded Quinton Monolithic (mu = 1/tau; K is bounded, not the starts)
        model = bounded_Quinton_monolithic.BoundedQuintonMonolithic("BQM_" + graph, data)
        model.setup_model("")
        model.solve(False, time_limit=TIME_LIMIT)
        model.write_solution(instance_dir, graph + "_Bounded_Quinton_Monolithic")
        v = verify_solution("Bounded_Quinton_Monolithic", DATA_DIR, graph, results_root)
        results[(graph, "BQM")] = {"mu": 1 / model.get_objective(), "verify": v,
                                   "starts": None, "B": B}

        # Bounded Quinton Benders (mu = 1/tau)
        driver = bounded_Quinton_Benders.BoundedQuintonMagnantiWongBenders("BQB_" + graph, data)
        driver.set_iterations(100000)
        driver.solve(False, log_verbosity=0, time_out=TIME_LIMIT)
        driver.write_solution(instance_dir, graph + "_Bounded_Quinton_Benders")
        v = verify_solution("Bounded_Quinton_Benders", DATA_DIR, graph, results_root)
        results[(graph, "BQB")] = {"mu": 1 / driver.get_objective(), "verify": v,
                                   "starts": None, "B": B}

    return results


def models_for(graph):
    return [m for m in MODELS if not (m == "BB" and graph in BB_KNOWN_LIMITATIONS)]


@pytest.mark.parametrize("graph", list(INSTANCES))
@pytest.mark.parametrize("model", MODELS)
def test_verifier_passes(solutions, graph, model):
    """Requirement 1: the verifier passes on every solution in the suite."""
    if model not in models_for(graph):
        pytest.skip("documented BB paper limitation on this instance "
                    "(pinned by its own test below)")
    assert solutions[(graph, model)]["verify"] == "verification successful", \
        f"{model} on {graph}: {solutions[(graph, model)]['verify']}"


@pytest.mark.parametrize("graph", list(INSTANCES))
def test_models_agree(solutions, graph):
    """Requirement 2: all models agree on mu (tolerance 1e-6)."""
    mus = {m: solutions[(graph, m)]["mu"] for m in models_for(graph)}
    ref = mus["BM"]
    for m, mu in mus.items():
        assert mu == pytest.approx(ref, rel=1e-6), f"{graph}: {mus}"


@pytest.mark.parametrize("graph", list(INSTANCES))
def test_ground_truth_optimum(solutions, graph):
    """Requirement 5: mu equals the hand-proven optimum (proofs in the module
    docstring)."""
    expected = INSTANCES[graph]
    for m in models_for(graph):
        assert solutions[(graph, m)]["mu"] == pytest.approx(expected, rel=1e-6), \
            f"{m} on {graph}: mu = {solutions[(graph, m)]['mu']}, proven optimum {expected}"


@pytest.mark.parametrize("graph", list(INSTANCES))
@pytest.mark.parametrize("model", ["BM", "BB"])
def test_theorem1_start_time_bound(solutions, graph, model):
    """Requirement 3: eqs. (18)/(31): s_t <= B*mu for BM and BB.

    For BB the start times come from the primal subproblem solved for the best
    trial values, whose own objective equals the reported mu up to epsilon.
    """
    if model not in models_for(graph):
        pytest.skip("documented BB paper limitation on this instance")
    r = solutions[(graph, model)]
    mu = r["mu"]
    for t, s in r["starts"].items():
        assert s <= r["B"] * mu + 1e-6, \
            f"{model} on {graph}: s[{t}] = {s} > B*mu = {r['B']}*{mu}"


def test_bb_example1_master_infeasible():
    """Pins the documented paper limitation (module docstring): on Example 1
    the Bounded Benders master with the paper's predefined cuts and constants
    (M = max_e*|T|*B = 45, mu_hat = 8, B = 5) is integer infeasible, while
    BM, BQM and BQB all find mu* = 5 (test_ground_truth_optimum). If a future
    change makes BB solve this instance, this test fails so the change gets
    reviewed against the paper.
    """
    data = APSPModelData(DATA_DIR, "example1.xml")
    driver = bounded_Benders.BoundedMagnantiWongBenders("BB_ex1_infeasible", data)
    driver.set_iterations(100)
    driver.solve(False, log_verbosity=0, time_out=60)
    assert driver.get_status() == "master problem integer infeasible", \
        f"expected the documented master infeasibility, got: {driver.get_status()}"


def test_bb_example1_scaled_abort_or_optimum():
    """Pins documented paper limitation 2 (module docstring): on
    example1_scaled the master may propose an ordering (y_(a0,a2,k) = 1 at
    k = B) whose subproblem is infeasible; Algorithm 1 has no feasibility
    cuts, so the run aborts with 'dual sub problem unbounded'. Which master
    optimum CPLEX picks is path dependent, so a run that does converge must
    converge to the true optimum mu* = 32; any other outcome is a real bug.
    """
    data = APSPModelData(DATA_DIR, "example1_scaled.xml")
    driver = bounded_Benders.BoundedMagnantiWongBenders("BB_ex1s", data)
    driver.set_iterations(10000)
    driver.solve(False, log_verbosity=0, time_out=120)
    status = driver.get_status()
    if status == "optimal solution":
        assert driver.get_objective() == pytest.approx(32, rel=1e-6), \
            f"BB converged to {driver.get_objective()}, proven optimum 32"
    else:
        assert status == "dual sub problem unbounded", \
            f"expected the documented abort or convergence to 32, got: {status}"
