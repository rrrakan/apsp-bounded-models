"""Correction 3: constraints (17) and (30) must hold for ALL t in T.

Paper, Sec. VI eq. (17) and Sec. VII-C eq. (30): the same-task no-overlap
constraint  s_t + sum_r x_{t,r} e(t,r) <= s_t + mu  (i.e. mu >= e(t, alpha(t)))
is quantified over all t in T. The authors' code skips it for every task that
has execution time 0 on ANY resource (`if 0 not in self.data.TASKS[t].values()`),
and correspondingly restricts the dual variable pi^nd in DSP (32)-(35), the
auxiliary subproblem ASP (36)-(37) and the Benders cut.

tests/data/zero_exec.xml: a0 has e(proc_0)=5, e(proc_1)=0; a1 has e(proc_1)=3.
If a0 is allocated to proc_0 the paper requires mu >= 5. With the skip, nothing
enforces that, and mu can fall to 3 (a1's own (17)).

These tests FAIL on authors_original (mu = 3) and PASS on corrected apsp_models
(mu = 5).
"""
from types import SimpleNamespace

import pytest

from conftest import DATA_DIR

from models.utilities import APSPModelData

TOL = 1e-6


@pytest.fixture()
def data():
    return APSPModelData(DATA_DIR, "zero_exec.xml")


def pinned_extract():
    """Master trial values fixing a0 -> proc_0 (e=5), a1 -> proc_1 (e=3).

    All y orderings may be 0: the only overlap pair (a0, a1) shares proc_1 and
    a0 is not allocated there, so every (28)/(29) big-M term is slack.
    """
    x_bar = {("a0", "proc_0"): 1, ("a0", "proc_1"): 0, ("a1", "proc_1"): 1}
    y_bar = {}
    B = 2 + 2  # |T| + |R|
    for t1, t2 in [("a0", "a1"), ("a1", "a0")]:
        for k in range(B + 1):
            y_bar[(t1, t2, k)] = 0
    return SimpleNamespace(x=x_bar, y=y_bar)


def test_monolithic_17_applies_to_zero_exec_task(data):
    """(17): with a0 pinned to proc_0, the monolithic MILP must give mu >= 5."""
    from models.bounded.bounded_monolithic import BoundedMonolithic

    model = BoundedMonolithic("corr3_monolithic", data)
    model.setup_model("")
    # Test-only pin: force the allocation the master could legitimately choose.
    x = model.decision_vars.x
    model.model.add(x[("a0", "proc_0")] == 1)
    model.solve(False, time_limit=60)

    mu = model.get_objective()
    assert mu >= 5 - TOL, \
        f"(17) violated: a0 on proc_0 has e=5 but mu = {mu}"
    assert mu == pytest.approx(5, abs=1e-4)


def test_primal_sub_30_applies_to_zero_exec_task(data):
    """(30): SP with x_bar fixing a0 -> proc_0 must give mu >= 5."""
    from models.bounded.bounded_Benders.bounded_primal_sub import BoundedPrimalSub

    sp = BoundedPrimalSub("corr3_sp", data)
    sp.setup_model(pinned_extract())
    sp.solve(False, time_limit=60)

    mu = sp.get_objective()
    assert mu >= 5 - TOL, \
        f"(30) violated: a0 on proc_0 has e=5 but SP mu = {mu}"
    assert mu == pytest.approx(5, abs=1e-4)


def test_dual_sub_strong_duality(data):
    """DSP (32)-(35) must match SP: objective 5, not 3.

    This also catches an inconsistent pi^nd: if (30) held for all tasks but
    pi^nd stayed restricted, DSP's optimum would stay at 3 < 5 and strong
    duality would fail.
    """
    from models.bounded.bounded_Benders.bounded_primal_sub import BoundedPrimalSub
    from models.bounded.bounded_Benders.bounded_dual_sub import BoundedDualSub

    extract = pinned_extract()

    sp = BoundedPrimalSub("corr3_sp2", data)
    sp.setup_model(extract)
    sp.solve(False, time_limit=60)

    dsp = BoundedDualSub("corr3_dsp", data)
    dsp.setup_model(extract)
    dsp.solve(False, time_limit=60)

    assert sp.get_objective() == pytest.approx(5, abs=1e-4)
    assert dsp.get_objective() == pytest.approx(sp.get_objective(), abs=1e-4), \
        "strong duality violated between SP and DSP"
