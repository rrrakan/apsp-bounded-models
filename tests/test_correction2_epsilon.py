"""Correction 2: Benders convergence tolerance epsilon in Algorithm 1.

Paper, Sec. VIII-B: "Solutions assume an optimality tolerance of 1e-6. ... for
Benders decompositions, this is reflected in the value of epsilon in
Algorithm 1." The engine's default is epsilon = 1e-9 and neither driver
overrides it.

Project decision (CLAUDE.md): pass epsilon = 1e-6 in Bounded Benders ONLY.
Bounded Quinton Benders keeps 1e-9, because its objective is tau = 1/mu ~ 1e-4
and an absolute gap of 1e-6 would allow period errors of tens of units.

This test FAILS on authors_original (Bounded Benders epsilon = 1e-9) and PASSES
on the corrected apsp_models.
"""
from conftest import DATA_DIR

from models.utilities import APSPModelData
from models.bounded.bounded_Benders.bounded_magnanti_benders import BoundedMagnantiWongBenders
from models.bounded.bounded_Quinton_Benders.bounded_quinton_magnanti_benders import (
    BoundedQuintonMagnantiWongBenders,
)


def test_bounded_benders_epsilon_is_1e6():
    """Sec. VIII-B: Algorithm 1 epsilon = 1e-6 for Bounded Benders."""
    data = APSPModelData(DATA_DIR, "zero_exec.xml")
    driver = BoundedMagnantiWongBenders("corr2_bb", data)
    assert driver.epsilon == 1e-6, \
        f"Bounded Benders epsilon is {driver.epsilon}, paper Sec. VIII-B requires 1e-6"


def test_bounded_quinton_benders_epsilon_stays_1e9():
    """CLAUDE.md decision: BQB keeps 1e-9 (tau-scale objective)."""
    data = APSPModelData(DATA_DIR, "zero_exec.xml")
    driver = BoundedQuintonMagnantiWongBenders("corr2_bqb", data)
    assert driver.epsilon == 1e-9, \
        f"Bounded Quinton Benders epsilon is {driver.epsilon}, decision says keep 1e-9"
