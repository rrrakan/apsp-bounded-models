"""Correction 1: Papadakos core-point update in both Benders drivers.

Paper, Sec. VII-D: the core point v'_k for a variable v in iteration k is
    v'_0 = 1/2 v_0            for k = 0,
    v'_k = 1/2 v'_{k-1} + 1/2 v_k   for k > 0,
where v_k is the master trial value of iteration k (a running average).

The authors' code (a) overwrites the previous core point with the CURRENT trial
point before averaging, so the core point equals the trial point from iteration
1 on, and (b) has a typo `0.5 + prev + 0.5 * cur` (instead of `0.5 * prev`) in
the mu/z component.

These tests drive __calculate_core_point directly with known trial values over
several iterations and compare against the paper's recurrence. They FAIL on
authors_original and PASS on the corrected apsp_models.
"""
from types import SimpleNamespace

import pytest

from conftest import DATA_DIR

from models.utilities import APSPModelData
from models.bounded.bounded_Benders.bounded_magnanti_benders import BoundedMagnantiWongBenders
from models.bounded.bounded_Quinton_Benders.bounded_quinton_magnanti_benders import (
    BoundedQuintonMagnantiWongBenders,
)

CALC = "_MagnantiWongPapadakosBendersDecomposition__calculate_core_point"


def check_recurrence(driver, trials, fields):
    """Feed trial points for iterations 0..n and check v'_k after each update.

    trials: list of dicts {field: value-dict or scalar} per iteration.
    fields: mapping {field_name: is_dict}.
    """
    expected = None  # expected core point per field
    for k, trial in enumerate(trials):
        driver.master_problem_extract = SimpleNamespace(**trial)
        getattr(driver, CALC)(k)

        if k == 0:
            # v'_0 = 1/2 v_0  (paper Sec. VII-D)
            expected = {
                f: ({i: 0.5 * v for i, v in trial[f].items()} if is_dict else 0.5 * trial[f])
                for f, is_dict in fields.items()
            }
        else:
            # v'_k = 1/2 v'_{k-1} + 1/2 v_k  (paper Sec. VII-D)
            expected = {
                f: ({i: 0.5 * expected[f][i] + 0.5 * trial[f][i] for i in trial[f]}
                    if is_dict else 0.5 * expected[f] + 0.5 * trial[f])
                for f, is_dict in fields.items()
            }

        for f, is_dict in fields.items():
            got = getattr(driver.core_point, f)
            if is_dict:
                for i in expected[f]:
                    assert got[i] == pytest.approx(expected[f][i]), (
                        f"iteration {k}, field {f}[{i}]: core point {got[i]} "
                        f"!= paper running average {expected[f][i]}"
                    )
            else:
                assert got == pytest.approx(expected[f]), (
                    f"iteration {k}, field {f}: core point {got} "
                    f"!= paper running average {expected[f]}"
                )


def test_bounded_benders_core_point_running_average():
    data = APSPModelData(DATA_DIR, "zero_exec.xml")
    driver = BoundedMagnantiWongBenders("corr1_bb", data)

    # Distinct trial points per iteration; after iteration 0 the paper's core
    # point differs from the trial point, which exposes the trial-point aliasing.
    trials = [
        {"x": {"A": 1, "B": 0}, "y": {"C": 1}, "mu": 10.0},
        {"x": {"A": 0, "B": 1}, "y": {"C": 0}, "mu": 6.0},
        {"x": {"A": 1, "B": 1}, "y": {"C": 1}, "mu": 4.0},
        {"x": {"A": 0, "B": 0}, "y": {"C": 1}, "mu": 2.0},
    ]
    check_recurrence(driver, trials, {"x": True, "y": True, "mu": False})


def test_bounded_quinton_benders_core_point_running_average():
    data = APSPModelData(DATA_DIR, "zero_exec.xml")
    driver = BoundedQuintonMagnantiWongBenders("corr1_bqb", data)

    trials = [
        {"K": {("A", "B"): 2}, "m": {("A", "r"): 1}, "z": 8.0},
        {"K": {("A", "B"): -1}, "m": {("A", "r"): 0}, "z": 4.0},
        {"K": {("A", "B"): 3}, "m": {("A", "r"): 1}, "z": 2.0},
    ]
    check_recurrence(driver, trials, {"K": True, "m": True, "z": False})
