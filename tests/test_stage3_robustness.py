"""Stage 3 robustness fixes in the Benders engine.

(a) solve(time_out=None) must mean "no time limit". The original engine computed
    `time_out - (time.time() - start_time)` unconditionally and crashed with a
    TypeError on None (silently caught as "Error during solve", then crashing
    again inside write_error_log).
(b) write_error_log must not crash when best_master_problem_extract is None
    (i.e. when the error occurred before the first completed iteration).

Both tests FAIL on authors_original and PASS on apsp_models.
"""
from conftest import DATA_DIR

from models.utilities import APSPModelData
from models.bounded.bounded_Benders.bounded_magnanti_benders import BoundedMagnantiWongBenders


def test_solve_with_time_out_none(tmp_path):
    """time_out=None runs to optimality instead of crashing.

    zero_exec.xml optimum: a0 -> proc_1 (e=0), a1 -> proc_1 (e=3) gives mu = 3;
    allocating a0 to proc_0 would give max(5, 3) = 5, so 3 is optimal.
    """
    data = APSPModelData(DATA_DIR, "zero_exec.xml")
    driver = BoundedMagnantiWongBenders("rob_timeout_none", data)
    driver.set_iterations(1000)
    driver.solve(False, log_verbosity=0, time_out=None,
                 log_destination=str(tmp_path) + "\\", log_file_name="rob")

    assert driver.get_status() == "optimal solution", \
        f"expected optimal solution with time_out=None, got: {driver.get_status()}"
    assert abs(driver.get_objective() - 3) < 1e-4, \
        f"expected mu = 3 on zero_exec.xml, got {driver.get_objective()}"


def test_write_error_log_with_no_best_extract(tmp_path):
    """write_error_log with best_master_problem_extract=None must not raise.

    This is the engine state when an exception occurs before any iteration
    completes; the original crashed on primal_sub.setup_constraints(None).
    """
    data = APSPModelData(DATA_DIR, "zero_exec.xml")
    driver = BoundedMagnantiWongBenders("rob_errorlog", data)
    assert driver.best_master_problem_extract is None  # state before any iteration

    # Must not raise (original: AttributeError 'NoneType' object has no attribute 'x').
    driver.write_error_log(str(tmp_path) + "\\", "rob")
