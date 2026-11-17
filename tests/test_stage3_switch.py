"""Stage 3 paper_faithful switch: --legacy reproduces the authors' behaviour.

Runs solve_instance.py in a subprocess with PYTHONHASHSEED=0 in both modes and
checks mu, verification and the Benders iteration counts against the recorded
seed-0 references for THIS machine (CPLEX 22.2, PYTHONHASHSEED=0):

  --legacy  : BB 55, BQB 92  (= authors' original, logs/baseline_seed0_log.txt)
  default   : BB 29, BQB 55  (= corrected,          logs/stage2_final_log.txt)

The counts are an observation on this one instance and machine, not a general
claim about convergence. Objectives (mu = 6494) are seed-independent.
"""
import os
import re
import subprocess
import sys

import pytest

APSP_MODELS = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
INSTANCE = os.path.join("example", "graph__10-1000-0.10-02_0.xml")

RESULT_RE = re.compile(
    r"RESULT model=(?P<model>\w+) mu=(?P<mu>[\d.eE+-]+) iterations=(?P<its>\S+) "
    r"status=(?P<status>.*?) time=(?P<time>[\d.]+) verification=(?P<verif>.+)"
)


def run_solver(extra_args):
    env = dict(os.environ, PYTHONHASHSEED="0")
    proc = subprocess.run(
        [sys.executable, "solve_instance.py", "--model", "all",
         "--instance", INSTANCE, "--time-limit", "3600"] + extra_args,
        cwd=APSP_MODELS, env=env, capture_output=True, text=True, timeout=1800,
    )
    assert proc.returncode == 0, f"solve_instance.py failed:\n{proc.stdout}\n{proc.stderr}"
    results = {m.group("model"): m.groupdict() for m in RESULT_RE.finditer(proc.stdout)}
    assert set(results) == {"BM", "BB", "BQM", "BQB"}, f"missing RESULT lines:\n{proc.stdout}"
    return results


def check_common(results):
    for model, r in results.items():
        assert float(r["mu"]) == pytest.approx(6494, abs=1e-3), \
            f"{model}: mu = {r['mu']} != 6494"
        assert "verification successful" in r["verif"], \
            f"{model}: {r['verif']}"


def test_default_mode_counts():
    results = run_solver([])
    check_common(results)
    assert int(results["BB"]["its"]) == 29, \
        f"default BB iterations {results['BB']['its']} != 29 (seed-0 reference)"
    assert int(results["BQB"]["its"]) == 55, \
        f"default BQB iterations {results['BQB']['its']} != 55 (seed-0 reference)"


def test_legacy_mode_reproduces_seed0_baseline():
    results = run_solver(["--legacy"])
    check_common(results)
    assert int(results["BB"]["its"]) == 55, \
        f"--legacy BB iterations {results['BB']['its']} != 55 (authors' seed-0 baseline)"
    assert int(results["BQB"]["its"]) == 92, \
        f"--legacy BQB iterations {results['BQB']['its']} != 92 (authors' seed-0 baseline)"
