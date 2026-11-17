"""Single-instance run script for the four bounded APSP models.

Usage (from apsp_models/, inside the venv):
    python solve_instance.py --model all --instance example/graph__10-1000-0.10-02_0.xml --time-limit 3600
    python solve_instance.py --model BB --legacy ...

--legacy sets paper_faithful=False: the authors' release-1.0.2 core-point logic
and epsilon = 1e-9 in Bounded Benders (corrections 3-4 stay on in both modes).
For run-to-run reproducible Benders iteration counts set PYTHONHASHSEED=0
(resource sets iterate in hash order); objectives do not depend on it.
"""
import argparse
import os
import sys

import pandas as pd

import models.bounded.bounded_monolithic as bounded_monolithic
import models.bounded.bounded_Benders as bounded_Benders
import models.bounded.bounded_Quinton_monolithic as bounded_Quinton_monolithic
import models.bounded.bounded_Quinton_Benders as bounded_Quinton_Benders
from models.utilities import APSPModelData
from verify_instance import verify_solution

MODELS = ["BM", "BB", "BQM", "BQB"]

parser = argparse.ArgumentParser(description="Solve one APSP instance with the bounded models (CPLEX).")
parser.add_argument("--model", choices=MODELS + ["all"], default="all")
parser.add_argument("--instance", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "example", "graph__10-1000-0.10-02_0.xml"))
parser.add_argument("--time-limit", type=float, default=3600, help="time limit in seconds per model (default 3600)")
parser.add_argument("--legacy", action="store_true", help="paper_faithful=False: authors' core point and epsilon=1e-9 in Bounded Benders")
args = parser.parse_args()

if args.time_limit <= 0:
    sys.exit("--time-limit must be positive (the Benders models always need a time limit)")

if os.environ.get("PYTHONHASHSEED") != "0":
    print("Note: PYTHONHASHSEED is not 0; Benders iteration counts may vary between runs (results do not).")

paper_faithful = not args.legacy
selected = MODELS if args.model == "all" else [args.model]

instance = os.path.abspath(args.instance)
filepath = os.path.dirname(instance) + os.sep
filename = os.path.basename(instance)
graph = filename[:-4]

timeout = args.time_limit
iterations = 1000000

# Per-instance output folder for solution exports (consumed by the verifier).
results_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results") + os.sep
instance_results_path = results_path + graph + os.sep
if not os.path.isdir(instance_results_path):
    os.makedirs(instance_results_path)

# Parse instance data from XML.
data = APSPModelData(filepath, filename)

rows = {}


def report(key, mu, status, time_s, its, verification):
    rows[key] = (mu, status, time_s, its, verification)
    its_str = "NA" if its is None else its
    print(f"RESULT model={key} mu={mu} iterations={its_str} status={status} "
          f"time={time_s:.2f} verification={verification}")


if "BM" in selected:
    print('Bounded Monolithic')
    model = bounded_monolithic.BoundedMonolithic('Bounded Monolithic', data)
    model.setup_model("")
    model.solve(False, time_limit=timeout)
    model.write_solution(instance_results_path, graph + '_Bounded_Monolithic')
    v = verify_solution("Bounded_Monolithic", filepath, graph, results_path)
    report("BM", model.get_objective(), model.get_status(), model.get_solve_time(), None, v)

if "BB" in selected:
    print('Bounded Benders')
    model = bounded_Benders.BoundedMagnantiWongBenders('Bounded Benders', data, paper_faithful=paper_faithful)
    model.set_iterations(iterations)
    model.solve(False, log_verbosity=0, time_out=timeout)
    model.write_solution(instance_results_path, graph + '_Bounded_Benders')
    v = verify_solution("Bounded_Benders", filepath, graph, results_path)
    report("BB", model.get_objective(), model.get_status(), model.get_solve_time(), model.iteration, v)

if "BQM" in selected:
    print('Bounded Quinton Monolithic')
    model = bounded_Quinton_monolithic.BoundedQuintonMonolithic('Bounded Quinton Monolithic', data)
    model.setup_model("")
    model.solve(False, time_limit=timeout)
    model.write_solution(instance_results_path, graph + '_Bounded_Quinton_Monolithic')
    v = verify_solution("Bounded_Quinton_Monolithic", filepath, graph, results_path)
    # The Quinton models optimise tau = 1/mu, so the period is 1/get_objective().
    report("BQM", 1 / model.get_objective(), model.get_status(), model.get_solve_time(), None, v)

if "BQB" in selected:
    print('Bounded Quinton Benders')
    model = bounded_Quinton_Benders.BoundedQuintonMagnantiWongBenders('Bounded Quinton Benders', data, paper_faithful=paper_faithful)
    model.set_iterations(iterations)
    model.solve(False, log_verbosity=0, time_out=timeout)
    model.write_solution(instance_results_path, graph + '_Bounded_Quinton_Benders')
    v = verify_solution("Bounded_Quinton_Benders", filepath, graph, results_path)
    # The Quinton models optimise tau = 1/mu, so the period is 1/get_objective().
    report("BQB", 1 / model.get_objective(), model.get_status(), model.get_solve_time(), model.iteration, v)

# Summary table (one row per model).
RESULTS = pd.DataFrame(
    [(k,) + rows[k] for k in MODELS if k in rows],
    columns=["Model", "Objective (mu)", "Status", "Time (s)", "Iterations", "Verification"],
)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)
print()
print(f"Instance: {filename} | mode: {'legacy (paper_faithful=False)' if args.legacy else 'paper_faithful'} | time limit: {timeout}")
print(RESULTS.to_string(index=False))
