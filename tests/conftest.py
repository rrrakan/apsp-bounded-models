"""Test configuration.

The package under test is selected with the APSP_ROOT environment variable so the
same tests can be run against the authors' unmodified release (expected to FAIL on
the correction tests) and against apsp_models (expected to PASS):

    pytest tests\test_correction4_lowest_depth.py                 # tests apsp_models
    $env:APSP_ROOT='C:\apsp_project\authors_original'; pytest ... # tests the original

Default: the apsp_models folder that contains this tests/ directory.
"""
import os
import sys

ROOT = os.environ.get(
    "APSP_ROOT",
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..")),
)
sys.path.insert(0, ROOT)

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data") + os.sep
