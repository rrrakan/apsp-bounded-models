"""Correction 4: duplicate dependencies between the same task pair.

Paper, Sec. V (proof of Theorem 1): "If multiple dependencies with different
depths are present, then the one with the lowest depth is dominant."

The authors' parser (sdf3_data_conversion.get_channels_dict) stores channels in a
plain dict keyed by (src, dst), so the LAST channel in the XML wins. This test
FAILS on authors_original (parses depth 3 / 2) and PASSES on the corrected
apsp_models (parses the minimum depth 0 / 1).

tests/data/dup_dependency.xml has:
  a0 -> a1: depth 0 first, depth 3 last  => minimum 0
  a1 -> a0: depth 1 first, depth 2 last  => minimum 1
"""
from conftest import DATA_DIR

import sdf3_python_utilities as sdf3


def test_lowest_depth_dominates():
    data = sdf3.parse_application_graph_from_xml(DATA_DIR + "dup_dependency.xml")
    channels = sdf3.get_channels_dict(data)

    # Lowest depth must dominate (paper Sec. V, Theorem 1 proof).
    assert channels[("a0", "a1")] == 0, \
        f"expected lowest depth 0 for (a0, a1), got {channels[('a0', 'a1')]}"
    assert channels[("a1", "a0")] == 1, \
        f"expected lowest depth 1 for (a1, a0), got {channels[('a1', 'a0')]}"
