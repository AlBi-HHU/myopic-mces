from unittest.mock import MagicMock
import pytest

from myopic_mces.filter_MCES import ComputationMode
from myopic_mces.prepare_input import _classify_db_chunk

# checks if logic of recomputation is correct
# not if DB is correct (this is assumed)
@pytest.mark.parametrize("db_row, always_stronger, expected_mode, is_cached", [
    # 1) Exact distance below threshold -> ComputationMode.EXACT
    ({"mces_exact": 3.0, "best_lower_bound": 3.0, "mces_match2_lower": 3.0}, True, ComputationMode.EXACT, True),
    # 2) Stronger bound above threshold -> ComputationMode.STRONGEST_BOUND
    ({"mces_exact": None, "best_lower_bound": 15.0, "mces_match2_lower": 15.0}, True, ComputationMode.STRONGEST_BOUND, True),
    # 3) Dynamic bound above threshold -> ComputationMode.UNKNOWN
    ({"mces_exact": None, "best_lower_bound": 15.0, "mces_match2_lower": None}, False, ComputationMode.UNKNOWN, True),
    # 4) Inconclusive (e.g. exact distance exceeds threshold) -> must be recomputed
    ({"mces_exact": 12.0, "best_lower_bound": 5.0, "mces_match2_lower": 5.0}, True, None, False),
])
def test_db_modes_emulate_computation_modes(db_row, always_stronger, expected_mode, is_cached):
    s1, s2 = "CC", "CCC"
    ik1, ik2 = "IK1", "IK2"
    client = MagicMock()
    client.get_pair_bounds.return_value = [{"inchikey_a": ik1, "inchikey_b": ik2, **db_row}]

    chunk = [[0, 0, 1]]
    filtered, pre_idx, pre_mces = _classify_db_chunk(
        chunk=chunk,
        smiles_list=[s1, s2],
        smiles_to_inchikey={s1: ik1, s2: ik2},
        client=client,
        threshold=10.0,
        always_stronger_bound=always_stronger,
    )

    if is_cached:
        assert len(pre_mces) == 1
        assert pre_mces[0][3] == expected_mode.value
        assert len(filtered) == 0
    else:
        assert len(pre_mces) == 0
        assert len(filtered) == 1
