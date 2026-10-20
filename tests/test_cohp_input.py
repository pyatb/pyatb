from pyatb.io.default_input import INPUT, function_switch, need_rR_matrix
from pyatb.io.input import get_general_parameter
from pyatb.main import _cohp_parameters_with_default_input


def test_cohp_input_block_is_registered_without_rr_requirement():
    assert "COHP" in function_switch
    assert "COHP" in INPUT
    assert "COHP" not in need_rR_matrix

    cohp = INPUT["COHP"]
    assert cohp["stru_file"][-1] is None
    assert cohp["method"][-1] == "COHP"
    assert cohp["spin"][-1] == "sum"
    assert cohp["input_file"][-1] == ""
    assert cohp["orbital_dir"][-1] == ""
    assert cohp["kpoint_mode"][-1] is None


def test_cohp_main_defaults_input_file_to_input_path(tmp_path):
    parameters = {"stru_file": "STRU", "input_file": ""}

    result = _cohp_parameters_with_default_input(parameters, str(tmp_path))

    assert result == {"stru_file": "STRU", "input_file": str(tmp_path / "Input")}
    assert parameters == {"stru_file": "STRU", "input_file": ""}


def test_cohp_main_keeps_explicit_input_file(tmp_path):
    parameters = {"stru_file": "STRU", "input_file": "custom/INPUT"}

    result = _cohp_parameters_with_default_input(parameters, str(tmp_path))

    assert result["input_file"] == "custom/INPUT"


def test_variable_length_orbital_selector_stops_at_next_known_parameter():
    data = ["atom_i_orbs", "2s", "2p", "atom_j_orbs", "all", "de", "0.02"]
    known = {"atom_i_orbs", "atom_j_orbs", "de"}

    assert get_general_parameter("atom_i_orbs", [str, -1, "all"], data, known) == "2s,2p"
    assert get_general_parameter("atom_j_orbs", [str, -1, "all"], data, known) == "all"
    assert get_general_parameter("de", [float, 1, 0.05], data, known) == 0.02


def test_variable_length_kpoint_weights_are_parsed_as_floats():
    data = ["kpoint_weights", "0.25", "0.75", "de", "0.02"]
    known = {"kpoint_weights", "de"}

    assert get_general_parameter("kpoint_weights", [float, -1, None], data, known) == [0.25, 0.75]


def test_omitted_cohp_kpoint_weights_use_uniform_default(tmp_path):
    # Follow the actual parser -> main adapter path; None means required to the parser.
    parsed = get_general_parameter('kpoint_weights', INPUT['COHP']['kpoint_weights'],
                                   ['de', '0.02'], INPUT['COHP'].keys())
    import numpy as np
    for weights in (parsed, np.asarray(parsed)):
        params = _cohp_parameters_with_default_input({'kpoint_weights': weights}, str(tmp_path))
        assert params.get('kpoint_weights') is None


def test_explicit_empty_kpoint_weights_are_not_silently_defaulted():
    import pytest
    with pytest.raises((ValueError, KeyError), match='kpoint_weights'):
        get_general_parameter('kpoint_weights', INPUT['COHP']['kpoint_weights'],
                              ['kpoint_weights', 'de', '0.02'], INPUT['COHP'].keys())
