import numpy as np
import json
import pytest
import inspect
import os
import types

import pyatb.fermi.cohp as cohp_module
from pyatb.fermi.cohp import COHP, cohp_values_one_k


def test_cohp_values_one_k_matches_standalone_formula():
    matrix = np.array(
        [
            [1.0, 2.0 + 0.5j, 0.0],
            [2.0 - 0.5j, 3.0, 4.0],
            [0.0, 4.0, 5.0],
        ],
        dtype=np.complex128,
    )
    eigenvalues = np.array([-1.0, 0.5], dtype=float)
    eigenvectors = np.array(
        [
            [1.0 + 0.0j, 0.0 + 1.0j],
            [0.5 + 0.5j, 1.0 + 0.0j],
            [0.25 + 0.0j, 0.5 - 0.25j],
        ],
        dtype=np.complex128,
    )

    energies, values = cohp_values_one_k(
        matrix, eigenvalues, eigenvectors, [0, 1], [2], pair_factor=2.0
    )

    expected = []
    for ib in range(eigenvectors.shape[1]):
        total = 0.0
        for iorb in [0, 1]:
            total += (eigenvectors[iorb, ib].conjugate() * matrix[iorb, 2] * eigenvectors[2, ib]).real
        expected.append(2.0 * total)

    np.testing.assert_allclose(energies, eigenvalues)
    np.testing.assert_allclose(values, np.array(expected))


def test_atom_blocks_reconstruct_hamiltonian_expectation_one_k():
    matrix = np.array([[1.0, 0.3 - 0.2j], [0.3 + 0.2j, 2.0]], dtype=np.complex128)
    eigenvectors = np.array([[0.8 + 0.1j], [0.4 - 0.3j]], dtype=np.complex128)
    eigenvectors /= np.linalg.norm(eigenvectors)
    eigenvalues = np.array([0.0])

    onsite_i = cohp_values_one_k(matrix, eigenvalues, eigenvectors, [0], [0])[1][0]
    onsite_j = cohp_values_one_k(matrix, eigenvalues, eigenvectors, [1], [1])[1][0]
    offsite = cohp_values_one_k(
        matrix, eigenvalues, eigenvectors, [0], [1], pair_factor=2.0
    )[1][0]
    expected = (eigenvectors[:, 0].conj() @ matrix @ eigenvectors[:, 0]).real
    np.testing.assert_allclose(onsite_i + onsite_j + offsite, expected)


def test_print_plot_script_only_writes_script(tmp_path):
    cohp = COHP.__new__(COHP)
    cohp.output_path = str(tmp_path)

    cohp.print_plot_script()

    assert (tmp_path / "plot_cohp.py").is_file()
    assert not (tmp_path / "cohp.pdf").exists()
    assert "subprocess.run" not in inspect.getsource(COHP.print_plot_script)


class _SingleKGenerator:
    total_kpoint_num = 1

    def __iter__(self):
        return iter([np.array([[0.0, 0.0, 0.0]], dtype=float)])


class _UnitSolver:
    def diago_H(self, kpoints):
        eigenvectors = np.ones((kpoints.shape[0], 1, 1), dtype=np.complex128)
        eigenvalues = np.zeros((kpoints.shape[0], 1), dtype=float)
        return eigenvectors, eigenvalues

    def get_Hk(self, kpoints):
        return np.ones((kpoints.shape[0], 1, 1), dtype=np.complex128)


class _TimerStub:
    """Stand-in for the runtime timer (upstream 7b4832e defers it to initialize_runtime)."""

    def start(self, *_args, **_kwargs):
        return None

    def end(self, *_args, **_kwargs):
        return None


def _cohp_with_fake_solver(nspin):
    cohp = COHP.__new__(COHP)
    cohp.nspin = nspin
    cohp._COHP__k_generator = _SingleKGenerator()
    solver = _UnitSolver()
    cohp._COHP__tb_solver = (solver, solver) if nspin == 2 else (solver,)
    return cohp


def _calculate_fake_spectrum(cohp, monkeypatch):
    captured = {}
    selection = types.SimpleNamespace(
        atom_i_orbitals=[0],
        atom_j_orbitals=[0],
        orbital_map=types.SimpleNamespace(total_orbitals=1, atoms=[]),
    )
    monkeypatch.setattr(cohp_module, "resolve_cohp_orbitals", lambda **kwargs: selection)
    monkeypatch.setattr(COHP, "print_data", lambda self, output_prefix, energy, spectrum, *args: captured.update(spectrum=spectrum.copy()))
    monkeypatch.setattr(cohp_module, "timer", _TimerStub())
    monkeypatch.setattr(cohp_module, "RUNNING_LOG", os.devnull)

    cohp.calculate_cohp(
        fermi_energy=0.0,
        stru_file="STRU",
        atom_i_index=1,
        atom_j_index=1,
        method="COHP",
        spin="sum",
        e_range=[-1.0, 1.0],
        de=0.5,
        sigma=0.1,
        invert=0,
    )
    return captured["spectrum"]


def test_nspin1_sum_cohp_includes_spin_degeneracy(monkeypatch):
    nspin1_spectrum = _calculate_fake_spectrum(_cohp_with_fake_solver(nspin=1), monkeypatch)
    nspin2_sum_spectrum = _calculate_fake_spectrum(_cohp_with_fake_solver(nspin=2), monkeypatch)

    np.testing.assert_allclose(nspin1_spectrum, nspin2_sum_spectrum)


class _TwoKGenerator:
    total_kpoint_num = 2

    def __iter__(self):
        return iter([np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]], dtype=float)])


class _WeightedSolver(_UnitSolver):
    def get_Hk(self, kpoints):
        values = np.where(kpoints[:, 0] == 0.0, 1.0, 3.0)
        return values[:, None, None].astype(np.complex128)


def test_calculate_cohp_uses_explicit_nonuniform_kpoint_weights(monkeypatch):
    cohp = COHP.__new__(COHP)
    cohp.nspin = 2
    cohp._COHP__k_generator = _TwoKGenerator()
    cohp._COHP__tb_solver = (_WeightedSolver(), _WeightedSolver())
    captured = {}
    selection = types.SimpleNamespace(
        atom_i_orbitals=[0],
        atom_j_orbitals=[0],
        orbital_map=types.SimpleNamespace(total_orbitals=1, atoms=[]),
    )
    monkeypatch.setattr(cohp_module, "resolve_cohp_orbitals", lambda **kwargs: selection)
    monkeypatch.setattr(
        COHP,
        "print_data",
        lambda self, output_prefix, energy, spectrum, *args: captured.update(spectrum=spectrum.copy()),
    )
    monkeypatch.setattr(cohp_module, "timer", _TimerStub())
    monkeypatch.setattr(cohp_module, "RUNNING_LOG", os.devnull)

    cohp.calculate_cohp(
        fermi_energy=0.0,
        stru_file="STRU",
        atom_i_index=1,
        atom_j_index=1,
        method="COHP",
        spin="up",
        e_range=[-1.0, 1.0],
        de=0.5,
        sigma=0.1,
        invert=0,
        kpoint_weights=[0.25, 0.75],
    )

    expected = np.zeros(5)
    cohp._accumulate_one_k(
        expected,
        np.array([[2.5]], dtype=np.complex128),
        np.array([0.0]),
        np.array([[1.0]], dtype=np.complex128),
        [0],
        [0],
        -1.0,
        0.5,
        0.1,
        False,
        pair_factor=1.0,
        kpoint_weight=1.0,
    )
    np.testing.assert_allclose(captured["spectrum"], expected)

# End-to-end output regressions use analytic eigenpairs, not legacy curves.
class _AnalyticSolver:
    def diago_H(self, kpoints):
        vectors = np.array([[1., 1.], [-1., 1.]]) / np.sqrt(2.)
        return (np.broadcast_to(vectors, (len(kpoints), 2, 2)),
                np.broadcast_to([-1., 1.], (len(kpoints), 2)))

    def get_Hk(self, kpoints):
        return np.broadcast_to([[0., 1.], [1., 0.]], (len(kpoints), 2, 2))

    def get_Sk(self, kpoints):
        return np.broadcast_to(np.eye(2), (len(kpoints), 2, 2))


def _analytic_output(tmp_path, monkeypatch, selection=([0], [1]), **overrides):
    cohp = _cohp_with_fake_solver(overrides.pop('nspin', 2))
    solver = overrides.pop("solver", _AnalyticSolver())
    cohp._COHP__tb_solver = (solver, solver)
    cohp._COHP__k_generator = overrides.pop("generator", _SingleKGenerator())
    cohp.output_path = str(tmp_path)
    resolved = types.SimpleNamespace(
        atom_i_orbitals=selection[0], atom_j_orbitals=selection[1],
        orbital_map=types.SimpleNamespace(total_orbitals=2, atoms=[]))
    monkeypatch.setattr(cohp_module, 'resolve_cohp_orbitals', lambda **kw: resolved)
    monkeypatch.setattr(cohp_module, 'timer', _TimerStub())
    monkeypatch.setattr(cohp_module, 'RUNNING_LOG', os.devnull)
    params = dict(fermi_energy=0., stru_file='STRU', atom_i_index=1,
                  atom_j_index=1, spin='up', e_range=[-2., 2.], de=.01,
                  sigma=.1, invert=0)
    params.update(overrides)
    cohp.calculate_cohp(**params)
    return cohp, np.loadtxt(tmp_path / 'COHP.dat'), json.loads((tmp_path / 'COHP.meta.json').read_text())


def test_same_atom_disjoint_channels_have_full_hermitian_weight(tmp_path, monkeypatch):
    _, data, meta = _analytic_output(tmp_path, monkeypatch)
    # The occupied bonding eigenstate has C0*C1=-1/2; both directions sum to -1.
    assert np.sum(data[data[:, 0] < 0, 1]) * .01 == pytest.approx(-1., abs=1e-7)
    assert meta['pair_factor'] == 2.


def test_equal_orbital_sets_are_order_independent(tmp_path, monkeypatch):
    _, data, meta = _analytic_output(tmp_path, monkeypatch, selection=([0, 1], [1, 0]))
    assert np.sum(data[data[:, 0] < 0, 1]) * .01 == pytest.approx(-1., abs=1e-7)
    assert meta['pair_factor'] == 1.


@pytest.mark.parametrize('selection', [([0], [0, 1]), ([], [1])])
def test_ambiguous_or_empty_orbital_selections_fail(tmp_path, monkeypatch, selection):
    with pytest.raises(ValueError, match='overlap|non-empty'):
        _analytic_output(tmp_path, monkeypatch, selection=selection)


@pytest.mark.parametrize('method', ['COHP', 'COOP'])
@pytest.mark.parametrize('invert', [0, 1])
def test_output_records_sign_and_discrete_integral(tmp_path, monkeypatch, method, invert):
    cohp, data, meta = _analytic_output(tmp_path, monkeypatch, method=method,
                                       invert=invert, selection=([0, 1], [0, 1]))
    label = ('-' if invert else '') + method
    assert (tmp_path / 'COHP.dat').read_text().splitlines()[0] == '# energy_eV ' + label
    assert meta['invert'] is bool(invert)
    assert meta['spectrum_label'] == label
    integral = meta['occupied_integral']
    assert integral['method'] == 'I' + method
    assert integral['occupation'] == 'E <= EF (zero temperature)'
    assert integral['unit'] == ('eV' if method == 'COHP' else 'dimensionless')
    expected = -1. if method == 'COHP' else 1.
    assert integral['raw_value'] == pytest.approx(expected)
    assert integral['value'] == pytest.approx(expected * (-1 if invert else 1))
    cohp.print_plot_script()
    # Execute generated plotting code and inspect the actual label it supplies.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.axes
    labels = []
    original = matplotlib.axes.Axes.set_xlabel
    def capture(self, label, *args, **kwargs):
        labels.append(label)
        return original(self, label, *args, **kwargs)
    monkeypatch.setattr(matplotlib.axes.Axes, 'set_xlabel', capture)
    monkeypatch.chdir(tmp_path)
    exec(compile((tmp_path / 'plot_cohp.py').read_text(), 'plot_cohp.py', 'exec'), {})
    assert label in labels
    assert (tmp_path / 'cohp.pdf').is_file()


@pytest.mark.parametrize('nspin,spin,factor', [(1, 'sum', 2), (2, 'sum', 2), (2, 'up', 1), (2, 'down', 1)])
def test_integral_ignores_plot_window_and_includes_fermi_boundary(tmp_path, monkeypatch, nspin, spin, factor):
    _, _, meta = _analytic_output(tmp_path, monkeypatch, nspin=nspin, spin=spin,
                                  fermi_energy=-1., e_range=[2., 3.], de=.2, sigma=.3)
    assert meta['occupied_integral']['raw_value'] == pytest.approx(-factor)


@pytest.mark.parametrize('energy', [0., 1., -.05, 1.05])
def test_gaussian_includes_boundary_states_and_external_tails(energy):
    spectrum = np.zeros(11)
    cohp = COHP.__new__(COHP)
    cohp._accumulate_one_k(spectrum, np.ones((1, 1)), np.array([energy]),
                           np.ones((1, 1)), [0], [0], 0., .1, .1, False)
    grid = np.arange(11) * .1
    expected = np.exp(-((grid-energy)/.1)**2) / (np.sqrt(np.pi)*.1)
    np.testing.assert_allclose(spectrum, expected, atol=1e-40)


@pytest.mark.parametrize('params', [dict(de=0), dict(de=-.1), dict(de=np.nan),
    dict(sigma=0), dict(sigma=np.inf), dict(e_range=[1., -1.]),
    dict(e_range=[np.nan, 1.]), dict(e_range=[0., 1., 2.]),
    dict(fermi_energy=np.inf), dict(de=5.)])
def test_invalid_energy_grid_fails_clearly(tmp_path, monkeypatch, params):
    with pytest.raises(ValueError, match='finite|positive|range|grid'):
        _analytic_output(tmp_path, monkeypatch, **params)


def test_kpoint_weights_are_stable_for_large_finite_inputs():
    np.testing.assert_allclose(COHP._normalized_kpoint_weights(2, [1e308, 1e308]), [.5, .5])
    with pytest.raises(ValueError, match='positive'):
        COHP._normalized_kpoint_weights(0)


class _ChunkedKGenerator:
    total_kpoint_num = 2

    def __iter__(self):
        return iter([np.array([[0., 0., 0.]]), np.array([[.5, 0., 0.]])])


class _WeightedAnalyticSolver(_AnalyticSolver):
    def diago_H(self, kpoints):
        vectors, energies = super().diago_H(kpoints)
        scale = np.where(kpoints[:, 0] == 0., 1., 3.)
        return vectors, energies * scale[:, None]

    def get_Hk(self, kpoints):
        scale = np.where(kpoints[:, 0] == 0., 1., 3.)
        return super().get_Hk(kpoints) * scale[:, None, None]


def test_integral_and_spectrum_preserve_weights_across_chunks(tmp_path, monkeypatch):
    results = []
    for generator in (_TwoKGenerator(), _ChunkedKGenerator()):
        _, data, meta = _analytic_output(tmp_path, monkeypatch, generator=generator,
            solver=_WeightedAnalyticSolver(), kpoint_weights=[1., 3.], e_range=[-4., 4.])
        assert meta['occupied_integral']['raw_value'] == pytest.approx(-2.5)
        results.append(data)
    np.testing.assert_allclose(*results)
    assert results[0][results[0][:, 0] < 0, 1].sum() * .01 == pytest.approx(-2.5, abs=1e-7)
