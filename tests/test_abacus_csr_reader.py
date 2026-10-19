"""Reader contracts derived from ABACUS legacy and HContainer CSR layouts."""
import numpy as np
import pytest

from pyatb.constants import Ry_to_eV, Ang_to_Bohr
from pyatb.io.abacus_read_xr import abacus_readHR, abacus_readSR, abacus_readrR


def matrix_text(modern, complex_values=False):
    values = '(1,0) (2, 3)\n(2,-3) (4,0)' if complex_values else '1 2\n3 4'
    if not modern:
        values = values.replace('\n', ' ')
        # Legacy zero-nnz blocks have no CSR payload.
        return ('STEP: 0\nMatrix Dimension of H(R): 2\nMatrix number of H(R): 2\n'
                '0 0 0 0\n1 0 -1 4\n' + values + '\n0 1 0 1\n0 2 4\n')
    return (' --- Ionic Step 1 ---\n # print H matrix in real space H(R)\n'
            ' 1 # number of spin directions\n 1 # spin index\n'
            ' 2 # number of localized basis\n 2 # number of Bravais lattice vector R\n\n'
            'user_defined_lattice\n1\n1 0 0\n0 1 0\n0 0 1\nSi\n2\nDirect\n0 0 0\n0.25 0.25 0.25\n'
            '# CSR Format\n# explanation\n\n0 0 0 0\n# CSR values\n\n'
            '# CSR column indices\n\n# CSR row pointers\n0 0 0\n'
            '1 0 -1 4\n# CSR values\n' + values + '\n# CSR column indices\n0 1\n0 1\n'
            '# CSR row pointers\n0 2\n4\n')


@pytest.mark.parametrize('modern', [False, True])
@pytest.mark.parametrize('nspin', [1, 2, 4])
@pytest.mark.parametrize('kind', ['H', 'S'])
def test_reads_both_formats_with_units_spin_and_zero_block(tmp_path, modern, nspin, kind):
    path = tmp_path / 'arbitrary-name.csr'
    path.write_text(matrix_text(modern, nspin == 4))
    if kind == 'H':
        result = abacus_readHR(nspin, str(path), 'Ry')
        factor = Ry_to_eV
    else:
        result = abacus_readSR(nspin, str(path))
        factor = 1
    assert result.des == kind and result.basis_num == 2 and result.R_num == 2
    np.testing.assert_array_equal(result.R_direct_coor, [[0, 0, 0], [1, 0, -1]])
    expected = np.array([[0, 0, 0], [1, 2+3j if nspin == 4 else 2, 4]]) * factor
    np.testing.assert_allclose(result.XR.toarray(), expected, rtol=0, atol=0)


@pytest.mark.parametrize('change', ['empty', 'truncated', 'extra_values', 'bad_columns', 'bad_pointers'])
def test_rejects_malformed_modern_matrix(tmp_path, change):
    text = matrix_text(True)
    if change == 'empty': text = ''
    if change == 'truncated': text = text.rsplit('4', 1)[0]
    if change == 'extra_values': text = text.replace('1 2\n3 4', '1 2 3 4 5')
    if change == 'bad_columns': text = text.replace('0 1\n0 1', '0 2\n0 1')
    if change == 'bad_pointers': text = text.replace('0 2\n4', '0 3\n2')
    path = tmp_path / 'bad.csr'
    path.write_text(text)
    with pytest.raises(ValueError):
        abacus_readHR(1, str(path), 'eV')


def test_legacy_position_matrix_units_and_empty_component(tmp_path):
    path = tmp_path / 'rr.csr'
    path.write_text('STEP: 0\nMatrix Dimension of r(R): 2\nMatrix number of r(R): 1\n'
                    '0 0 0\n2\n1 2\n0 1\n0 1 2\n0\n2\n3 4\n0 1\n0 1 2\n')
    result = abacus_readrR(str(path), 'Bohr')
    for matrix, expected in zip(result, [[1, 0, 2], [0, 0, 0], [3, 0, 4]]):
        np.testing.assert_allclose(matrix.XR.toarray()[0], np.array(expected)/Ang_to_Bohr)


@pytest.mark.parametrize('kind', ['H', 'S'])
def test_real_si_both_formats_preserve_every_matrix_element(tmp_path, kind):
    import gzip
    from pathlib import Path

    fixtures = Path(__file__).parent / 'fixtures/abacus-csr'
    names = (['develop-hrs1_nao.csr.gz', 'lts-data-HR-sparse_SPIN0.csr.gz'] if kind == 'H'
             else ['develop-sr_nao.csr.gz', 'lts-data-SR-sparse_SPIN0.csr.gz'])
    matrices = []
    for name in names:
        path = tmp_path / name.removesuffix('.gz')
        path.write_bytes(gzip.decompress((fixtures / name).read_bytes()))
        result = abacus_readHR(1, str(path), 'Ry') if kind == 'H' else abacus_readSR(1, str(path))
        assert result.basis_num == 26
        matrices.append({tuple(r): data for r, data in zip(result.R_direct_coor, result.XR.toarray())})
    # Develop omits 30 all-zero S blocks retained by the legacy producer.
    assert [len(m) for m in matrices] == ([93, 93] if kind == 'H' else [63, 93])
    for r in matrices[0].keys() | matrices[1].keys():
        np.testing.assert_array_equal(matrices[0].get(r, np.zeros(351)), matrices[1].get(r, np.zeros(351)))


@pytest.mark.parametrize('kind', ['H', 'S'])
def test_legacy_header_allows_leading_whitespace(tmp_path, kind):
    path = tmp_path / 'indented.csr'
    path.write_text(matrix_text(False).replace('Matrix', '   Matrix'))
    result = abacus_readHR(1, str(path), 'eV') if kind == 'H' else abacus_readSR(1, str(path))
    np.testing.assert_array_equal(result.XR.toarray()[1], [1, 2, 4])


@pytest.mark.parametrize('nspin', [1, 2, 4])
@pytest.mark.parametrize('sparse', [False, True])
def test_init_tb_aligns_abacus_R_union_without_losing_values(tmp_path, monkeypatch, nspin, sparse):
    import importlib
    module = importlib.import_module('pyatb.init_tb')
    solver_module = importlib.import_module('pyatb.tb.solver')

    class Capture:
        def __init__(self, *args): self.hsr = None
        def set_HSR(self, *args): self.hsr = args
        def set_HSR_sparse(self, *args): self.hsr = args
        def set_rR(self, *args): self.rr = args
        def set_rR_sparse(self, *args): self.rr = args

    monkeypatch.setattr(solver_module, 'tb_solver_', Capture)
    hr = tmp_path / 'h.csr'
    sr = tmp_path / 's.csr'
    rr = tmp_path / 'rr.csr'
    hr.write_text(matrix_text(True, nspin == 4))
    # The S block has a distinct R to exercise the union, not only padding to H.
    sr.write_text('Matrix Dimension of S(R): 2\nMatrix number of S(R): 1\n'
                  '2 0 0 2\n' + ('(5,0) (6,0)' if nspin == 4 else '5 6') + '\n0 1\n0 1 2\n')
    rr.write_text('Matrix Dimension of r(R): 2\nMatrix number of r(R): 1\n'
                  '3 0 0\n1\n7\n0\n0 1 1\n0\n0\n')
    model = module.init_tb(nspin=nspin, HR_route=[str(hr), str(hr)] if nspin == 2 else str(hr),
                           SR_route=str(sr), need_rR=True, rR_route=str(rr), isSparse=sparse)
    engines = [model.tb_solver_up, model.tb_solver_dn] if nspin == 2 else [model.tb_solver]
    for engine in engines:
        count, coords, basis, h, s = engine.tb_solver.hsr
        assert count == 4 and basis == 2
        np.testing.assert_array_equal(coords, [[0, 0, 0], [1, 0, -1], [2, 0, 0], [3, 0, 0]])
        dense = lambda value: value.toarray() if sparse else value
        np.testing.assert_array_equal(dense(h), [[0,0,0],[1,2+3j if nspin==4 else 2,4],[0,0,0],[0,0,0]])
        np.testing.assert_array_equal(dense(s), [[0,0,0],[0,0,0],[5,0,6],[0,0,0]])
        np.testing.assert_array_equal(dense(engine.tb_solver.rr[0]), [[0,0,0],[0,0,0],[0,0,0],[7,0,0]])


def test_alignment_keeps_inputs_and_rejects_ambiguous_coordinates():
    from scipy.sparse import csc_matrix
    from pyatb.tb.multixr import multiXR
    from pyatb.io.abacus_read_xr import _align_abacus_R

    h, s = multiXR('H'), multiXR('S')
    h.set_XR(1, np.array([[0,0,0]]), 1, csc_matrix([[2.0]]))
    s.set_XR(1, np.array([[1,0,0]]), 1, csc_matrix([[3.0]]))
    aligned = _align_abacus_R(h, s)
    assert h.R_num == s.R_num == 1
    np.testing.assert_array_equal(h.XR.toarray(), [[2]])
    np.testing.assert_array_equal(s.XR.toarray(), [[3]])
    np.testing.assert_array_equal(aligned[0].XR.toarray(), [[2],[0]])
    np.testing.assert_array_equal(aligned[1].XR.toarray(), [[0],[3]])
    s.set_XR(2, np.array([[0,0,0],[0,0,0]]), 1, csc_matrix([[3.0],[4.0]]))
    with pytest.raises(ValueError, match='duplicate R'):
        _align_abacus_R(h, s)
    s.set_XR(1, np.array([[0,0,0]]), 2, csc_matrix([[1,0,1]]))
    with pytest.raises(ValueError, match='basis dimensions'):
        _align_abacus_R(h, s)
