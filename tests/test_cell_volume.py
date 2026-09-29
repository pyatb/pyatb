"""Physical volume and native optical response must not depend on basis handedness."""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

import numpy as np
import pytest

from pyatb.tb.tb import tb


@pytest.mark.parametrize("left", [False, True])
@pytest.mark.parametrize("scale", [1.0, 2.0])
def test_physical_volume_preserves_lattice_and_reciprocal_basis(left, scale):
    lattice = np.diag([2.0, 3.0, 4.0])
    if left:
        lattice = lattice[[0, 2, 1]]
    original = lattice.copy()
    model = tb(1, scale, lattice)
    assert model.unit_cell_volume == pytest.approx(24.0 * scale**3)
    np.testing.assert_array_equal(model.lattice_vector, original)
    np.testing.assert_allclose(model.reciprocal_vector @ original.T, np.eye(3))


def _native_libraries():
    # Reuse existing BLAS/LAPACKE; do not install or download during tests.
    if os.environ.get("PYATB_NATIVE_TEST_LIBS"):
        return shlex.split(os.environ["PYATB_NATIVE_TEST_LIBS"])
    wheel_libs = Path(sys.prefix) / "lib" / f"python{sys.version_info.major}.{sys.version_info.minor}" / "site-packages" / "pyatb.libs"
    if wheel_libs.is_dir():
        libraries = [next(wheel_libs.glob(pattern), None)
                     for pattern in ("libopenblas*.so*", "liblapacke*.so*")]
        if all(libraries):
            return [*(str(p) for p in libraries),
                    "-Wl,--disable-new-dtags", f"-Wl,-rpath,{wheel_libs}", f"-Wl,-rpath-link,{wheel_libs}"]
    return ["-lopenblas", "-llapacke"]


def test_native_optical_and_static_response_are_handedness_invariant(tmp_path):
    compiler = shutil.which(os.environ.get("CXX", "g++"))
    if compiler is None:
        pytest.skip("native optical regression requires a C++ compiler")
    root = Path(__file__).resolve().parents[1]
    if not (root / "eigen" / "Eigen" / "Eigen").is_file():
        pytest.skip("initialize the pinned Eigen submodule for native optical tests")
    core = root / "src" / "cpp" / "core"
    binary = tmp_path / "optical_handedness"
    command = [
        compiler, "-std=c++11", "-O1", "-fopenmp", "-ffunction-sections",
        "-fdata-sections", "-Wl,--gc-sections", "-I" + str(core),
        "-I" + str(root / "eigen"),
        str(root / "tests" / "native" / "optical_handedness.cpp"),
        *(str(core / (name + ".cpp")) for name in (
            "base_data", "cell_atom", "tools", "xr_operation",
            "band_structure_solver", "velocity_solver", "optical_conductivity_solver")),
        *_native_libraries(), "-o", str(binary),
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    result = subprocess.run([str(binary)], check=True, capture_output=True, text=True,
                            env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
    from io import StringIO
    rows = np.loadtxt(StringIO(result.stdout))
    right, left = rows[rows[:, 0] == 0], rows[rows[:, 0] == 1]
    np.testing.assert_allclose(right[:, 4], 24.0)
    np.testing.assert_allclose(left[:, 4], 24.0)
    np.testing.assert_allclose(left[:, 1:], right[:, 1:], rtol=1e-11, atol=1e-10)
    assert np.all(rows[:, [5, 9, 13]] > 0)  # dissipative diagonal conductivity
    assert np.all(rows[:, [7, 11, 15]] >= -1e-12)  # Im epsilon, including omega=0
    assert np.all(rows[:, [8, 12, 16]] > 0)  # static epsilon minus identity
