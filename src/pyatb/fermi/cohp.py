from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.kpt import kpoint_generator
from pyatb.parallel import op_sum
from pyatb.tb import tb
from pyatb.tools.smearing import gauss
from pyatb.fermi.orbital_selection import resolve_cohp_orbitals

import json
import os
import shutil
import numpy as np


def cohp_values_one_k(
    matrix,
    eigenvalues,
    eigenvectors,
    atom_i_orbs,
    atom_j_orbs,
    pair_factor=1.0,
):
    atom_i_orbs = np.asarray(atom_i_orbs, dtype=int)
    atom_j_orbs = np.asarray(atom_j_orbs, dtype=int)
    block = matrix[np.ix_(atom_i_orbs, atom_j_orbs)]
    coeff_i = eigenvectors[atom_i_orbs, :]
    coeff_j = eigenvectors[atom_j_orbs, :]
    values = pair_factor * np.einsum(
        "ib,ij,jb->b", coeff_i.conjugate(), block, coeff_j, optimize=True
    ).real
    return np.asarray(eigenvalues, dtype=float), values


class COHP:
    def __init__(
        self,
        tb: tb,
        **kwarg
    ):
        self.__tb = tb
        self.__max_kpoint_num = tb.max_kpoint_num
        if tb.nspin != 2:
            self.__tb_solver = (tb.tb_solver, )
        else:
            self.__tb_solver = (tb.tb_solver_up, tb.tb_solver_dn)
        self.__k_generator = None
        self.nspin = tb.nspin

        output_path = os.path.join(OUTPUT_PATH, "COHP")
        if RANK == 0:
            if os.path.exists(output_path):
                shutil.rmtree(output_path)
            os.mkdir(output_path)

        self.output_path = output_path

        if RANK == 0:
            with open(RUNNING_LOG, "a") as f:
                f.write("\n")
                f.write("\n------------------------------------------------------")
                f.write("\n|                                                    |")
                f.write("\n|                        COHP                        |")
                f.write("\n|                                                    |")
                f.write("\n------------------------------------------------------")
                f.write("\n\n")

    def set_k_generator(self, kpoint_mode=None, **kwarg):
        if kpoint_mode is None:
            kpoint_mode = "direct"
            kwarg["kpoint_direct_coor"] = np.array([[0.0, 0.0, 0.0]], dtype=float)

        if kpoint_mode == "mp":
            self.__k_generator = kpoint_generator.mp_generator(
                self.__max_kpoint_num,
                kwarg.get("k_start", np.array([0.0, 0.0, 0.0], dtype=float)),
                kwarg.get("k_vect1", np.array([1.0, 0.0, 0.0], dtype=float)),
                kwarg.get("k_vect2", np.array([0.0, 1.0, 0.0], dtype=float)),
                kwarg.get("k_vect3", np.array([0.0, 0.0, 1.0], dtype=float)),
                kwarg["mp_grid"],
            )
        elif kpoint_mode == "line":
            self.__k_generator = kpoint_generator.line_generator(
                self.__max_kpoint_num,
                kwarg["high_symmetry_kpoint"],
                kwarg["kpoint_num_in_line"],
            )
        elif kpoint_mode == "direct":
            self.__k_generator = kpoint_generator.array_generater(
                self.__max_kpoint_num,
                np.asarray(kwarg["kpoint_direct_coor"], dtype=float),
            )
        else:
            raise ValueError("unknown COHP kpoint_mode: %s" % kpoint_mode)

        if RANK == 0:
            with open(RUNNING_LOG, "a") as f:
                f.write("\nParameter setting of COHP kpoints : \n")
                f.write(" >> kpoint_mode : %s\n" % kpoint_mode)

    def _spin_indices(self, spin):
        spin = spin.lower()
        if self.nspin == 4:
            raise ValueError("COHP currently supports nspin = 1 or 2; non-collinear nspin = 4 is not implemented")
        if spin not in ("sum", "up", "down"):
            raise ValueError("COHP spin must be one of: sum, up, down")
        if self.nspin != 2:
            if spin in ("up", "down"):
                raise ValueError("COHP spin = up/down requires nspin = 2")
            return [0]
        if spin == "up":
            return [0]
        if spin == "down":
            return [1]
        return [0, 1]

    def _accumulate_one_k(self, spectrum, matrix, eigenvalues, eigenvectors, atom_i_orbs,
                          atom_j_orbs, energy_min, de, sigma, invert,
                          pair_factor=1.0, kpoint_weight=1.0, fermi_energy=None):
        e_num = spectrum.shape[0]
        energies, values = cohp_values_one_k(
            matrix,
            eigenvalues,
            eigenvectors,
            atom_i_orbs,
            atom_j_orbs,
            pair_factor=pair_factor,
        )
        # Integrate the discrete states before any windowing or broadening.
        # Reuse the contraction for both outputs; inversion affects display only.
        occupied = None
        if fermi_energy is not None:
            occupied = float(kpoint_weight * values[energies <= fermi_energy].sum())
        if invert:
            values = -values
        for energy, value in zip(energies, values):
            # Include boundary states and tails of states outside the plot window.
            # The continuous Gaussian is truncated at 10 sigma, without rescaling
            # the part that falls within the requested energy grid.
            start_index = max(0, int(np.ceil((energy - 10 * sigma - energy_min) / de)))
            end_index = min(e_num, int(np.floor((energy + 10 * sigma - energy_min) / de)) + 1)
            if start_index >= end_index:
                continue
            delta_E = energy_min + np.arange(start_index, end_index, dtype=float) * de - energy
            spectrum[start_index:end_index] += kpoint_weight * value * gauss(sigma, delta_E)
        return occupied

    @staticmethod
    def _pair_factor(atom_i_orbs, atom_j_orbs):
        first, second = set(atom_i_orbs), set(atom_j_orbs)
        if not first or not second:
            raise ValueError("COHP orbital selections must be non-empty")
        if first == second:
            return 1.0
        if first.isdisjoint(second):
            return 2.0
        raise ValueError("COHP orbital selections partially overlap; use identical or disjoint sets")

    @staticmethod
    def _normalized_kpoint_weights(total_kpoint_num, weights=None):
        if total_kpoint_num <= 0:
            raise ValueError("COHP requires a positive number of k points")
        if weights is None:
            return np.full(total_kpoint_num, 1.0 / total_kpoint_num, dtype=float)
        weights = np.asarray(weights, dtype=float)
        if weights.shape != (total_kpoint_num,):
            raise ValueError(
                "kpoint_weights must contain one value per generated k point "
                "(%d expected, %d received)" % (total_kpoint_num, weights.size)
            )
        if not np.all(np.isfinite(weights)) or np.any(weights < 0.0):
            raise ValueError("kpoint_weights must be finite and non-negative")
        scale = float(weights.max())
        if scale <= 0.0:
            raise ValueError("kpoint_weights must have a positive sum")
        weights = weights / scale
        return weights / weights.sum()

    def _matrix_for_method(self, solver, kpoints, method):
        method = method.upper()
        if method == "COHP":
            return solver.get_Hk(kpoints)
        if method == "COOP":
            return solver.get_Sk(kpoints)
        raise ValueError("COHP method must be COHP or COOP")

    def _spin_degeneracy_factor(self, spin):
        if self.nspin == 1 and spin.lower() == "sum":
            return 2.0
        return 1.0

    def calculate_cohp(self, fermi_energy, stru_file, atom_i_index, atom_j_index,
                       atom_i_orbs="all", atom_j_orbs="all", method="COHP", spin="sum",
                       e_range=None, de=0.05, sigma=0.15, invert=1, shift_to_efermi=1,
                       output_prefix="COHP", **kwarg):
        COMM.Barrier()
        timer.start("cohp", "calculate COHP")

        fermi_energy = float(fermi_energy)
        if not np.isfinite(fermi_energy):
            raise ValueError("fermi_energy must be finite")
        if e_range is None:
            e_range = [fermi_energy - 10.0, fermi_energy + 10.0]
        limits = np.asarray(e_range, dtype=float)
        if limits.shape != (2,) or not np.all(np.isfinite(limits)):
            raise ValueError("e_range must contain two finite energies")
        energy_min, energy_max = limits
        de, sigma = float(de), float(sigma)
        if not np.isfinite(de) or not np.isfinite(sigma) or de <= 0 or sigma <= 0:
            raise ValueError("de and sigma must be finite and positive")
        with np.errstate(over="ignore"):
            intervals = (energy_max - energy_min) / de
        if not np.isfinite(intervals) or intervals < 1:
            raise ValueError("e_range must increase and span at least one grid step de")
        e_num = int(np.floor(np.nextafter(intervals, np.inf))) + 1
        spectrum = np.zeros(e_num, dtype=float)
        occupied_integral = 0.0
        method = method.upper()
        if method not in ("COHP", "COOP"):
            raise ValueError("COHP method must be COHP or COOP")

        if self.__k_generator is None:
            self.set_k_generator(**kwarg)

        selection = resolve_cohp_orbitals(
            stru_file=stru_file,
            atom_i_index=atom_i_index,
            atom_j_index=atom_j_index,
            atom_i_orbs=atom_i_orbs,
            atom_j_orbs=atom_j_orbs,
            input_file=kwarg.get("input_file"),
            orbital_dir=kwarg.get("orbital_dir"),
        )
        spin_indices = self._spin_indices(spin)
        pair_factor = self._pair_factor(selection.atom_i_orbitals, selection.atom_j_orbitals)
        supplied_kpoint_weights = kwarg.get("kpoint_weights")
        kpoint_weights = self._normalized_kpoint_weights(
            self.__k_generator.total_kpoint_num,
            supplied_kpoint_weights,
        )
        self._last_pair_factor = pair_factor
        self._last_kpoint_weight_mode = "explicit_normalized" if supplied_kpoint_weights is not None else "uniform"

        if RANK == 0:
            with open(RUNNING_LOG, "a") as f:
                f.write("\nEnter the COHP calculation module ==> \n")
                f.write(" >> atom_i_index : %d\n" % atom_i_index)
                f.write(" >> atom_j_index : %d\n" % atom_j_index)
                f.write(" >> atom_i_orbs  : %s\n" % selection.atom_i_orbitals)
                f.write(" >> atom_j_orbs  : %s\n" % selection.atom_j_orbitals)
                f.write(" >> method       : %s\n" % method.upper())
                f.write(" >> spin         : %s\n" % spin)
                f.write(" >> e_range      : %.6f %.6f\n" % (energy_min, energy_max))
                f.write(" >> de           : %.6f\n" % de)
                f.write(" >> sigma        : %.6f\n" % sigma)

        chunk_start = 0
        for ik in self.__k_generator:
            chunk_weights = kpoint_weights[chunk_start:chunk_start + ik.shape[0]]
            chunk_start += ik.shape[0]
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, ik)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]
            if not kpoint_num:
                continue
            local_start = ik_process.ik_start_index
            local_weights = chunk_weights[local_start:local_start + kpoint_num]
            for ispin in spin_indices:
                solver = self.__tb_solver[ispin]
                eigenvectors, eigenvalues = solver.diago_H(ik_process.k_direct_coor_local)
                matrices = self._matrix_for_method(solver, ik_process.k_direct_coor_local, method)
                for single_k in range(kpoint_num):
                    occupied_integral += self._accumulate_one_k(
                        spectrum,
                        matrices[single_k],
                        eigenvalues[single_k],
                        eigenvectors[single_k],
                        selection.atom_i_orbitals,
                        selection.atom_j_orbitals,
                        energy_min,
                        de,
                        sigma,
                        bool(invert),
                        pair_factor=pair_factor,
                        kpoint_weight=local_weights[single_k],
                        fermi_energy=fermi_energy,
                    )

        spectrum = COMM.reduce(spectrum, root=0, op=op_sum)
        occupied_integral = COMM.reduce(occupied_integral, root=0, op=op_sum)

        if RANK == 0:
            degeneracy = self._spin_degeneracy_factor(spin)
            spectrum = spectrum * degeneracy
            occupied_integral *= degeneracy
            energy_grid = np.array([energy_min + i * de for i in range(e_num)], dtype=float)
            if shift_to_efermi:
                output_energy = energy_grid - fermi_energy
            else:
                output_energy = energy_grid
            self.print_data(
                output_prefix,
                output_energy,
                spectrum,
                selection,
                method,
                spin,
                fermi_energy,
                bool(shift_to_efermi),
                bool(invert),
                occupied_integral,
                {"sigma_eV": sigma, "de_eV": de,
                 "requested_energy_range_eV": limits.tolist(),
                 "sampled_energy_range_eV": [float(energy_grid[0]), float(energy_grid[-1])]},
            )

        COMM.Barrier()
        timer.end("cohp", "calculate COHP")

    def print_data(self, output_prefix, energy, spectrum, selection, method, spin,
                   fermi_energy, shifted, invert=False, occupied_integral=None, broadening=None):
        label = ("-" if invert else "") + method.upper()
        data_filename = "%s.dat" % output_prefix
        meta_filename = "%s.meta.json" % output_prefix
        np.savetxt(
            os.path.join(self.output_path, data_filename),
            np.c_[energy, spectrum],
            fmt="%0.8f",
            header="energy_eV %s" % label,
        )

        metadata = {
            "method": method.upper(),
            "spin": spin,
            "invert": bool(invert),
            "sign_convention": "negative population" if invert else "raw population",
            "spectrum_label": label,
            "spectrum_unit": "dimensionless" if method.upper() == "COHP" else "1/eV",
            "fermi_energy_eV": float(fermi_energy),
            "energy_shifted_to_efermi": shifted,
            "files": {
                "spectrum": data_filename,
            },
            "atom_i_orbitals": selection.atom_i_orbitals,
            "atom_j_orbitals": selection.atom_j_orbitals,
            "total_orbitals": selection.orbital_map.total_orbitals,
            "pair_convention": "unordered Hermitian orbital-set pair (identical: 1, disjoint: 2)",
            "pair_factor": float(getattr(self, "_last_pair_factor", 1.0)),
            "kpoint_weight_mode": getattr(self, "_last_kpoint_weight_mode", "uniform"),
            "atoms": [
                {
                    "index": atom.index,
                    "element": atom.element,
                    "orbital_file": atom.orbital_file,
                    "shell_counts": atom.shell_counts,
                    "start": atom.start,
                    "stop": atom.stop,
                }
                for atom in selection.orbital_map.atoms
            ],
        }
        if occupied_integral is not None:
            metadata["occupied_integral"] = {
                "method": "I" + method.upper(),
                "occupation": "E <= EF (zero temperature)",
                "unit": "eV" if method.upper() == "COHP" else "dimensionless",
                "raw_value": float(occupied_integral),
                "value": float(-occupied_integral if invert else occupied_integral),
                "label": ("-" if invert else "") + "I" + method.upper(),
                "spin_degeneracy": self._spin_degeneracy_factor(spin),
            }
        if broadening is not None:
            metadata["broadening"] = dict(broadening, kernel="exp(-(dE/sigma)^2)/(sqrt(pi)*sigma)",
                                          cutoff_sigma=10., window_renormalized=False)
        with open(os.path.join(self.output_path, meta_filename), "w") as file_obj:
            json.dump(metadata, file_obj, indent=4)

    def print_plot_script(self):
        script_path = os.path.join(self.output_path, "plot_cohp.py")
        with open(script_path, "w") as f:
            plot_script = """
import glob
import json
import os
import numpy as np
import matplotlib.pyplot as plt

work_path = os.getcwd()
meta_files = sorted(glob.glob(os.path.join(work_path, "*.meta.json")))
if not meta_files:
    raise FileNotFoundError("no COHP metadata file found")
with open(meta_files[0]) as file_obj:
    metadata = json.load(file_obj)

data = np.loadtxt(os.path.join(work_path, metadata["files"]["spectrum"]))
energy, cohp = data[:, 0], data[:, 1]

fig, ax = plt.subplots(1, 1, tight_layout=True)
ax.axvline(0.0, color="0.7", linewidth=0.8)
ax.axhline(0.0, color="0.7", linewidth=0.8)
ax.plot(cohp, energy)
ax.set_xlabel(metadata.get("spectrum_label", metadata["method"]))
ax.set_ylabel(r"$E-E_F$ (eV)" if metadata["energy_shifted_to_efermi"] else "Energy (eV)")
ax.set_title(metadata.get("spectrum_label", metadata["method"]))
plt.savefig(os.path.join(work_path, "cohp.pdf"))
plt.close("all")
"""
            f.write(plot_script)
