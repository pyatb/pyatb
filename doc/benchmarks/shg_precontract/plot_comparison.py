"""Compare full raw SHG tensors before plotting representative GaAs spectra.

The overlay establishes equivalence to the original method=0 API, not absolute
physical accuracy, material convergence, or a validated SI susceptibility.
"""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def read_result(path):
    with np.load(path, allow_pickle=False) as result:
        return {key: result[key] for key in ['omega_eV', 'kpoints', 'chi_raw']}


def component_name(index):
    return ''.join('xyz'[i] for i in (index // 9, index // 3 % 3, index % 3))


def plot_comparison(original, candidate, output):
    original, candidate, output = Path(original), Path(candidate), Path(output)
    a, b = read_result(original), read_result(candidate)
    for key in ['omega_eV', 'kpoints']:
        np.testing.assert_array_equal(a[key], b[key])
    x, y, omega = a['chi_raw'], b['chi_raw'], a['omega_eV']
    if x.shape != (27, len(omega)) or y.shape != x.shape:
        raise ValueError('Expected raw C++ tensor shape (27, number of energies)')
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Nonfinite response')
    np.testing.assert_allclose(y, x, rtol=1e-10, atol=1e-12)
    peak = np.maximum(np.max(abs(x), axis=1), 1e-30)
    error = np.max(abs(y-x), axis=1) / peak
    dominant = int(np.argmax(peak))
    # Display xyz plus the largest raw component; selection uses only baseline
    # magnitudes, never which component gives the smallest comparison error.
    selected = [5, dominant] if dominant != 5 else [5, 0]
    fig, axes = plt.subplots(3, 2, figsize=(7.3, 7.8), constrained_layout=True)
    for row, component in enumerate(selected):
        for col, transform in enumerate([np.real, np.imag]):
            ax = axes[row, col]
            ax.plot(omega, transform(x[component]), color='#243b53', lw=1.6, label='Original')
            stride = max(1, len(omega)//35)
            ax.plot(omega[::stride], transform(y[component])[::stride], 'o',
                    color='#d46a32', mfc='none', ms=3.8, mew=.9, label='Precontracted')
            ax.set_title(('Real' if col == 0 else 'Imaginary') + ' part: ' + component_name(component))
            ax.set_ylabel('Raw response (no SI prefactor)')
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
            if row == 0: ax.legend(frameon=False, fontsize=9)
    axes[2, 0].semilogy(omega, np.maximum(np.max(abs(y-x), axis=0)/np.max(abs(x)), 1e-18),
                         color='#243b53', lw=1.)
    axes[2, 0].set_ylabel('Global-peak-normalized error')
    axes[2, 1].semilogy(np.arange(27), np.maximum(error, 1e-18), 'o', color='#243b53', ms=3.5)
    axes[2, 1].set_xlabel('Cartesian component index')
    axes[2, 1].set_ylabel('Component-peak-normalized error')
    for i, ax in enumerate(axes.flat):
        if i != 5: ax.set_xlabel('Photon energy (eV)')
        ax.text(-.15, 1.05, '('+chr(97+i)+')', transform=ax.transAxes, weight='bold')
        ax.grid(alpha=.15)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output)+'.pdf', bbox_inches='tight')
    fig.savefig(str(output)+'.png', dpi=180, bbox_inches='tight')
    plt.close(fig)
    summary = dict(nk=len(a['kpoints']), nomega=len(omega), rtol=1e-10, atol=1e-12,
                   full_complex_allclose=True, displayed_components=[component_name(i) for i in selected],
                   max_abs=float(np.max(abs(y-x))), max_component_peak_relative=float(np.max(error)),
                   component_peak_relative=error.tolist(),
                   original_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),
                   candidate_sha256=hashlib.sha256(candidate.read_bytes()).hexdigest(),
                   scope='Equivalence to the original raw API, not absolute physical validation')
    Path(str(output)+'.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='output prefix, without extension')
    args = parser.parse_args()
    print(json.dumps(plot_comparison(args.original, args.candidate, args.output), indent=2))


if __name__ == '__main__':
    main()
