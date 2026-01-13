from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import kwant


@dataclass(slots=True)
class HRData:
    num_wann: int
    nrpts: int
    R: np.ndarray
    weight: np.ndarray
    H_R: np.ndarray


def load_hr(path: str | Path) -> HRData:
    """Load Wannier90 real-space Hamiltonian (wannier90_hr.dat).

    Returns
    -------
    HRData
        Parsed Wannier90 real-space Hamiltonian data. The arrays are:
        - R: integer lattice vectors with shape (3, nrpts)
        - weight: integer weights with shape (nrpts,)
        - H_R: complex hopping blocks with shape (nw, nw, nrpts)
    """
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        handle.readline()  # header
        num_wann = int(handle.readline())
        nrpts = int(handle.readline())

        weight: list[int] = []
        while len(weight) < nrpts:
            line = handle.readline()
            if not line:
                break
            line = line.strip()
            if line:
                weight.extend(int(x) for x in line.split())
        weight_arr = np.array(weight[:nrpts], dtype=int)

        R = np.zeros((3, nrpts), dtype=int)
        H_R = np.zeros((num_wann, num_wann, nrpts), dtype=np.complex128)

        R_index_map: dict[tuple[int, int, int], int] = {}
        idx_counter = 0

        for line in handle:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 7:
                continue

            Rx, Ry, Rz = map(int, parts[0:3])
            i, j = map(int, parts[3:5])
            real, imag = map(float, parts[5:7])

            key = (Rx, Ry, Rz)
            idx = R_index_map.get(key)
            if idx is None:
                if idx_counter >= nrpts:
                    raise ValueError(
                        "Parsed more R points than declared in header; "
                        "check the Wannier90 file."
                    )
                idx = idx_counter
                R_index_map[key] = idx
                R[:, idx] = [Rx, Ry, Rz]
                idx_counter += 1

            H_R[i - 1, j - 1, idx] = real + 1j * imag

        nr_eff = idx_counter
        return HRData(
            num_wann=num_wann,
            nrpts=nr_eff,
            R=R[:, :nr_eff],
            weight=weight_arr[:nr_eff],
            H_R=H_R[:, :, :nr_eff],
        )


def iter_hoppings(
    hr: HRData, tol_block: float | None = 1e-12
) -> Iterable[tuple[tuple[int, int, int], np.ndarray]]:
    """Yield non-onsite hopping blocks as (R, block)."""
    for idx in range(hr.nrpts):
        R_vec = tuple(int(v) for v in hr.R[:, idx])
        if R_vec == (0, 0, 0):
            continue

        block = hr.H_R[:, :, idx] / hr.weight[idx]
        if tol_block is not None and np.max(np.abs(block)) < tol_block:
            continue
        yield R_vec, block


def build_supercell_from_hr(
    hr: HRData,
    Lx: int,
    Ly: int,
    disorder_strength: float = 0.0,
    tol_block: float | None = 1e-12,
    seed: int | None = None,
) -> kwant.system.FiniteSystem:
    """Build a 2D supercell from Wannier90 HR data with PBC in x/y."""
    norb = hr.num_wann
    lat = kwant.lattice.square(norbs=norb)
    syst = kwant.Builder()

    onsite_block = np.zeros((norb, norb), dtype=complex)
    for idx in range(hr.nrpts):
        if np.all(hr.R[:, idx] == 0):
            onsite_block = hr.H_R[:, :, idx].copy()
            break

    rng = np.random.default_rng(seed)

    for x in range(Lx):
        for y in range(Ly):
            onsite = onsite_block.copy()
            if disorder_strength:
                shift = disorder_strength * (rng.random(norb) - 0.5)
                onsite += np.diag(shift)
            syst[lat(x, y)] = onsite

    for (Rx, Ry, Rz), hop_block in iter_hoppings(hr, tol_block=tol_block):
        if Rz != 0:
            continue
        for x in range(Lx):
            for y in range(Ly):
                x2 = (x + Rx) % Lx
                y2 = (y + Ry) % Ly
                s1 = lat(x, y)
                s2 = lat(x2, y2)
                if (s1, s2) in syst:
                    syst[s1, s2] = syst[s1, s2] + hop_block
                else:
                    syst[s1, s2] = hop_block

    return syst.finalized()


def report_system_size(fsys: kwant.system.FiniteSystem) -> None:
    H = fsys.hamiltonian_submatrix(sparse=True).tocsr()
    n = H.shape[0]
    nnz = H.nnz
    bytes_total = H.data.nbytes + H.indices.nbytes + H.indptr.nbytes

    print("------ Kwant system report ------")
    print(f"Number of Kwant sites: {len(fsys.sites)}")
    print(f"Orbitals per site: {fsys.sites[0].family.norbs}")
    print(f"Hilbert-space dimension N: {n}")
    print(f"Number of non-zeros (nnz): {nnz}")
    print(f"Sparse H memory: {bytes_total / 1024**2:.2f} MB")
    print(f"Dense H would be: {n * n * 16 / 1024**3:.2f} GB")
    print("---------------------------------")


def estimate_dos(
    fsys: kwant.system.FiniteSystem,
    energy_grid: np.ndarray,
    num_moments: int = 2000,
    rng_seed: int | None = 0,
) -> np.ndarray:
    """Compute DOS via KPM."""
    rng = np.random.default_rng(rng_seed)
    spectrum = kwant.kpm.SpectralDensity(
        fsys,
        num_moments=num_moments,
        rng=rng,
    )
    return spectrum(energy_grid)


if __name__ == "__main__":
    hr = load_hr("wannier90_hr_r0.dat")

    fsys = build_supercell_from_hr(
        hr,
        Lx=10,
        Ly=10,
        disorder_strength=0.0,
        tol_block=1e-3,
        seed=1234,
    )
    report_system_size(fsys)

    energies = np.linspace(-2.0, 2.0, 400)
    dos = estimate_dos(fsys, energies)
    print("DOS sample:", dos[:5])
