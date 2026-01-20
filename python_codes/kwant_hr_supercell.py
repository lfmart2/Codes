from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import h5py
import kwant
import numpy as np


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


def _build_hydrogen_coupling_map() -> dict[str, tuple[int, np.ndarray]]:
    mapping: dict[str, tuple[int, np.ndarray]] = {}
    layer_specs = {
        "first_d_layer": (range(2, 12, 2), range(3, 12, 2)),
        "second_d_layer": (range(12, 22, 2), range(13, 22, 2)),
        "third_d_layer": (range(22, 32, 2), range(23, 32, 2)),
        "first_sp_layer": (range(280, 288, 2), range(281, 288, 2)),
        "second_sp_layer": (range(288, 296, 2), range(289, 296, 2)),
        "third_sp_layer": (range(296, 304, 2), range(297, 304, 2)),
    }
    for prefix, (even_range, odd_range) in layer_specs.items():
        even_idx = np.fromiter(even_range, dtype=int)
        odd_idx = np.fromiter(odd_range, dtype=int)
        for spin_index, spin_label in enumerate(("up", "dn")):
            for target_label, target_idx in (("up", even_idx), ("dn", odd_idx)):
                for suffix in ("", "_im"):
                    label = f"{prefix}_{spin_label}_{target_label}{suffix}"
                    mapping[label] = (spin_index, target_idx)
    return mapping


HYDROGEN_COUPLING_MAP = _build_hydrogen_coupling_map()


def coupling_hydrogen_slab_mapping(label: str) -> tuple[int, np.ndarray] | None:
    """Map hydrogen coupling labels to integer indices."""
    return HYDROGEN_COUPLING_MAP.get(label)


@lru_cache(maxsize=32)
def _load_hydrogen_block_keys(
    file_path: str, r_vec: tuple[int, int, int]
) -> list[str]:
    with h5py.File(file_path, "r") as file:
        block_interaction = file[f"[{r_vec[0]},{r_vec[1]},{r_vec[2]}]"]
        return list(block_interaction.keys())


def coupling_hydrogen_slab(
    R_vec: tuple[int, int, int],
    num_slab: int,
    num_molecules: int,
    mol_distance: np.ndarray,
    file_path: str | Path,
) -> np.ndarray:
    """Compute hydrogen coupling block for a given lattice vector."""
    final_block = np.zeros(
        (num_slab + num_molecules, num_slab + num_molecules), dtype=complex
    )
    rows_up = np.arange(0, mol_distance.size, 2, dtype=int)
    rows_dn = rows_up + 1
    dist_up = mol_distance[0::2]
    dist_dn = mol_distance[1::2]
    list_Ri = _load_hydrogen_block_keys(str(file_path), R_vec)
    with h5py.File(file_path, "r") as file:
        block_interaction = file[f"[{R_vec[0]},{R_vec[1]},{R_vec[2]}]"]
        for label in list_Ri:
            if label.endswith("_im"):
                continue
            mapping = coupling_hydrogen_slab_mapping(label)
            if mapping is None:
                continue
            spin_index, interaction_indices = mapping
            data_re = block_interaction[label]
            data_im = block_interaction[f"{label}_im"]
            coeffs_re = np.vstack([data_re[:, 1], data_re[:, 0]])
            coeffs_im = np.vstack([data_im[:, 1], data_im[:, 0]])
            if spin_index == 0:
                distances = dist_up
                rows = rows_up
            else:
                distances = dist_dn
                rows = rows_dn
            values = np.polyval(coeffs_re, distances) + 1j * np.polyval(
                coeffs_im, distances
            )
            cols = interaction_indices + num_molecules
            final_block[np.ix_(rows, cols)] += values.T
    return final_block


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


def iter_hoppings_with_hydrogen(
    hr: HRData,
    num_molecules: int,
    mol_distance: np.ndarray,
    coupling_file: str | Path,
    tol_block: float | None = 1e-12,
) -> Iterable[tuple[tuple[int, int, int], np.ndarray]]:
    """Yield non-onsite hopping blocks with hydrogen coupling added."""
    block = np.zeros((hr.num_wann + num_molecules, hr.num_wann + num_molecules), dtype=complex)
    for idx in range(hr.nrpts):
        R_vec = tuple(int(v) for v in hr.R[:, idx])
        if R_vec == (0, 0, 0):
            continue
        block[:, :] = 0.0
        block[2 * num_molecules :, 2 * num_molecules :] = (
            hr.H_R[:, :, idx] / hr.weight[idx]
        )
        if tol_block is not None and np.max(np.abs(block)) < tol_block:
            continue
        hydrogen_block = coupling_hydrogen_slab(
            R_vec,
            num_slab=hr.num_wann,
            num_molecules=num_molecules,
            mol_distance=mol_distance,
            file_path=coupling_file,
        )
        yield R_vec, block + hydrogen_block


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


def build_supercell_with_hydrogen(
    hr: HRData,
    Lx: int,
    Ly: int,
    coupling_file: str | Path,
    poly_coeffs_up: list[float],
    poly_coeffs_dn: list[float],
    tol_block: float | None = 1e-12,
    seed: int | None = None,
) -> kwant.system.FiniteSystem:
    """Build a 2D supercell including hydrogen coupling data."""
    num_molecules = Lx * Ly // 2
    norb = hr.num_wann + num_molecules
    lat = kwant.lattice.square(norbs=norb)
    syst = kwant.Builder()

    rng = np.random.default_rng(seed)
    mol_distance = rng.uniform(-0.6, 0.6, size=2 * num_molecules)

    onsite_block = np.zeros((norb, norb), dtype=complex)
    for idx in range(hr.nrpts):
        if np.all(hr.R[:, idx] == 0):
            onsite_block[2 * num_molecules :, 2 * num_molecules :] = hr.H_R[
                :, :, idx
            ].copy()
            break

    mol_idx = 0
    for x in range(0, Lx, 2):
        for y in range(0, Ly, 2):
            onsite = onsite_block.copy()
            onsite[0, 0] = np.polyval(poly_coeffs_up, mol_distance[mol_idx])
            onsite[1, 1] = np.polyval(poly_coeffs_dn, mol_distance[mol_idx + 1])
            syst[lat(x, y)] = onsite
            mol_idx += 2

    for (Rx, Ry, Rz), hop_block in iter_hoppings_with_hydrogen(
        hr,
        num_molecules=num_molecules,
        mol_distance=mol_distance,
        coupling_file=coupling_file,
        tol_block=tol_block,
    ):
        if Rz != 0:
            continue
        for x in range(Lx):
            for y in range(Ly):
                x2 = (x + Rx) % Lx
                y2 = (y + Ry) % Ly
                s1 = lat(x, y)
                s2 = lat(x2, y2)
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
