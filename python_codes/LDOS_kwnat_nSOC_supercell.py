from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import h5py
import kwant
import numpy as np
import sys
import pandas as pd

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

def _build_hydrogen_coupling_map() -> dict[str, tuple[int, int, np.ndarray]]:
    mapping: dict[str, tuple[int, int, np.ndarray]] = {}
    layer_specs = {
        "first_d_layer":  range(0,   5),
        "second_d_layer": range(5,  10),
        "third_d_layer":  range(10, 15), 
        "first_sp_layer": range(140, 144),
        "second_sp_layer": range(144, 148),
        "third_sp_layer": range(148, 152),
    }
    for prefix, (interaction_range) in layer_specs.items():
        interaction_idx = np.fromiter(interaction_range, dtype=int)
        for suffix in ("", "_im"):
            label = f"{prefix}{suffix}"
            mapping[label] = (interaction_idx)
    return mapping


HYDROGEN_COUPLING_MAP = _build_hydrogen_coupling_map()


def coupling_hydrogen_slab_mapping(
    label: str,
) -> tuple[int, int, np.ndarray] | None:
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
    mol_distance: np.ndarray,
    file_path: str | Path,
    x: int,
    y: int,
    Lx: int,
    Ly: int,
    mol_sites: dict[tuple[int, int], int],
) -> np.ndarray:
    """Compute hydrogen coupling block for a given lattice vector."""
    hydrogen_orbitals = 1
    final_block = np.zeros(
        (num_slab + hydrogen_orbitals, num_slab + hydrogen_orbitals), dtype=complex
    )
    rows = np.array([0], dtype=int)
    src_idx = mol_sites.get((x, y))
    dst_idx = mol_sites.get(((x + R_vec[0]) % Lx, (y + R_vec[1]) % Ly))
    has_src = src_idx is not None
    has_dst = dst_idx is not None
    if not has_src and not has_dst:
        return final_block
    distances = (
        np.array([mol_distance[src_idx]], dtype=float) if has_src else np.array([])
    )
    distances_conj = (
        np.array([mol_distance[dst_idx]], dtype=float) if has_dst else np.array([])
    )

    list_Ri = _load_hydrogen_block_keys(str(file_path), R_vec)
    with h5py.File(file_path, "r") as file:
        block_interaction = file[f"[{R_vec[0]},{R_vec[1]},{R_vec[2]}]"]
        for label in list_Ri:
            if label.endswith("_im"):
                continue
            mapping = coupling_hydrogen_slab_mapping(label)
            if mapping is None:
                continue
            interaction_indices = mapping
            data_re = block_interaction[label]
            data_im = block_interaction[f"{label}_im"]
            cols = interaction_indices + hydrogen_orbitals
            if has_src:
                values = (
                    data_re[:, 3][None, :] * distances[:, None] ** 3 +
                    data_re[:, 2][None, :] * distances[:, None] ** 2 +
                    data_re[:, 1][None, :] * distances[:, None] + 
                    data_re[:, 0][None, :]
                ) + 1j * (
                    data_im[:, 3][None, :] * distances[:, None] ** 3 +
                    data_im[:, 2][None, :] * distances[:, None] ** 2 +
                    data_im[:, 1][None, :] * distances[:, None] + 
                    data_im[:, 0][None, :]
                )
                final_block[np.ix_(rows, cols)] += values
            if has_dst:
                values_conj = (
                    data_re[:, 3][None, :] * distances_conj[:, None] ** 3 +
                    data_re[:, 2][None, :] * distances_conj[:, None] ** 2 +
                    data_re[:, 1][None, :] * distances_conj[:, None] + 
                    data_re[:, 0][None, :]
                ) + 1j * (
                    data_im[:, 3][None, :] * distances_conj[:, None] ** 3 +
                    data_im[:, 2][None, :] * distances_conj[:, None] ** 2 +
                    data_im[:, 1][None, :] * distances_conj[:, None] + 
                    data_im[:, 0][None, :]
                )
                final_block[np.ix_(cols, rows)] += values_conj.T
    return final_block


def _rvec_is_positive(r_vec: tuple[int, int, int]) -> bool:
    Rx, Ry, Rz = r_vec
    return (Rz > 0) or (Rz == 0 and Ry > 0) or (Rz == 0 and Ry == 0 and Rx > 0)

def iter_hoppings_with_hydrogen(
    hr: HRData,
    mol_distance: np.ndarray,
    coupling_file: str | Path,
    x: int,
    y: int,
    Lx: int,
    Ly: int,
    mol_sites: dict[tuple[int, int], int],
) -> Iterable[tuple[tuple[int, int, int], np.ndarray]]:
    """Yield non-onsite hopping blocks with hydrogen coupling added."""
    hydrogen_orbitals = 1
    block = np.zeros(
        (hr.num_wann + hydrogen_orbitals, hr.num_wann + hydrogen_orbitals),
        dtype=complex,
    )
    for idx in range(hr.nrpts):
        R_vec = tuple(int(v) for v in hr.R[:, idx])
        if R_vec == (0, 0, 0):
            continue
        if not _rvec_is_positive(R_vec):
            continue
        if R_vec[2] != 0:
            continue
        block[:, :] = 0.0
        block[hydrogen_orbitals:, hydrogen_orbitals:] = (
            hr.H_R[:, :, idx] / hr.weight[idx]
        )
        dst_site = ((x + R_vec[0]) % Lx, (y + R_vec[1]) % Ly)
        if (x, y) not in mol_sites and dst_site not in mol_sites:
            yield R_vec, block
            continue
        hydrogen_block = coupling_hydrogen_slab(
            R_vec,
            num_slab=hr.num_wann,
            mol_distance=mol_distance,
            file_path=coupling_file,
            x=x,
            y=y,
            Lx=Lx,
            Ly=Ly,
            mol_sites=mol_sites,
        )
        yield R_vec, block + hydrogen_block


def build_supercell_with_hydrogen(
    hr: HRData,
    Lx: int,
    Ly: int,
    coupling_file: str | Path,
    poly_coeffs: list[float],
) -> kwant.system.FiniteSystem:
    """Build a 2D supercell including hydrogen coupling data."""
    num_molecules = (Lx * Ly) // 2
    hydrogen_orbitals = 1
    norb = hr.num_wann + hydrogen_orbitals
    lat = kwant.lattice.square(norbs=norb)
    syst = kwant.Builder()

    mol_distance = np.random.uniform(
        -0.6, 0.6, size=hydrogen_orbitals * num_molecules
    )

    onsite_block = np.zeros((norb, norb), dtype=complex)
    for idx in range(hr.nrpts):
        if np.all(hr.R[:, idx] == 0):
            onsite_block[hydrogen_orbitals:, hydrogen_orbitals:] = hr.H_R[
                :, :, idx
            ].copy()
            break

    mol_sites: dict[tuple[int, int], int] = {}
    mol_idx = 0
    for x in range(Lx):
        for y in range(Ly):
            if (x + y) % 2 != 0:
                continue
            if mol_idx >= mol_distance.size:
                continue
            mol_sites[(x, y)] = mol_idx
            mol_idx += hydrogen_orbitals

    for x in range(Lx):
        for y in range(Ly):
            onsite = onsite_block.copy()
            mol_idx = mol_sites.get((x, y))
            if mol_idx is not None:
                onsite[0, 0] = np.polyval(poly_coeffs, mol_distance[mol_idx])
            syst[lat(x, y)] = onsite

    for x in range(Lx):
        for y in range(Ly):
            for (Rx, Ry, Rz), hop_block in iter_hoppings_with_hydrogen(
                hr,
                mol_distance=mol_distance,
                coupling_file=coupling_file,
                x=x,
                y=y,
                Lx=Lx,
                Ly=Ly,
                mol_sites=mol_sites,
            ):
                x2 = (x + Rx) % Lx
                y2 = (y + Ry) % Ly
                s1 = lat(x, y)
                s2 = lat(x2, y2)
                syst[s1, s2] = hop_block

    return syst.finalized(), mol_sites

def compute_hydrogen_pdos_kpm(
    fsys: kwant.system.FiniteSystem,
    mol_sites: Dict[Tuple[int, int], int],
    *,
    orbital_H: int = 0,
    energy_grid: np.ndarray | None = None,
    num_moments: int = 1500,
    num_vectors: int = 30,
) -> tuple[np.ndarray, list[tuple[int, int]], np.ndarray]:
    """Compute hydrogen PDOS for two orbitals using KPM + custom operator.

    Returns
    -------
    energies : (NE,)
    coords : list of (x, y) in mol_sites order
    pdos_up : (NE, Ns)
    pdos_down : (NE, Ns)
    """
    norb = fsys.sites[0].family.norbs
    coords = [xy for xy, _ in sorted(mol_sites.items(), key=lambda kv: kv[1])]
    coords_set = set(coords)
    Ns = len(coords)

    site_to_id = {site: i for i, site in enumerate(fsys.sites)}
    dof_H: list[int] = []
    for site in fsys.sites:
        xy = tuple(site.tag)
        if xy in coords_set:
            sid = site_to_id[site]
            dof_H.append(sid * norb + orbital_H)

    dof_H = np.asarray(dof_H, dtype=int)

    if len(dof_H) != len(coords):
        missing = [xy for xy in coords if xy not in {tuple(s.tag) for s in fsys.sites}]
        raise ValueError(
            f"Requested {len(coords)} sites, matched {len(dof_H)}. Missing? {missing[:10]}"
        )

    def op(bra: np.ndarray, ket: np.ndarray, **_params: object) -> np.ndarray:
        v_H = np.conjugate(bra[dof_H]) * ket[dof_H]
        return v_H

    spectrum = kwant.kpm.SpectralDensity(
        fsys,
        operator=op,
        num_moments=num_moments,
        num_vectors=num_vectors,
        mean=False,
    )
    if energy_grid is None:
        result = spectrum()
    else:
        result = spectrum(np.asarray(energy_grid))
    if len(result) == 2:
        energies, dens = result
    elif len(result) == 3:
        energies, dens, _errors = result
    else:
        raise RuntimeError(f"Unexpected spectrum return values: {len(result)}")
    dens = np.asarray(dens)

    # Robust shape check (do not reshape silently)
    if dens.ndim != 2 or dens.shape[1] != Ns:
        raise RuntimeError(f"Unexpected dens shape {dens.shape}; expected (NE, {Ns}).")

    return energies, coords, dens

if __name__ == "__main__":
    if len(sys.argv) < 5:
        raise SystemExit("Usage: python DOS_kwant.py <wannier90_hr_up.dat> <output_up.csv> <wannier90_hr_dn.dat> <output_dn.csv>")
    hr = load_hr(sys.argv[1])

    E_F = 5.58170362

    coupling_file_up = "nSOC_linregress_up.h5"
    coupling_file_dn = "nSOC_linregress_dn.h5"
    
    poly_coeffs_up = [
        -0.6888020535087722,
        -0.11799012444444447 + E_F,
    ]
    poly_coeffs_dn = [
        -0.7410776675438597,
        -0.1485784577777777 + E_F,
    ]

    fsys, mol_sites = build_supercell_with_hydrogen(
        hr,
        Lx=10,
        Ly=10,
        coupling_file=coupling_file_up,
        poly_coeffs=poly_coeffs_up,
    )

    filename_up = sys.argv[3]
    emin_up = -1.28
    emax = E_F + 1
    energies_up = np.linspace(emin_up, emax, 800)

    energies, coords, pdos_up = compute_hydrogen_pdos_kpm(
        fsys,
        mol_sites,
        energy_grid=energies_up,
        num_moments=1500,
        num_vectors=30,
    )

    df = pd.DataFrame({"Energies": energies})
    for idx, (x, y) in enumerate(coords):
        df[f"H_{x}_{y}_up"] = pdos_up[:, idx]
    df.to_csv(filename_up, index=False)

    # Spin-down calculation
    hr = load_hr(sys.argv[2])
    fsys, mol_sites = build_supercell_with_hydrogen(
        hr,
        Lx=10,
        Ly=10,
        coupling_file=coupling_file_dn,
        poly_coeffs=poly_coeffs_dn,
    )

    filename_dn = sys.argv[4]
    emin_dn = -1.16
    energies_dn = np.linspace(emin_dn, emax, 800)
    dos_dn = rho_dn(energies_dn)
    
    energies, coords, pdos_dn = compute_hydrogen_pdos_kpm(
        fsys,
        mol_sites,
        energy_grid=energies_dn,
        num_moments=1500,
        num_vectors=30,
    )
    df = pd.DataFrame({"Energies": energies})
    for idx, (x, y) in enumerate(coords):
        df[f"H_{x}_{y}_dn"] = pdos_dn[:, idx]
    df.to_csv(filename_dn, index=False)
