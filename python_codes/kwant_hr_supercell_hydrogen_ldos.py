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
    """Load Wannier90 real-space Hamiltonian (wannier90_hr.dat)."""
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
        (num_slab + 2 * num_molecules, num_slab + 2 * num_molecules), dtype=complex
    )
    rows_up = np.arange(0, 2 * num_molecules, 2, dtype=int)
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
            cols = interaction_indices + 2 * num_molecules
            final_block[np.ix_(rows, cols)] += values.T
    return final_block


def iter_hoppings_with_hydrogen(
    hr: HRData,
    num_molecules: int,
    mol_distance: np.ndarray,
    coupling_file: str | Path,
    tol_block: float | None = 1e-12,
) -> Iterable[tuple[tuple[int, int, int], np.ndarray]]:
    """Yield non-onsite hopping blocks with hydrogen coupling added."""
    block = np.zeros(
        (hr.num_wann + 2 * num_molecules, hr.num_wann + 2 * num_molecules),
        dtype=complex,
    )
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


def build_supercell_with_hydrogen(
    hr: HRData,
    Lx: int,
    Ly: int,
    coupling_file: str | Path,
    poly_coeffs_up: list[float],
    poly_coeffs_dn: list[float],
    tol_block: float | None = 1e-12,
    seed: int | None = None,
) -> tuple[kwant.system.FiniteSystem, kwant.lattice.Monatomic, dict[tuple[int, int], int]]:
    """Build a 2D supercell including hydrogen coupling data (SOC)."""
    num_molecules = Lx * Ly // 2
    norb = hr.num_wann + 2 * num_molecules
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

    mol_sites: dict[tuple[int, int], int] = {}
    mol_idx = 0
    for x in range(0, Lx, 2):
        for y in range(0, Ly, 2):
            if mol_idx >= num_molecules:
                continue
            mol_sites[(x, y)] = mol_idx
            mol_idx += 1

    for x in range(Lx):
        for y in range(Ly):
            onsite = onsite_block.copy()
            mol_idx = mol_sites.get((x, y))
            if mol_idx is not None:
                onsite[2 * mol_idx, 2 * mol_idx] = np.polyval(
                    poly_coeffs_up, mol_distance[2 * mol_idx]
                )
                onsite[2 * mol_idx + 1, 2 * mol_idx + 1] = np.polyval(
                    poly_coeffs_dn, mol_distance[2 * mol_idx + 1]
                )
            syst[lat(x, y)] = onsite

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

    return syst.finalized(), lat, mol_sites


def hydrogen_orbital_map(
    lat: kwant.lattice.Monatomic, mol_sites: dict[tuple[int, int], int]
) -> dict[kwant.builder.Site, tuple[int, int]]:
    """Return hydrogen up/down orbital indices for each molecule site."""
    orbital_map: dict[kwant.builder.Site, tuple[int, int]] = {}
    for (x, y), mol_idx in mol_sites.items():
        orbital_map[lat(x, y)] = (2 * mol_idx, 2 * mol_idx + 1)
    return orbital_map


def compute_hydrogen_ldos(
    fsys: kwant.system.FiniteSystem,
    hydrogen_orbitals: dict[kwant.builder.Site, tuple[int, int]],
    energy_grid: np.ndarray,
    num_moments: int = 2000,
    rng_seed: int | None = 0,
) -> dict[kwant.builder.Site, dict[str, np.ndarray]]:
    """Compute LDOS for hydrogen up/down orbitals at each molecule site."""
    rng = np.random.default_rng(rng_seed)
    sites = list(hydrogen_orbitals.keys())
    density_operator = kwant.operator.Density(fsys, where=sites, sum=False)
    ldos = kwant.kpm.LocalDensityOfStates(
        fsys,
        num_moments=num_moments,
        rng=rng,
        operator=density_operator,
    )
    ldos_values = np.atleast_2d(ldos(energy_grid))
    norb = fsys.sites[0].family.norbs
    results: dict[kwant.builder.Site, dict[str, np.ndarray]] = {}

    if ldos_values.shape[1] == len(sites):
        for site_idx, site in enumerate(sites):
            results[site] = {
                "up": ldos_values[:, site_idx],
                "down": ldos_values[:, site_idx],
            }
        return results

    if ldos_values.shape[1] != len(sites) * norb:
        raise ValueError(
            "Unexpected LDOS output shape; expected per-site or per-orbital data."
        )

    for site_idx, site in enumerate(sites):
        start = site_idx * norb
        up_idx, dn_idx = hydrogen_orbitals[site]
        results[site] = {
            "up": ldos_values[:, start + up_idx],
            "down": ldos_values[:, start + dn_idx],
        }

    return results


def save_hydrogen_ldos_csv(
    ldos_by_site: dict[kwant.builder.Site, dict[str, np.ndarray]],
    energy_grid: np.ndarray,
    output_path: str | Path,
) -> None:
    """Save hydrogen LDOS (up/down) for each site to a CSV file."""
    import pandas as pd

    data: dict[str, np.ndarray] = {"Energies": np.asarray(energy_grid)}
    for site, spin_ldos in ldos_by_site.items():
        x, y = site.tag
        data[f"H_{x}_{y}_up"] = spin_ldos["up"]
        data[f"H_{x}_{y}_down"] = spin_ldos["down"]
    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)


def compute_hydrogen_pdos_kpm(
    fsys: kwant.system.FiniteSystem,
    mol_sites: dict[tuple[int, int], int],
    *,
    orbital_up: int = 0,
    orbital_down: int = 1,
    energy_grid: np.ndarray | None = None,
    num_moments: int = 1500,
    num_vectors: int = 30,
    rng_seed: int | None = 0,
) -> tuple[np.ndarray, list[tuple[int, int]], np.ndarray, np.ndarray]:
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

    site_to_id = {site: i for i, site in enumerate(fsys.sites)}
    dof_up: list[int] = []
    dof_down: list[int] = []
    for site in fsys.sites:
        xy = tuple(site.tag)
        if xy in coords_set:
            sid = site_to_id[site]
            dof_up.append(sid * norb + orbital_up)
            dof_down.append(sid * norb + orbital_down)

    dof_up = np.asarray(dof_up, dtype=int)
    dof_down = np.asarray(dof_down, dtype=int)

    if len(dof_up) != len(coords):
        missing = [xy for xy in coords if xy not in {tuple(s.tag) for s in fsys.sites}]
        raise ValueError(
            f"Requested {len(coords)} sites, matched {len(dof_up)}. Missing? {missing[:10]}"
        )

    def op(bra: np.ndarray, ket: np.ndarray, **_params: object) -> np.ndarray:
        v_up = np.conjugate(bra[dof_up]) * ket[dof_up]
        v_down = np.conjugate(bra[dof_down]) * ket[dof_down]
        return np.concatenate([v_up, v_down])

    rng = np.random.default_rng(rng_seed)
    spectrum = kwant.kpm.SpectralDensity(
        fsys,
        operator=op,
        num_moments=num_moments,
        num_vectors=num_vectors,
        rng=rng,
        mean=True,
    )
    if energy_grid is None:
        energies, dens = spectrum()
    else:
        energies, dens = spectrum(np.asarray(energy_grid))
    dens = np.asarray(dens)

    num_sites = len(coords)
    dens_2 = dens.reshape(len(energies), 2, num_sites)
    pdos_up = dens_2[:, 0, :]
    pdos_down = dens_2[:, 1, :]

    return energies, coords, pdos_up, pdos_down


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


if __name__ == "__main__":
    hr = load_hr("wannier90_hr_r0.dat")

    fsys, lat, mol_sites = build_supercell_with_hydrogen(
        hr,
        Lx=10,
        Ly=10,
        coupling_file="hydrogen_interaction_data.h5",
        poly_coeffs_up=[0.0, 0.0, 0.0],
        poly_coeffs_dn=[0.0, 0.0, 0.0],
        seed=1234,
    )
    report_system_size(fsys)

    energies = np.linspace(-2.0, 2.0, 400)
    hydrogen_orbitals = hydrogen_orbital_map(lat, mol_sites)
    ldos_by_site = compute_hydrogen_ldos(
        fsys,
        hydrogen_orbitals=hydrogen_orbitals,
        energy_grid=energies,
    )
    sample_site = next(iter(ldos_by_site))
    print("Sample hydrogen LDOS (up):", ldos_by_site[sample_site]["up"][:5])
    print("Sample hydrogen LDOS (down):", ldos_by_site[sample_site]["down"][:5])
