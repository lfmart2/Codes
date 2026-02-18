# Mathematica Codes (PhD Portfolio)

This folder is a **curated set of Julia** from my PhD and graduate research training.
It shows how I use Julia to move from theory to reproducible computational results in condensed matter physics. Here, I have added a plethora of post-processing tools that takes Tight-Binding modles, usually extracted from [Wannier90](https://wannier.org/), i.e., electronic structure calculations, and plots the corresponding outputs. Also, I have added the codes to calculate the electronic coupling and electron's wavefunction in a semiconducr Quantum Dot. For more codes and applications, please take a look in the [RibeiroGroup](https://github.com/RibeiroGroup) to explore all my codes and its applications.

Julia is a relatively new code language alternative capable of high-level coding for multiple purposes. I adventured to learn this new coding for its versatile and high efficiency. Take a [look](https://discourse.julialang.org/t/why-is-julia-so-great/94718/6) on this new language

---
 - **TB_Bandstructure:** Jupyter notebook that calculates the Band structure from the <seedname>_hr.dat (Wannier90) file. This includes the calculation of the collinear and non-collinear Spin. Notebook requieres the julia script Bandstructure.jl.
 - **TB_DOS:** Calculates the Density of States (DOS) from a full-diagonalized grid of the tight-binding model extacted from the files <seendname>_hr.dat from Wannier90. The Jupyter notebook requires the julia script TB_DOS.jl
 - **TB_FS:** Jupyter Notebook that calculates the 2D Fermi Surface of a Tight_Binding model and plots the contour lines arounf the Fermi Energy. Notebook runs with the <seedname>_hr.eig file and the Julia script FS_TB.jl
 - **TB-Greens.jl**: A very interesting script where calculates the surface greens function of a Tight-Binding model through principal layers.
 - **VAS_DOS**: An example of the post-processing calculations of the output files of VASP. This Notebook analyses the file PROCAR and calculates the corresponding DOS and PDOS. Notebook requires the Julia script VASP_DOS.jl
 - **VASP_Fermi_surface:** A Julia script that takes the VASP output file EIGENCAR and calculates the contour plots of the energy at fixed energies (Fermi surfaces).
 - **Electronic_couplings.jl:** Script that contains the codes to calculates the electornic couplings and wavefunctions of a semiconductor Quantum Dot. Here, I have implemented such interaction with a tunable QD size and at different Core/Shells radius.
 - **Fano_IR_1S.jl:** Applications of Julia codes to the exploration of assymetry lineshape lines in a pump-probe scattering scheme on semiconductor Quantum Dots.