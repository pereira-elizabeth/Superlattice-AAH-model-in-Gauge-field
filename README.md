# Superlattice AAH — Non-Hermitian Grid (B × v)

This repository contains a **non-Hermitian superlattice Aubry–André–Harper (AAH) model** implemented on a 1D multimodal chain.  
It demonstrates Hamiltonian construction, non-Hermitian diagonalization, and HPC-scale parameter sweeps.

---

## 🔍 What this project does
- **Builds the Hamiltonian from scratch**  
  - Superlattice chain with multimodal “islands” of width W = 3 or 4  
  - Complex Peierls phases (gauge field retained) on the hopping terms  
  - Imaginary onsite AAH potential (gain/loss) making the system **non-Hermitian**

- **Diagonalizes using biorthogonal eigenvectors**  
  - Solved via `scipy.linalg.eig(left=True, right=True)`  
  - Computes the **biorthogonal inverse participation ratio (IPR)** as a localization metric and computes the percentage of states in the p-manifold that have IPR above a limit.

- **Parallel computation on HPC (SLURM)**  
  - Sweeps a $(B_{micro}, v_{I})$ parameter grid  
  - Each job array element solves a subset of parameters in parallel  
  - Results are written to `results/imag_even_parallel_array_ALL_{job_id}.dat` with rows `(B, v, IPR%)`

---

## 📊 Example Output
Example: heatmap of the percentage of energy states in the p-manifold with IPR values across the $(B_{micro}, v_{I})$ grid that lie above a certain limit.

<p align="center">
  <img src="figure.png" alt="Example IPR heatmap" width="100%">
  <br><em>Figure 1 — Percentage of localized states in p-manifold across parameter space</em>
</p>

---
### What the figure shows
We reproduce **mobility-edge phase diagrams** for the full microscopic model over the \((B, v)\) grid.  
Color encodes the **percentage of localized states** (biorthogonal IPR threshold \(>\!1/L\)), not raw IPR values.  
Panels correspond to different island widths: **(a) \(W=3\)** and **(c) \(W=4\)**.  
As the gauge field \(B\) increases, orbital mixing shifts the **localization threshold** versus \(v\).

*Reference:* Pereira *et al.*, *Topology and criticality in non-Hermitian multimodal optical resonators through engineered losses*, arXiv:2509.05163. :contentReference[oaicite:1]{index=1}


---

## 🚀 How to Run

### Local (small grids)
```bash
pip install -r requirements.txt
python examples/superlattice_demo.py

