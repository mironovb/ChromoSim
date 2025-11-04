# ChromoSim

Open-Source chromatography modeling software - 1D ion-exchange chromatography column simulator (mechanistic model; multi-species; with gradient; FPLC-style outputs)

**ChromoSim** simulates axial convection–dispersion with competitive Langmuir binding. Affinity can be constant (e.g., proteins) or mapped to a per-species property (e.g., surface potential) with salt-dependent screening. Outputs: breakthrough curves, in-column heatmaps, and an FPLC-style trace (mL / pseudo-mAU) with a %B overlay.

## Highlights

* Multi-species breakthrough + in-column heatmaps
* FPLC-style plot (elution volume in mL, pseudo-mAU) with %B overlay
* Fast stiff-ODE solve with a sparse Jacobian
* Configurable for proteins, nanoparticles, or EVs

## Install

```bash
# from the repo root
python -m venv .venv
source .venv/bin/activate            # Windows: .venv\Scripts\activate
pip install -U pip setuptools wheel
pip install -r requirements.txt      # installs ChromoSim in editable mode
```

## Run a demo

```bash
# CLI entry point (if available)
chromosim-run

# or run the example script directly
python -m scripts.main_ev_aex
```

This builds derivative matrices, integrates the model, and opens:

* Breakthrough/elution with %B overlay
* In-column concentration heatmap
* FPLC-style trace (mL / pseudo-mAU) with %B overlay

## Configure a run

Edit parameters in **`scripts/main_ev_aex.py`** (or create your own script):

* **Column & grid:** `L`, `N` (with `dz = L/N`), `eps`
* **Hydrodynamics:** `u` (via flow), `Dax`
* **Species:** `Nsp`, `frac`, optional per-species property array
* **Capacity:** `qmax`
* **Kinetics:** `kd0`, `K0`, `gamma`, and `fI(I)` (or set `gamma=0` for constant affinity)
* **Programs:** `grad_start`, `grad_end`, `I_load`, `I_elute` and `salt_profile(t)`
* **Feed:** `Cfeed_total`, `frac`, `feed_profile(t)` (e.g., 5 min load at 1 mL/min)

## What the model solves

# Chromatography Column Model

## Core Equations

### Mobile Phase Mass Balance
**For species i:**
$$\partial_t C_i = -u\,\partial_z C_i + D_{\mathrm{ax}}\,\partial_{zz} C_i - \frac{1-\varepsilon}{\varepsilon}\,\partial_t Q_i$$

### Stationary Phase Kinetics (Competitive Langmuir)
**For species i:**
$$\partial_t Q_i = k_{a,i}\,C_i\left(Q_{\max}-\sum_j Q_j\right) - k_{d,i}\,Q_i$$

### Optional Salt-Dependent Affinity
**For species i:**
$$K_i = K^\circ \exp\left(-\frac{\gamma\,|P_i|}{I_\*}(I - I_{\mathrm{ref}})\right)$$

Use constant $K_i$ by setting $\gamma=0$.

## Key Parameters
- **Mobile phase**: $C_i$ (concentration), $u$ (velocity), $D_{\mathrm{ax}}$ (dispersion), $\varepsilon$ (porosity)
- **Stationary phase**: $Q_i$ (loading), $Q_{\max}$ (capacity), $k_{a,i}$, $k_{d,i}$ (rate constants)
- **Salt effects**: $I$ (ionic strength), $I_{\mathrm{ref}}$ (reference), $\gamma$ (modulation)

[See full documentation for detailed parameter explanations and theory.]

## Repo layout

* `chromosim/` core library

  * `numerics/` derivative matrices, sparsity, ODE system
  * `models/` parameter helpers (profiles, defaults)
* `scripts/` runnable examples
* `docs/` user docs & public figures (MIT)
* `requirements.txt`, `pyproject.toml` packaging

## Tips

* Coarser grids (`N`) run faster; refine only to resolve peak shape.
* If plots don’t show on macOS, try: `export MPLBACKEND=TkAgg` before running.
* To compare with instrument data, adjust `eps_det` in the example to match the mAU scale visually.

## Contributing

Issues and PRs welcome. See `CONTRIBUTING.md` and `CODE_OF_CONDUCT.md`.

## License & citation

MIT for code and assets. If you use ChromoSim in a publication, please cite the repo (`CITATION.cff`).
