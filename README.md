# ChromoSim
1D ion-exchange chromatography column simulator (mechanistic; multi-species; gradient; FPLC-style outputs)

**ChromoSim** is a general, mechanistic simulator for ion-exchange chromatography
(anionic or cationic). It models 1D convection–dispersion with competitive binding
and supports user-defined affinity maps (e.g., surface-potential dependent) or 
constant-affinity protein models.

## Highlights
- Multi-species breakthrough & in-column heatmaps
- FPLC-style output (mL / pseudo-mAU) with %B overlay
- Fast stiff-ODE integration with sparse Jacobian
- Works for proteins, nanoparticles, EVs (by configuration)

## Repo layout
- `src/` library code (added later)
- `examples/` runnable demos (added later)
- `assets/` figures/schematics (MIT-licensed)
- `data/` small demo data (LFS)
- `docs/` user docs & images
- `tests/` lightweight checks

## Quick start
Examples and core library will be added in v0.1.0.  
For now, see `docs/` for the modeling overview.

## License
MIT for code and assets.
