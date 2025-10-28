# ChromoSim Documentation

ChromoSim is a 1D ion-exchange column simulator with:
- convection–dispersion transport
- competitive binding (Langmuir-style)
- configurable affinity models (constant-K or surface-potential mapped)
- gradient programming and FPLC-style outputs

## Getting started
- Install MATLAB (R2022b+) or Octave (experimental)
- Clone the repo, see `examples/` (added in v0.1.0)

## Concepts
- State variables: mobile concentration c(z,t), bound-phase q(z,t)
- Parameters: flow, dispersion, porosity, capacity, gradient
- Outputs: breakthrough, heatmaps, FPLC trace with %B overlay
