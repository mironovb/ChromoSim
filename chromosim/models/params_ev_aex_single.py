from __future__ import annotations
from types import SimpleNamespace as NS
from ..profiles.salt_and_feed import salt_profile_clamped_linear, feed_profile_by_load_volume

def params_ev_aex_single():
    """Single-species parameter set matching your MATLAB 'params_ev_aex_single.m'."""
    par = NS()
    # Geometry / grid
    par.L    = 0.025
    par.dID  = 0.007
    par.eps  = 0.40
    par.N    = 80
    par.dz   = par.L / par.N

    # Hydro
    par.u    = 1.0e-4
    par.Dax  = 9e-11

    # EV class (single)
    par.zeta_mV = [-25]
    par.frac    = [1.0]
    par.Nsp     = 1

    # Scaling / capacity
    par.Cfeed_total = 5e9
    par.qmax        = 4e14
    par.Cscale      = par.Cfeed_total
    par.Qscale      = par.qmax
    par.CplotScale  = 1e10
    par.CplotLabel  = r"10^{10} per m^{-3}"

    # Kinetics & zeta mapping
    par.kd0   = 1e-2
    par.K0    = 1e-22
    par.gamma = 0.750
    par.Imax  = 1.2
    par.fI    = lambda I: max(0.0, 1.0 - float(I)/par.Imax)

    # Gradient
    par.grad_start = 12.0*60.0
    par.grad_end   = par.grad_start + 20.0*60.0
    par.I_load     = 0.01
    par.I_elute    = 1.0
    par.salt_profile = lambda t: salt_profile_clamped_linear(t, par)

    # Feed (physical units)
    # flow & load for a clear 5 mL at 1 mL/min → 300 s load
    par.flow_mL_min = 1.0
    par.load_mL     = 5.0
    par.feed_profile = lambda t: feed_profile_by_load_volume(t, par)

    # Time / solver tolerances
    par.tspan = (0.0, 45.0*60.0)
    par.rtol  = 1e-6
    par.atol  = 1e-6
    par.max_step = 8.0

    # Plot helper
    par.B_start = 10.0
    par.B_end   = 100.0
    return par
