from __future__ import annotations
from types import SimpleNamespace as NS
from ..profiles.salt_and_feed import salt_profile_clamped_linear, feed_profile_by_load_volume

def params_ev_aex():
    """
    Multi-species parameter set matching your MATLAB 'params_ev_aex.m'.
    Returns a SimpleNamespace so attributes can be attached.
    """
    par = NS()
    # Geometry / grid
    par.L    = 0.025         # m
    par.dID  = 0.007         # m
    par.eps  = 0.40
    par.N    = 80
    par.dz   = par.L / par.N

    # Hydro
    par.u    = 1.0e-4        # m s^-1
    par.Dax  = 8e-11         # m^2 s^-1

    # EV classes
    par.zeta_mV = [-12, -25, -45]
    par.frac    = [0.10, 0.15, 0.75]
    par.Nsp     = len(par.zeta_mV)

    # Scaling / capacity
    par.Cfeed_total = 7.5e9        # particles m^-3
    par.qmax        = 4e13         # particles m^-3 wall
    par.Cscale      = par.Cfeed_total
    par.Qscale      = par.qmax
    par.CplotScale  = 1e10
    par.CplotLabel  = r"10^{10} per m^{-3}"

    # Kinetics & zeta mapping
    par.kd0   = 2e-3
    par.K0    = 1e-22
    par.gamma = 0.53
    par.Imax  = 1.2
    par.fI    = lambda I: max(0.0, 1.0 - float(I)/par.Imax)

    # Gradient
    par.grad_start = 460.0
    par.grad_end   = par.grad_start + 1200.0
    par.I_load     = 0.05
    par.I_elute    = 1.0
    par.salt_profile = lambda t: salt_profile_clamped_linear(t, par)

    # Feed (physical units)
    # flow & load (optional but enables volume-based load windows)
    par.flow_mL_min = 1.0
    par.load_mL     = 5.0
    par.feed_profile = lambda t: feed_profile_by_load_volume(t, par)

    # Time / solver tolerances (used by runner)
    par.tspan = (0.0, 3000.0)       # s
    par.rtol  = 1e-6
    par.atol  = 1e-6
    par.max_step = 8.0

    # Plot helper
    par.B_start = 10.0
    par.B_end   = 100.0
    return par
