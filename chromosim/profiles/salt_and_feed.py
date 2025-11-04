from __future__ import annotations
import numpy as np

def salt_profile_clamped_linear(t, par):
    """
    Clamped linear ramp:
        I(t) = I_load + (I_elute - I_load) * clamp((t - grad_start)/(grad_end - grad_start), 0, 1)
    Accepts scalar or array t. Returns same shape as t.
    """
    tarr = np.atleast_1d(np.asarray(t, dtype=float))
    frac = (tarr - par.grad_start) / (par.grad_end - par.grad_start)
    frac = np.clip(frac, 0.0, 1.0)
    I = par.I_load + (par.I_elute - par.I_load) * frac
    return I if np.ndim(t) else float(I.item())

def feed_profile_by_load_volume(t, par):
    """
    Physical inlet concentration vector for Ns species (length Ns).
    During the 'load' window, Cin = Cfeed_total * frac; otherwise zeros.
    The load window is computed as:
        load_duration_s = 60 * load_mL / flow_mL_min
    If those fields are missing, defaults to 300 s.
    Works with scalar or array t. For scalar t returns 1×Ns array.
    """
    Ns = par.Nsp
    tarr = np.atleast_1d(np.asarray(t, dtype=float))
    Cin = np.zeros((tarr.size, Ns), dtype=float)

    # derive load duration if not given
    if not hasattr(par, "load_duration_s"):
        if hasattr(par, "load_mL") and hasattr(par, "flow_mL_min"):
            par.load_duration_s = 60.0 * float(par.load_mL) / float(par.flow_mL_min)
        else:
            par.load_duration_s = 300.0  # sensible default (5 min)

    mask_load = tarr < par.load_duration_s
    if np.any(mask_load):
        Cin[mask_load, :] = (par.Cfeed_total * np.asarray(par.frac)).reshape(1, Ns)

    if np.ndim(t) == 0:
        return Cin[0, :]
    return Cin
