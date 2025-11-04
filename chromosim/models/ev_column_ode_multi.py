from __future__ import annotations
import numpy as np

def _get(par, key):
    try:
        return getattr(par, key)
    except AttributeError:
        return par[key]

def ev_column_ode_multi(t: float, y: np.ndarray, par) -> np.ndarray:
    """
    Multi-species 1D convection–dispersion + competitive Langmuir adsorption
    with zeta- and ionic-strength–dependent kinetics, integrated in *scaled* vars.

    y = [ c(:) ; qhat(:) ], with c,Q scaled by Cscale,Qscale respectively.
    """
    N   = int(_get(par, "N"))
    Ns  = int(_get(par, "Nsp"))
    dz  = float(_get(par, "dz"))
    u   = float(_get(par, "u"))
    Dax = float(_get(par, "Dax"))
    eps = float(_get(par, "eps"))

    Cscale = float(_get(par, "Cscale"))
    Qscale = float(_get(par, "Qscale"))
    qmax   = float(_get(par, "qmax"))

    zeta_mV = np.asarray(_get(par, "zeta_mV"), dtype=float).reshape(-1)  # (Ns,)
    NN = N * Ns

    c    = y[0:NN].reshape(N, Ns)            # dimensionless
    qhat = y[NN:2*NN].reshape(N, Ns)         # dimensionless

    C = c    * Cscale                         # physical
    Q = qhat * Qscale

    I  = float(_get(par, "salt_profile")(t))  # ionic strength (scalar)
    fI = float(_get(par, "fI")(I))

    K0    = float(_get(par, "K0"))
    gamma = float(_get(par, "gamma"))
    Ki_vec = K0 * np.exp(gamma * np.abs(zeta_mV) * fI)  # (Ns,)

    # desorption (constant or function of I)
    kd_fun = getattr(par, "kd_fun", None) if hasattr(par, "kd_fun") else _get(par, "kd0")
    if callable(kd_fun):
        kd_vec = np.asarray(kd_fun(I), dtype=float).reshape(-1)
        if kd_vec.size == 1:
            kd_vec = np.repeat(kd_vec, Ns)
    else:
        kd_vec = np.ones(Ns) * float(_get(par, "kd0"))

    ka_vec = Ki_vec * kd_vec  # (Ns,)

    ka = np.broadcast_to(ka_vec, (N, Ns))
    kd = np.broadcast_to(kd_vec, (N, Ns))

    Qsum      = np.sum(Q, axis=1)             # (N,)
    Qfree     = np.maximum(qmax - Qsum, 0.0)  # (N,)
    Qfree_mat = Qfree[:, None]                # (N,1) -> broadcast

    dQdt    = ka * C * Qfree_mat - kd * Q
    dqhatdt = dQdt / Qscale

    D1 = _get(par, "D1")
    D2 = _get(par, "D2")
    dcdz  = D1.dot(c)
    d2cdz = D2.dot(c)

    beta = ((1.0 - eps) / eps) * (Qscale / Cscale)
    dcdt = -u * dcdz + Dax * d2cdz - beta * dqhatdt

    # inlet boundary: Dirichlet via ghost/upwind (scaled)
    Cin_phys = np.asarray(_get(par, "feed_profile")(t), dtype=float).reshape(-1)
    if Cin_phys.size != Ns:
        raise ValueError(f"feed_profile(t) must return length {Ns}, got {Cin_phys.size}")
    cin = Cin_phys / Cscale
    dcdt[0, :] += (-u/dz) * (c[0, :] - cin) + (Dax/dz**2) * (cin - c[0, :])

    return np.concatenate([dcdt.reshape(-1), dqhatdt.reshape(-1)])
