from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# pick params (single or multi)
from chromosim.models.params_ev_aex_single import params_ev_aex_single as make_params
# from chromosim.models.params_ev_aex import params_ev_aex as make_params

from chromosim.numerics.build_matrices import build_matrices, attach_jacobian_sparsity
from chromosim.models.ev_column_ode_multi import ev_column_ode_multi

def run():
    par = make_params()
    build_matrices(par)
    attach_jacobian_sparsity(par)

    N, Ns = par.N, par.Nsp
    NN    = N * Ns

    y0 = np.zeros(2*NN)  # [c(:); qhat(:)]

    def rhs(t, y):
        return ev_column_ode_multi(t, y, par)

    sol = solve_ivp(
        rhs,
        t_span=par.tspan,
        y0=y0,
        method="BDF",
        vectorized=False,
        rtol=par.rtol,
        atol=par.atol,
        max_step=par.max_step,
        jac_sparsity=par.Jpattern
    )
    t = sol.t
    Y = sol.y.T

    # unpack scaled → physical
    cScal = Y[:, :NN].reshape(-1, N, Ns)
    qScal = Y[:, NN:].reshape(-1, N, Ns)
    Cphys = cScal * par.Cscale
    # Qphys = qScal * par.Qscale  # not used below, but available

    # outlet
    Cout_all = Cphys[:, -1, :]           # Nt × Ns
    Cout_tot = np.sum(Cout_all, axis=1)

    # gradient %B (mask until grad_start)
    I_vec = par.salt_profile(t)
    pctB  = par.B_start + (par.B_end - par.B_start) * \
            (I_vec - par.I_load)/(par.I_elute - par.I_load)
    pctB  = np.clip(pctB, par.B_start, par.B_end)
    pctB[t < par.grad_start] = np.nan  # hide before ramp

    # ======== Plots ========
    tMin = t/60
    z_edges_mm = np.linspace(0.0, par.L, par.N+1) * 1e3
    Ct_sum = np.sum(Cphys, axis=2)      # Nt × N

    # Figure 1: breakthrough + gradient, plus map
    plt.figure(figsize=(8,7))

    ax1 = plt.subplot(2,1,1)
    # left axis
    ax1.plot(tMin, Cout_tot/par.CplotScale, 'k-', lw=1.7, label='Total')
    if Ns > 1:
        cols = plt.cm.tab10.colors
        for s in range(Ns):
            ax1.plot(tMin, Cout_all[:,s]/par.CplotScale, color=cols[s%10], lw=1.2,
                     label=f"Zeta={par.zeta_mV[s]:.0f} mV")
    ax1.set_ylabel(f"C_out ({par.CplotLabel})")
    ax1.grid(True)
    # right axis for gradient
    ax1b = ax1.twinx()
    ax1b.plot(tMin, pctB, 'k--', lw=1.2, label='%B')
    ax1b.set_ylabel('%B (1 M NaCl)')
    ax1b.set_ylim(0,105)
    ax1.set_xlim(0, tMin[-1])
    ax1.set_title('Breakthrough & programmed salt gradient')
    # build one legend that includes both handles
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax1b.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc='upper left')

    # map
    ax2 = plt.subplot(2,1,2)
    T, Z = np.meshgrid(tMin, z_edges_mm, indexing='ij')
    # pcolormesh needs cell-centered to edge mapping: repeat last column for edges
    Ct_edges = np.hstack([Ct_sum, Ct_sum[:,[-1]]])
    pcm = ax2.pcolormesh(Z, T, Ct_edges/par.CplotScale, shading='auto')
    ax2.set_xlabel('Axial position z (mm)')
    ax2.set_ylabel('Time (min)')
    ax2.set_xlim(0, par.L*1e3)
    ax2.set_title('In-column total EV concentration')
    cbar = plt.colorbar(pcm, ax=ax2)
    cbar.set_label(par.CplotLabel)

    plt.tight_layout()

    # Figure 2: FPLC-style (mL, mAU-equivalent) with gradient
    plt.figure(figsize=(7,4))
    flow_mL_min = 1.0
    vol_mL = tMin * flow_mL_min
    eps_det = 1e-13
    signal_mAU = Cout_tot * eps_det

    axL = plt.gca()
    axL.plot(vol_mL, signal_mAU, 'b-', lw=1.5)
    axL.set_xlabel('Elution volume (mL)')
    axL.set_ylabel('Simulated detector (pseudo-mAU)')
    axL.grid(True)

    axR = axL.twinx()
    axR.plot(vol_mL, pctB, 'k--', lw=1.2)
    axR.set_ylabel('%B (1 M NaCl)')
    axR.set_ylim(0,105)
    plt.title('Model-predicted FPLC trace with programmed gradient')

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run()
