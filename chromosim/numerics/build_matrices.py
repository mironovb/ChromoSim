from __future__ import annotations
import numpy as np
from scipy.sparse import diags, csc_matrix

def build_matrices(par):
    """
    Assemble finite-volume derivative matrices in z and attach to `par`.
      par.N, par.dz must exist.
    Adds:
      par.D1  -- first derivative (backward/upwind)  [N x N] sparse
      par.D2  -- second derivative (central)         [N x N] sparse
    """
    N  = int(getattr(par, "N", par["N"]))
    dz = float(getattr(par, "dz", par["dz"]))

    e = np.ones(N)
    # Backward/upwind first derivative
    D1 = diags([-e, e], offsets=[-1, 0], shape=(N, N), format="csc") / dz
    D1 = D1.tolil()
    D1[0, :] = 0  # inlet row handled via boundary condition
    D1 = D1.tocsc()

    # Central second derivative
    D2 = diags([e, -2*e, e], offsets=[-1, 0, 1], shape=(N, N), format="csc") / (dz**2)
    D2 = D2.tolil()
    D2[0, 0] = -2.0 / (dz**2)  # ghost for inlet Dirichlet
    D2[-1, -1] = 1.0 / (dz**2) # Danckwerts outlet
    D2 = D2.tocsc()

    setattr(par, "D1", D1)
    setattr(par, "D2", D2)
    return par

# Optional helper to attach Jacobian sparsity later
from .sparsity import jacobian_sparsity_multi

def attach_jacobian_sparsity(par):
    N  = int(getattr(par, "N", par["N"]))
    Ns = int(getattr(par, "Nsp", par["Nsp"]))
    J  = jacobian_sparsity_multi(N, Ns)
    setattr(par, "Jpattern", J)
    return par
