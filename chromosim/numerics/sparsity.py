from __future__ import annotations
import numpy as np
from scipy.sparse import diags, eye, bmat, csc_matrix

def jacobian_sparsity_multi(N: int, Ns: int) -> csc_matrix:
    """
    2*N*Ns x 2*N*Ns logical Jacobian sparsity for multi-species column model.

    State layout (row-major):
      y = [ c(0,0..Ns-1), ..., c(N-1,0..Ns-1),
            qhat(0,0..Ns-1), ..., qhat(N-1,0..Ns-1) ]
    """
    e = np.ones(N)
    T = diags([e, e, e], offsets=[-1, 0, 1], shape=(N, N), format="csc")

    # replicate tri-diagonal axial stencil per species
    Ac = bmat([[T if s == r else None for s in range(Ns)]
               for r in range(Ns)], format="csc")
    Nc = N * Ns
    I  = eye(Nc, format="csc")

    S = bmat([[Ac, I],
              [I,  I]], format="csc")
    S.data[:] = 1
    return S
