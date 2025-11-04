from __future__ import annotations
import numpy as np
from scipy.sparse import spdiags, kron, eye, bmat, csr_matrix

def sparsity_pattern_multi(N: int, Ns: int):
    """
    Boolean sparse pattern for the Jacobian of size 2*N*Ns × 2*N*Ns.
    State ordering:
        [ c(1..N, species 1..Ns) ; qhat(1..N, species 1..Ns) ]
    Blocks:
        dc/dc  ~ tri-diagonal in z for each species
        dc/dq  ~ diagonal (local kinetics)
        dq/dc  ~ diagonal (local kinetics)
        dq/dq  ~ diagonal
    """
    e = np.ones(N)
    T = spdiags([e, e, e], [-1, 0, 1], N, N, format='csr')  # tri-diagonal
    # replicate per species
    Ac = kron(eye(Ns, format='csr'), (T != 0), format='csr')  # Nc×Nc, Nc=N*Ns

    Nc = N * Ns
    B  = eye(Nc, format='csr', dtype=bool)    # diagonal
    Iq = eye(Nc, format='csr', dtype=bool)

    # assemble block matrix
    S = bmat([[Ac.astype(bool), B],
              [B,               Iq]], format='csr')
    return S
