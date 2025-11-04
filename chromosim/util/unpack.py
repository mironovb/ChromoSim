from __future__ import annotations
import numpy as np

def unpack_state(Y: np.ndarray, par):
    """
    Reshape solution (Nt × (2*N*Ns)) into:
        C : Nt × N × Ns   (mobile phase, physical or scaled depending on caller)
        q : Nt × N × Ns   (bound phase)
    """
    Nt, total = Y.shape
    N, Ns = par.N, par.Nsp
    NN = N * Ns
    C = Y[:, :NN].reshape(Nt, N, Ns)
    q = Y[:, NN:].reshape(Nt, N, Ns)
    return C, q
