import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import root_scalar
import time

def return_probability(epsilon, n):
    """
    Calculate return probability for a 2D grid with biased random walk.
    
    Parameters:
    -----------
    epsilon : float
        Bias parameter for east-west movement
    n : int
        Grid size parameter (resulting grid will be of size 2n+1 x 2n+1)
        
    Returns:
    --------
    p : float
        Return probability to the center
    """
    pE = 1/4 + epsilon; pW = 1/4 - epsilon; pN = 1/4; pS = 1/4
    
    m = 2*n+1; ctr = (n+1)*m + (n+1)
    A_EW = sparse.spdiags(np.array([[pW]*m, [pE]*m]), [-1, 1], m, m)
    A_NS = sparse.spdiags(np.array([[pS]*m, [pN]*m]), [-1, 1], m, m)
    A = sparse.kron(A_EW, sparse.eye(m)) + sparse.kron(sparse.eye(m), A_NS)
    
    # Extract the column before modifying the matrix
    r = A[:, ctr].toarray().flatten()
    
    # Use lil_matrix format which is efficient for modifying the sparsity structure
    A_lil = A.tolil()
    A_lil[:, ctr] = 0
    A = A_lil.tocsr()  # Convert back to CSR for efficient computations
    start_time = time.time()
    q = spsolve(sparse.eye(m*m)-A, r)
    # q, info = sparse.linalg.gmres(sparse.eye(m*m)-A, r, rtol=1e-13)
    end_time = time.time()
    print(f"Time taken linear system: {end_time - start_time} seconds")
    p = q[ctr]
    
    return p

def Numeria_delta(n=40):
    def objective(delta): return return_probability(delta, n) - 1/2
    result = root_scalar(objective, bracket=[0, 0.25], method='brentq')
    return result.root

if __name__ == "__main__":
    n = 40
    start_time = time.time()
    print(f"delta for n = {n}: {Numeria_delta(n)}")
    end_time = time.time()
    print(f"Time taken: {end_time - start_time} seconds")


