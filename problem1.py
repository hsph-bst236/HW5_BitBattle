import numpy as np
import time
def SummationFormula(n, alpha):
    """
    Calculates weights w and nodes c for a summation formula that
    approximates sum(f(k),k=1..infinity) by sum(w(k)*f(c(k)),k=1..n).
    Works well if f(k) is asymptotic to k^(-alpha) for large k.
    Put alpha = 'exp' to use an exponential formula.
    
    Parameters:
    -----------
    n : int
        Number of terms in the approximation
    alpha : float or str
        Power in the asymptotic behavior or 'exp' for exponential formula
        
    Returns:
    --------
    w : numpy.ndarray
        Weights for the summation formula
    c : numpy.ndarray
        Nodes for the summation formula
    """
    
    if alpha == 'exp':
        k = np.arange(n, 1, -1)  # n down to 2
        u = (k - 1) / n
        c1 = np.exp(2 / (1 - u)**2 - 1 / u**2 / 2)
        c = k + c1
        w = 1 + c1 * (4 / (1 - u)**3 + 1 / u**3) / n
        
        finite_indices = np.isfinite(c)
        c = c[finite_indices]
        w = w[finite_indices]
        
        c = np.append(c, 1)
        w = np.append(w, 1)
    else:
        alpha = float(alpha)
        n = int(np.ceil(n / 2))
        a1 = (alpha - 1) / 6
        a6 = 1 + 1 / a1
        
        k2 = np.arange(n - 1, -1, -1)  # (n-1) down to 0
        w2 = n**a6 / (n - k2)**a6 - a6 * k2 / n
        c2 = n + a1 * (-n + n**a6 / (n - k2)**(1 / a1)) - a6 * k2**2 / 2 / n
        
        k1 = np.arange(n - 1, 0, -1)  # (n-1) down to 1
        
        w = np.concatenate([w2, np.ones_like(k1)])
        c = np.concatenate([c2, k1])
    
    return w, c 

def sig1_M(m = 9) -> float:
    n = 2**m + 1
    lambda1 = 0
    lambda0 = 1
    tol = 1e-16
    alpha = 4 # 'exp' or 4
    w, c = SummationFormula(n, alpha)
    x = np.ones(len(w))
    y = np.zeros(len(w))
 
    
    while abs(lambda1 - lambda0) > tol:
        lambda0 = lambda1
        x_new = x / np.linalg.norm(x)
        for j in range(len(x)):
            y[j] = np.sum(w * (1/((c[j] + c - 1) * (c[j] + c)/2 - c + 1))  * x_new)
        for j in range(len(x)):
            x[j] = np.sum(w * (1/((c[j] + c - 1) * (c[j] + c)/2 - c[j] + 1)) * y)
        lambda1 = np.sum(x * x_new)
    s = np.sqrt(lambda1)
    return s

if __name__ == "__main__":
    TRUTH = 1.2742241528212282123
    start_time = time.time()
    s_max = sig1_M()
    end_time = time.time()
    print(f"s_max: {s_max:.16f}, error: {abs(s_max - TRUTH):.1e}")
    print(f"Time taken: {end_time - start_time} seconds")
 