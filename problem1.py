import numpy as np
import time

def sig1_M(m = 9) -> float:
    n = 2**m + 1
    x = np.ones(n)
    y = np.zeros(n)
    lambda1 = 0
    lambda0 = 1
    k = np.arange(1, n)
    tol = 1e-16
    alpha = 4
    n1 = int(np.ceil(n/2))
    a = (alpha - 1) / 6
    b = 1 + 1 / a
    k1 = np.arange(1, n1)
    k2 = np.arange(0, n1)
    w2 = n1 ** b / (n1 - k2) ** b - b * k2 / n1
    c2 = n1 + a * (-n1 + n1 ** b / (n1 - k2) ** (1/a)) - b * k2 ** 2 / (2 * n1)
    w = np.concatenate((np.ones(len(k1)), w2))
    c = np.concatenate((k1, c2))
 
    
    while abs(lambda1 - lambda0) > tol:
        lambda0 = lambda1
        x_new = x / np.linalg.norm(x)
        for j in range(n):
            y[j] = np.sum(w * (1/((c[j] + c - 1) * (c[j] + c)/2 - c + 1))  * x_new)
        for j in range(n):
            x[j] = np.sum(w * (1/((c[j] + c - 1) * (c[j] + c)/2 - c[j] + 1)) * y)
        lambda1 = np.sum(x * x_new)
    s = np.sqrt(lambda1)
    return s

if __name__ == "__main__":
    start_time = time.time()
    s_max = sig1_M()
    end_time = time.time()
    print(f"s_max: {s_max:.16f}")
    print(f"Time taken: {end_time - start_time} seconds")
    