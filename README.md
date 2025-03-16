# HW5 Bit Battle Leaderboard

## Performance Results

*Last updated: 2025-03-14 20:53:38 EDT*

| Rank | Error | Execution Time (seconds) |
|------|-------|-------------------------|
| 1 | 1e-13 | 46.0854 |
| 2 | 2e-12 | 5.5940 |
| 3 | 5e-12 | 1.0862 |
| 4 | 5e-12 | 18.2781 |
| 5 | 5e-9 | 0.1831 |

## Solutions

Both problems are from the [SIAM 100-Digit Challenge](https://epubs.siam.org/doi/10.1137/1.9780898717969). You can find all problems and solutions in Bornemann, Folkmar; Laurie, Dirk; Wagon, Stan; Waldvogel, Jörg (2004). The SIAM 100-digit challenge: A study in high-accuracy numerical computing.

The beauty of the two homework problems is that they can be solved both by straightforward linear algebra methods and can be accelerated by using some advanced mathematical tools like complex analysis. We strongly recommend you to read the solutions in the `problem1_answer.pdf` and `problem2_answer.pdf` files and we can promise you that you will learn a lot from them.

### Problem 1

You can find the multiple solutions to this problem in the `problem1_answer.pdf` file.

In `problem1.py`, we implemented the power method combined with the Strebel’s summation formula to accelerate the convergence (See Section 3.5 in `problem1_answer.pdf`). It achieves the accuracy of `2e-16` within 0.074 seconds only for a  `513` by `513` matrix.

Since we are implementing the power method, it involves the matrix-vector multiplication like

$$
(M x)_j  = \sum_{k=1}^{\infty} M_{jk} x_k
$$

Euler-Maclaurin formula aims to convert the infinite sum into a finite sum. Let us first consider the integral. We aim to find a one-to-one mapping $\phi: [1,n] \to [1,\infty)$ such that by change of the variable:

$$
\int_1^\infty f(x) dx =  \int_1^n \phi'(x) f(\phi(x)) dx.
$$

Our idea is to convert the infinite sum into a finite sum:

$$
\sum_{k=1}^\infty f(k) = \sum_{k=1}^{n-1} \phi_n'(k) f(\phi_n(k))
$$

So why using the right-hand side is better than computing the left-hand side truncated at $n$? We need to find a function $\phi_n$ such that the finite sum is closer to the integral on the right-hand side.

By the Euler-Maclaurin formula, we have

$$
\sum_{k=0}^{n} f(k) = \int_0^{n} f(x)dx + \frac{f(0) + f(n)}{2} + \sum_{j=1}^p \frac{B_{2j}}{(2j)!}[f^{(2j-1)}(n) - f^{(2j-1)}(0)] + R_{2p}
$$

where $B_{2j}$ are the Bernoulli numbers and $R_{2p}$ is the remainder term such that

$$
|R_{2p}| \leq \frac{B_{2p}}{(2p)!} \int_0^n |f^{(2p)}(x)| dx.
$$

Our goal is to find a function $\phi_n$ such that it has the smallest error in the Euler-Maclaurin formula. 

Lemma (R. Strebel). Let $n$ be a positive integer and $f_n(x) = (x + n)^{-\alpha}$, $\alpha > 1$. The function

$$
\phi_n(\xi) = \frac{n^{1+\beta} (n - \xi)^{-\beta}}{\beta} - \frac{n}{\beta} - \frac{(1 + \beta) \xi^2}{2n}, \quad \beta = \frac{6}{\alpha - 1},
$$
is strictly increasing and maps $[0,n)$ onto $[0,\infty)$. Then

$$
\sum_{k=0}^{\infty} f_n(k) = \sum_{k=0}^{n-1} \phi_n'(k) \cdot f_n(\phi_n(k)) + O(n^{-3-\alpha}).
$$
For problem 1, we roughly have $\alpha = 4$.


Actually, choosing the following 
$$
\phi_{\exp}(\xi) = \exp\left(\frac{2}{(1-u)^2} -\frac{1}{2u^2}\right)
$$
will make even faster convergence though it is hard to prove.












### Problem 2

You can find the multiple solutions to this problem in the `problem2_answer.pdf` file.

A key trick is that you can construct the matrix `A` by Kronecker product. See `problem2.py` for the implementation.



