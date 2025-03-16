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

Both problems are from the [SIAM 100-Digit Challenge](https://epubs.siam.org/doi/10.1137/1.9780898717969). 

The beauty of the two homework problems is that they can be solved both by straightforward linear algebra methods and can be accelerated by using some advanced mathematical tools like complex analysis. We strongly recommend you to read the solutions in the `problem1_answer.pdf` and `problem2_answer.pdf` files and we can promise you that you will learn a lot from them.

### Problem 1

You can find the multiple solutions to this problem in the `problem1_answer.pdf` file.

In `problem1.py`, we implemented the power method combined with the Euler-Maclaurin formula to accelerate the convergence (See Section 3.5 in `problem1_answer.pdf`).

It achieves the accuracy of `2e-16` within 0.074 seconds for a `2^9+1` by `2^9+1` matrix.

### Problem 2

You can find the multiple solutions to this problem in the `problem2_answer.pdf` file.

A key trick is that you can construct the matrix `A` by Kronecker product. See `problem2.py` for the implementation.



