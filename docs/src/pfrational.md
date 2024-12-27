# Prime Factorization representation of Rational

For any positive integer $N \in \mathbb{Z}^+$, it have a unique prime factorization representation
```math
N = \prod_{i}p_i^{e_i},
```
where $i = 1, 2, \dots, \infty$, $p_i$ are primes, and $e_i \ge 0$. Actually, if we modify this theorem a little, we have the prime factorization representation of rationals:

For any positive rational $q \in \mathbb{Q}^+$, it have a unique representation
```math
q = \prod_{i}p_i^{e_i},
```
where $i = 1, 2, \dots, \infty$ and $e_i \in \mathbb{N}$ can be positve, negative or zero.

Using this representation, we can quickly evaluate $q_1 \cdot q_2$ and $q_1 / q_2$:
```math
q_1 = \prod_{i}p_i^{e_i(1)}, \quad q_2 = \prod_{i}p_i^{e_i(2)}, \\
```
```math
q_1\cdot q_2 = \prod_{i}p_i^{e_i(1)+e_i(2)}, \quad q_1\cdot q_2 = \prod_{i}p_i^{e_i(1)-e_i(2)}.
```

We can also define
```math
\mathrm{gcd}(q_1, q_2) = \prod_{i}p_i^{\min\{e_i(1), e_i(2)\}}, \quad \mathrm{lcm}(q_1, q_2) = \prod_{i}p_i^{\max\{e_i(1),e_i(2)\}}
```

The prime factorization representation is very suitable for computing the factorial $n!$ and binomial $\tbinom{n}{m}$. If we known the maximum $n$ in all of the calculations, then the maximum $p_i$ needed in the calculation satisfy $p_i \le n \lt p_{i+1}$.

However, in our calculation of Wigner-3nj symbols, we need calculate some addition and subtraction of the numbers.
Considering to calculate $q_1 + q_2 - q_3$, we should first calculate
```math
q_g = gcd(q_1, q_2, q_3).
```
Then $q_1/q_g, q_2/q_g, q_3/q_g$ become integers. Next, we convert them into `BigInt` and do addition and subtraction.

## The maximum exponent

To convert the prime factorization representation of Rational into `Rational{BigInt}`, we need to store the pow table of primes to quickly obtain the result of $p_i^{e_i}$.

[Legendre's formula](https://en.wikipedia.org/wiki/Legendre%27s_formula): For any prime number $p$ and any positive integer $n$, let $\nu_{p}(n)$ be the exponent of the largest power of $p$ that divides $n$ (that is, the $p$-adic valuation of $n$). Then
```math
\nu _{p}(n!)=\sum _{i=1}^{\infty }\left\lfloor {\frac {n}{p^{i}}}\right\rfloor
```

[Kummer's Theorem](https://en.wikipedia.org/wiki/Kummer%27s_theorem): For given integers $n \ge m \ge 0$ and a prime number $p$, the p-adic valuation $\nu_{p}\!{\tbinom{n}{m}}$ of the binomial coefficient $\tbinom{n}{m}$ is equal to the number of carries when m is added to n − m in base p.

If we use factorial, then for $p^a \le n \lt p^{a-1}$, we have
```math
\nu_p(n!) = 1+p+p^2+\cdots +p^{a-1} = \frac{p^a-1}{p-1}
```
For $n = p^a, p = 2$, the maximum exponent needed is $n - 1$.

If we use the binomial, we have
```math
\sup_{0\le k\le n} \nu_p(C_n^k) \le \lfloor log_p n \rfloor,
```
because the number of carries is never larger than the number of digits of the number $n$ in base $p$.

For $n = p^a, p = 2$, the maximum exponent needed is $\lfloor log_2 n \rfloor$.
