# CG coefficient and Wigner symbols

## CG coefficient

Ref [^1], P240, Section 8.2.4, Formula (20).

```math
C_{j_1m_1j_2m_2}^{j_3m_3} = \delta_{m_3, m_1+m_2}
\left[
    \dfrac{\begin{pmatrix}2j_1 \\ J-2j_2\end{pmatrix}\begin{pmatrix}2j_2\\J-2j_3\end{pmatrix}}{\begin{pmatrix}J+1\\J-2j_3\end{pmatrix}\begin{pmatrix}2j_1\\j_1-m_1\end{pmatrix}\begin{pmatrix}2j_2 \\ j_2 - m_2\end{pmatrix}\begin{pmatrix}2j_3 \\ j_3-m_3\end{pmatrix}}
\right]^{1/2} \\
\times \sum_{z} (-1)^{z}\begin{pmatrix}J-2j_3 \\ z\end{pmatrix}\begin{pmatrix}J-2j_2 \\ j_1-m_1-z\end{pmatrix}\begin{pmatrix}J-2j_1 \\ j_2 + m_2 - z\end{pmatrix}
```

where, $J = j_1+j_2+j_3$. It is already combination of binomials.

## 3j symbol

Ref [^1], P236, Section 8.1.2, Formula (11).

```math
\begin{aligned}
\begin{pmatrix}j_1 & j_2 & j_3 \\ m_1 & m_2 & m_3\end{pmatrix} &= (-1)^{j_3+m_3+2j_1}\dfrac{1}{\sqrt{2j_3+1}}C_{j_1-m_1j_2-m_2}^{j_3m_3} \\
&= (-1)^{j_3+m_3+2j_1} \dfrac{\delta_{-m_3, m_1+m_2}}{\sqrt{J+1}} \sum_{z} (-1)^{z}\begin{pmatrix}J-2j_3 \\ z\end{pmatrix}\begin{pmatrix}J-2j_2 \\ j_1 + m_1-z\end{pmatrix}\begin{pmatrix}J-2j_1 \\ j_2 - m_2 - z\end{pmatrix} \\
&\times \left[
    \dfrac{\begin{pmatrix}2j_1 \\ J-2j_2\end{pmatrix}\begin{pmatrix}2j_2\\J-2j_3\end{pmatrix}}{\begin{pmatrix}J\\J-2j_3\end{pmatrix}\begin{pmatrix}2j_1\\j_1 - m_1\end{pmatrix}\begin{pmatrix}2j_2 \\ j_2 - m_2\end{pmatrix}\begin{pmatrix}2j_3 \\ j_3-m_3\end{pmatrix}}\right]^{1/2}.
\end{aligned}
```

## 6j symbol

Ref [^1], P293, Section 9.2.1, Formula (1).

```math
\begin{Bmatrix}j_1 & j_2 & j_3 \\ j_4 & j_5 & j_6\end{Bmatrix} = \Delta(j_1j_2j_3)\Delta(j_4j_5j_3)\Delta(j_1j_5j_6)\Delta(j_4j_2j_6) \\
\times \sum\limits_x\dfrac{(-1)^x(x+1)!}{(x-j_{123})!(x-j_{453})!(x-j_{156})!(x-j_{426})!(j_{1245}-x)!(j_{1346}-x)!(j_{2356}-x)!}
```

Here, I use $j_{123} \equiv j_1+j_2+j_3$ for simplicity. The symbol $\Delta(abc)$ is defined as

```math
\Delta(abc) = \left[\dfrac{(a+b-c)!(a-b+c)!(-a+b+c)!}{(a+b+c+1)!}\right]^{\frac{1}{2}}.
```

Move $\Delta^2(j_1j_2j_3)$ into the summation, we can find that
```math
\begin{pmatrix}j_1+j_2-j_3 \\ x - j_{453}\end{pmatrix} = \dfrac{(j_1+j_2-j_3)!}{(x-j_{453})!(j_{1245}-x)!} \\
\begin{pmatrix}j_1-j_2+j_3 \\ x - j_{426}\end{pmatrix} = \dfrac{(j_1-j_2+j_3)!}{(x-j_{426})!(j_{1346}-x)!} \\
\begin{pmatrix}j_2+j_3-j_1 \\ x - j_{156}\end{pmatrix} = \dfrac{(j_2+j_3-j_1)!}{(x-j_{156})!(j_{2356}-x)!} \\
\begin{pmatrix}x+1 \\ j_{123}+1\end{pmatrix} = \dfrac{(x+1)!}{(x-j_{123})(j_{123}+1)!}.
```

So, we have

```math
\begin{Bmatrix}j_1 & j_2 & j_3 \\ j_4 & j_5 & j_6\end{Bmatrix} = \dfrac{\Delta(j_4j_5j_3)\Delta(j_1j_5j_6)\Delta(j_4j_2j_6)}{\Delta(j_1j_2j_3)} \\
\times \sum\limits_x (-1)^x \binom{x+1}{j_{123}+1} \binom{j_1+j_2-j_3}{x - j_{453}}
\binom{j_1-j_2+j_3}{x - j_{426}}\binom{j_2+j_3-j_1}{x - j_{156}}.
```

The $\Delta(abc)$ can also be rewrite with binomials,

```math
\Delta(abc) = \left[\dfrac{1}{\begin{pmatrix}a+b+c+1 \\ 2a + 1\end{pmatrix} \begin{pmatrix}2a \\ a + b - c\end{pmatrix}(2a+1)}\right]^{\frac{1}{2}}.
```

So

```math
\dfrac{\Delta(j_4j_5j_3)\Delta(j_1j_5j_6)\Delta(j_4j_2j_6)}{\Delta(j_1j_2j_3)} \\
= \dfrac{1}{2j_4+1} \left[\dfrac{\begin{pmatrix}j_{123}+1 \\ 2j_1+1\end{pmatrix}\begin{pmatrix}2j_1 \\ j_1 + j_2 - j_3\end{pmatrix}}{\begin{pmatrix}j_{156}+1 \\ 2j_1+1\end{pmatrix}\begin{pmatrix}2j_1 \\ j_1+j_5-j_6\end{pmatrix}\begin{pmatrix}j_{453}+1\\2j_4+1\end{pmatrix}\begin{pmatrix}2j_4\\j_4+j_5-j_3\end{pmatrix}\begin{pmatrix}j_{426}+1 \\ 2j_4+1\end{pmatrix}\begin{pmatrix}2j_4 \\ j_4+j_2-j_6\end{pmatrix}}\right]^{\frac{1}{2}}
```

## Racha coefficient

Ref [^1], P291, Section 9.1.2, Formula (11)

```math
W(j_1j_2j_3j_4, j_5j_6) = (-1)^{j_1+j_2+j_3+j_4} \begin{Bmatrix}
j_1 & j_2 & j_5 \\
j_4 & j_3 & j_6
\end{Bmatrix}
```

This package uses the 6j symbol above to calculate Racha coefficient. 

## 9j symbol

Ref [^1], P340, Section 10.2.4, Formula (20)

```math
\begin{aligned}
\begin{Bmatrix}j_1 & j_2 & j_3 \\ j_4 & j_5 & j_6 \\ j_7 & j_8 & j_9\end{Bmatrix} &= \sum\limits_{t}(-1)^{2t}(2t+1)\begin{Bmatrix}j_1 & j_2 & j_3 \\ j_6 & j_9 & t\end{Bmatrix} \begin{Bmatrix}j_4 & j_5 & j_6 \\ j_2 & t & j_8\end{Bmatrix} \begin{Bmatrix}j_7 & j_8 & j_9 \\ t & j_1 & j_4\end{Bmatrix} \\
&= \sum\limits_{t}(-1)^{2t}(2t+1)\begin{Bmatrix}j_1 & j_2 & j_3 \\ j_6 & j_9 & t\end{Bmatrix} \begin{Bmatrix}j_6 & j_4 & j_5 \\ j_8 & j_2 & t\end{Bmatrix} \begin{Bmatrix}j_8 & j_9 & j_7 \\ j_1 & j_4 & t\end{Bmatrix}.
\end{aligned}
```

To obtain a better formula for calculation, we take an alternative factorization of 6j symbol
```math
\begin{Bmatrix}j_1 & j_2 & j_3 \\ j_4 & j_5 & j_6\end{Bmatrix} = \frac{\Delta(j_1j_2j_3)\Delta(j_4j_2j_6)\Delta(j_4j_5j_3)}{\Delta(j_1j_5j_6)}\sum_{x}\binom{x+1}{j_{156}+1}\binom{j_1+j_5-j_6}{x - j_{426}}\binom{j_1+j_6-j_5}{x - j_{453}}\binom{j_5+j_6-j_1}{x - j_{123}}.
```

So, we can obtain
```math
\begin{aligned}
\begin{Bmatrix}j_1 & j_2 & j_3 \\ j_4 & j_5 & j_6 \\ j_7 & j_8 & j_9\end{Bmatrix}
&= P_0  \sum_{t}(-1)^{2t}(2t+1)\left(\sum_{x}A(t,x)\right)\left(\sum_{x}B(t,y)\right)\left(\sum_{x}C(t,z)\right)
\end{aligned}
```
where
```math
\begin{aligned}
P_0 &= \Delta(j_1j_2j_3)\Delta(j_4j_5j_6)\Delta(j_7j_8j_9)\Delta(j_1j_4j_7)\Delta(j_2j_5j_8)\Delta(j_3j_6j_9), \\
A(t, x) &= (-1)^x\binom{x+1}{j_{19t}+1}\binom{j_1+j_9-t}{x - j_{26t}}\binom{j_1+t-j_9}{x - j_{369}}\binom{t+j_9-j_1}{x - j_{123}}, \\
B(t, y) &= (-1)^y\binom{y+1}{j_{26t}+1}\binom{j_6+j_2-t}{y - j_{48t}}\binom{j_6+t-j_2}{y - j_{258}}\binom{t+j_2-j_6}{y - j_{456}}, \\
C(t, z) &= (-1)^z\binom{z+1}{j_{48t}+1}\binom{j_8+j_4-t}{z - j_{19t}}\binom{j_8+t-j_4}{z - j_{147}}\binom{t+j_4-j_8}{z - j_{789}}.
\end{aligned}
```

## Estimate the capacity

Assume we are doing a calculation for a system, we will not calculate the Winger Symbols with very large angular momentum, because usually we will take some trunction. If in such truncated system, the max angular momentum is $J_{max}$, now we can estimate how many binomial coefficients we need to store to compute those Wigner Symbols.

### With maximum $J_{max}$
####  CG & 3j
According to the formula, the max possible binomial is $\begin{pmatrix} J+1 \\ J-2j_3\end{pmatrix}$, where $J = j_1+j_2+j_3$. So we need at least store binomials to $\begin{pmatrix} n_{min} \\ k\end{pmatrix}$ where
```math
n_{min} \geq 3J_{max}+1
```

#### 6j & Racha
For binomials in `sqrt`，the maximum possible $n = 3J_{max}+1$, and for those in the summation, we have to calculate the boundary of $x + 1$. According to the formula
```math
x \leq \min\{j_1+j_2+j_4+j_5, j_1+j_3+j_4+j_6, j_2+j_3+j_5+j_6\} \leq 4J_{max}
```
So we have to at least store to
```math
n_{min} \geq 4J_{max} + 1
```

#### 9j

Similarly, according to $P_0$, we need $3J_{max}+1$. According to $P(t)$,
```math
t \leq \min\{j_2 + j_6, j_4 + j_8, j_1 + j_9\} \leq 2J_{max}
```
so we need $j_1 + j_9 + t + 1 \leq 4J_{max} + 1$. According to $A(t, x),\; B(t, y),\; C(t, z)$,
```math
x \leq \min\{j_1+j_2+j_6+j_9, j_1+j_3+j_6+t, j_2+j_3+j_9+t\} \leq 5J_{max}
```
so we need to at least store to
```math
n_{min} \geq 5J_{max} + 1
```

### With maximum single particle $j_{max}$

Considering a quantum many body calculation, we can get a more optimistic estimation of the calculation capacity. In a quantum many body calculation, we often truncate single particle orbits, and in this condition we assume that we only need to consider two-body coupled angular momentum.

The maximum angular momentum of the single particle orbits defined as $j_{max}$. If for each function, there are at least one of the parameters is single particle angular momentum, then we can get the following estimation.

#### CG & 3j

In this condition, the $j_1, j_2, j_3$ contains two single particle angular momentum, and the rest one is a coupled two body angular momentum. So
```math
n_{min} \geq 4j_{max} + 1
```

#### 6j & Racha

Similarly
```math
x \leq \min\{j_1+j_2+j_4+j_5, j_1+j_3+j_4+j_6, j_2+j_3+j_5+j_6\}
```
These summations of four angular momentum can be split into two kinds, all single particle angular momentum, and two single particle with two two body angular momentum. Consider the worst condition
```math
n_{min} \geq 6j_{max} + 1
```

### 9j

We assume the 9j symbol contains at least one single particle angular momentum, then according to the coupling rules, the 9j symbol can be split into kinds
```math
\begin{Bmatrix}
1 & 1 & 2 \\
1 & 1 & 2 \\
2 & 2 & 2
\end{Bmatrix}
\quad
\begin{Bmatrix}
1 & 1 & 2 \\
1 & 2 & 1 \\
2 & 1 & 1
\end{Bmatrix}
```
where $1$ represents single particle and $2$ for two body angular momentum.
In both cases, we can deduce that
```math
n_{min} \geq 8j_{max} + 1
```

----
Reference

[^1]: A. N. Moskalev D. A. Varshalovich and V. K. Khersonskii, *Quantum theory of angular momentum*.
