# Formula

## CG coefficient

Ref [^1], P240, Section 8.2.4, Formula (20).

```math
C_{j_1m_1j_2m_2}^{j_3m_3} = \delta_{m_3, m_1+m_2}\left[\dfrac{\begin{pmatrix}2j_1 \\ J-2j_3\end{pmatrix}\begin{pmatrix}2j_2\\J-2j_3\end{pmatrix}}{\begin{pmatrix}J+1\\J-2j_3\end{pmatrix}\begin{pmatrix}2j_1\\j_1-m_1\end{pmatrix}\begin{pmatrix}2j_2 \\ j_2 - m_2\end{pmatrix}\begin{pmatrix}2j_3 \\ j_3-m_3\end{pmatrix}}\right]^{1/2} \sum_{z} (-1)^{z}\begin{pmatrix}J-2j_3 \\ z\end{pmatrix}\begin{pmatrix}J-2j_2 \\ j_1-m_1-z\end{pmatrix}\begin{pmatrix}J-2j_1 \\ j_2 + m_2 - z\end{pmatrix}
```

where, $J = j_1+j_2+j_3$. It is already combination of binominals.

## 3j symbol

Ref [^1], P236, Section 8.1.2, Formula (11).

```math
\begin{pmatrix}j_1 & j_2 & j_3 \\ m_1 & m_2 & m_3\end{pmatrix} = (-1)^{j_3+m_3+2j_1}\dfrac{1}{\sqrt{2j_3+1}}C_{j_1-m_1j_2-m_2}^{j_3m_3}
```

This package use the CG coefficient above to calculate 3j symbol.

## 6j symbol

Ref [^1], P293, Section 9.2.1, Formula (1).

```math
\begin{Bmatrix}j_1 & j_2 & j_3 \\ j_4 & j_5 & j_6\end{Bmatrix} = \Delta(j_1j_2j_3)\Delta(j_4j_5j_3)\Delta(j_1j_5j_6)\Delta(j_4j_2j_6) \\
\times \sum\limits_n\dfrac{(-1)^n(n+1)!}{(n-j_1j_2j_3)!(n-j_4j_5j_3)!(n-j_1j_5j_6)!(n-j_4j_2j_6)!(j_1j_2j_4j_5-n)!(j_1j_3j_4j_6-n)!(j_2j_3j_5j_6-n)!}
```

Here, I use $j_1j_2j_3$ to represent $j_1+j_2+j_2$. The symbol $\Delta(abc)$ is defined as

```math
\Delta(abc) = \left[\dfrac{(a+b-c)!(a-b+c)!(-a+b+c)!}{(a+b+c+1)!}\right]^{\frac{1}{2}}.
```

We can find that

```math
\begin{pmatrix}j_1+j_2-j_3 \\ n - j_4j_5j_3\end{pmatrix} = \dfrac{(j_1+j_2-j_3)!}{(n-j_4j_5j_3)!(j_1j_2j_4j_5-n)!} \\
\begin{pmatrix}j_1-j_2+j_3 \\ n - j_4j_2j_6\end{pmatrix} = \dfrac{(j_1-j_2+j_3)!}{(n-j_4j_2j_6)!(j_1j_3j_4j_6-n)!} \\
\begin{pmatrix}j_2+j_3-j_1 \\ n - j_1j_5j_6\end{pmatrix} = \dfrac{(j_2+j_3-j_1)!}{(n-j_1j_5j_6)!(j_2j_3j_5j_6-n)!} \\
\begin{pmatrix}n+1 \\ j_1j_2j_3+1\end{pmatrix} = \dfrac{(n+1)!}{(n-j_1j_2j_3)(j_1j_2j_3+1)!}.
```

So, we have

```math
\begin{Bmatrix}j_1 & j_2 & j_3 \\ j_4 & j_5 & j_6\end{Bmatrix} = \dfrac{\Delta(j_4j_5j_3)\Delta(j_1j_5j_6)\Delta(j_4j_2j_6)}{\Delta(j_1j_2j_3)} \\
\times \sum\limits_n (-1)^n \begin{pmatrix}n+1 \\ j_1j_2j_3+1\end{pmatrix} \begin{pmatrix}j_1+j_2-j_3 \\ n - j_4j_5j_3\end{pmatrix} \begin{pmatrix}j_1-j_2+j_3 \\ n - j_4j_2j_6\end{pmatrix} \begin{pmatrix}j_2+j_3-j_1 \\ n - j_1j_5j_6\end{pmatrix}.
```

Rewrite $\Delta(abc)$ with binominals,

```math
\Delta(abc) = \left[\dfrac{1}{\begin{pmatrix}a+b+c+1 \\ 2a + 1\end{pmatrix} \begin{pmatrix}2a \\ a + b - c\end{pmatrix}(2a+1)}\right]^{\frac{1}{2}}.
```

So

```math
\dfrac{\Delta(j_4j_5j_3)\Delta(j_1j_5j_6)\Delta(j_4j_2j_6)}{\Delta(j_1j_2j_3)} \\
= \dfrac{1}{2j_4+1} \left[\dfrac{\begin{pmatrix}j_1j_2j_3+1 \\ 2j_1+1\end{pmatrix}\begin{pmatrix}2j_1 \\ j_1 + j_2 - j_3\end{pmatrix}}{\begin{pmatrix}j_1j_5j_6+1 \\ 2j_1+1\end{pmatrix}\begin{pmatrix}2j_1 \\ j_1+j_5-j_6\end{pmatrix}\begin{pmatrix}j_4j_5j_3+1\\2j_4+1\end{pmatrix}\begin{pmatrix}2j_4\\j_4+j_5-j_3\end{pmatrix}\begin{pmatrix}j_4j_2j_6+1 \\ 2j_4+1\end{pmatrix}\begin{pmatrix}2j_4 \\ j_4+j_2-j_6\end{pmatrix}}\right]^{\frac{1}{2}}
```

## Racha coefficient

Ref [^1], P291, Section 9.1.2, Formula (11)

```math
W(j_1j_2j_3j_4, j_5j_6) = (-1)^{j_1+j_2+j_3+j_4} \begin{Bmatrix}
j_1 & j_2 & j_5 \\
j_4 & j_3 & j_6
\end{Bmatrix}
```

This package use the 6j symbol above to calculate Racha coefficient. 

## 9j symbol

Ref [^1], P340, Section 10.2.4, Formula (20)

```math
\begin{Bmatrix}j_1 & j_2 & j_3 \\ j_4 & j_5 & j_6 \\ j_7 & j_8 & j_9\end{Bmatrix} = \sum\limits_{t}(-1)^{2t}(2t+1)\begin{Bmatrix}j_1 & j_2 & j_3 \\ j_6 & j_9 & t\end{Bmatrix} \begin{Bmatrix}j_4 & j_5 & j_6 \\ j_2 & t & j_8\end{Bmatrix} \begin{Bmatrix}j_7 & j_8 & j_9 \\ t & j_1 & j_4\end{Bmatrix}
```

Use the 6j symbol result above, we get

```math
\dfrac{\Delta(j_1j_9t)\Delta(j_6j_9j_3)\Delta(j_6j_2t)}{\Delta(j_1j_2j_3)} \dfrac{\Delta(j_4tj_8)\Delta(j_2tj_6)\Delta(j_2j_5j_8)}{\Delta(j_4j_5j_6)} \dfrac{\Delta(j_7j_1j_4)\Delta(tj_1j_9)\Delta(tj_8j_4)}{\Delta(j_7j_8j_9)} \\
 = \dfrac{\Delta(j_3j_6j_9)\Delta(j_2j_5j_8)\Delta(j_1j_4j_7)}{\Delta(j_1j_2j_3)\Delta(j_4j_5j_6)\Delta(j_7j_8j_9)} \Delta^2(j_1j_9t)\Delta^2(j_2j_6t)\Delta^2(j_4j_8t).
```

Define $j_{123} \equiv j_1+j_2+j_3$, and

```math
P_0 \equiv \dfrac{\Delta(j_3j_6j_9)\Delta(j_2j_5j_8)\Delta(j_1j_4j_7)}{\Delta(j_1j_2j_3)\Delta(j_4j_5j_6)\Delta(j_7j_8j_9)} \\
= \left[\dfrac{\begin{pmatrix}j_{123} + 1 \\ 2j_1+1\end{pmatrix}\begin{pmatrix}2j_1\\j_1+j_2-j_3\end{pmatrix}\begin{pmatrix}j_{456}+1\\2j_5+1\end{pmatrix}\begin{pmatrix}2j_5 \\ j_4+j_5-j_6\end{pmatrix}\begin{pmatrix}j_{789}+1 \\ 2j_9+1\end{pmatrix}\begin{pmatrix}2j_9 \\ j_7 + j_9 - j_8\end{pmatrix}}{\begin{pmatrix}j_{147} + 1\\ 2j_1 + 1\end{pmatrix}\begin{pmatrix}2j_1 \\ j_1+j_4-j_7\end{pmatrix}\begin{pmatrix}j_{258}+1 \\ 2j_5+1\end{pmatrix}\begin{pmatrix}2j_5 \\ j_2+j_5 - j_8\end{pmatrix}\begin{pmatrix}j_{369}+1 \\ 2j_9+1\end{pmatrix}\begin{pmatrix}2j_9 \\ j_3+j_9 - j_6\end{pmatrix}}\right]^{1/2}.
```

Define

```math
P(t) \equiv (-1)^{2t}(2t+1)\Delta^2(j_1j_9t)\Delta^2(j_2j_6t)\Delta^2(j_4j_8t)  \\
 = (-1)^{2t} \dfrac{1}{\begin{pmatrix}j_1+j_9+t+1 \\ 2t+1\end{pmatrix}\begin{pmatrix}2t \\ j_1+t-j_9\end{pmatrix}\begin{pmatrix}j_2+j_6+t+1 \\ 2t+1\end{pmatrix}\begin{pmatrix}2t \\ j_2+t-j_6\end{pmatrix}\begin{pmatrix}j_4+j_8+t \\ 2t+1\end{pmatrix}\begin{pmatrix}2t \\ j_4+t-j_8\end{pmatrix}(2t+1)^2}.
```

```math
A(t,x) \equiv (-1)^x\begin{pmatrix}x+1 \\ j_{123} + 1\end{pmatrix}\begin{pmatrix}j_1+j_2-j_3 \\ x - (j_6+j_9+j_3)\end{pmatrix} \begin{pmatrix}j_1+j_3-j_2 \\ x - (j_6+j_2+t)\end{pmatrix}\begin{pmatrix}j_2+j_3-j_1 \\ x - (j_1 + j_9 + t)\end{pmatrix}.
```

```math
B(t,y) \equiv (-1)^y\begin{pmatrix}y+1 \\ j_{456} + 1\end{pmatrix}\begin{pmatrix}j_4+j_5-j_6 \\ y - (j_2+t+j_6)\end{pmatrix}\begin{pmatrix}j_4+j_6-j_5 \\ y - (j_2+j_5+j_8)\end{pmatrix}\begin{pmatrix}j_5+j_6-j_4 \\ y - (j_4+t+j_8)\end{pmatrix}.
```

```math
C(t,z) \equiv (-1)^z \begin{pmatrix}z+1 \\ j_{789} + 1\end{pmatrix}\begin{pmatrix}j_7+j_8-j_9 \\ z - (t+j_1+j_9)\end{pmatrix}\begin{pmatrix}j_7+j_9-j_8 \\ z - (t+j_8+j_4)\end{pmatrix}\begin{pmatrix}j_8+j_9- j_7 \\ z - (j_7 + j_1 + j_4)\end{pmatrix}.
```

At last, we get

```math
\begin{Bmatrix}j_1 & j_2 & j_3 \\ j_4 & j_5 & j_6 \\ j_7 & j_8 & j_9\end{Bmatrix} = P_0\sum_t P(t) \left(\sum_x A(t,x)\right) \left(\sum_y B(t,y)\right) \left(\sum_z C(t,z)\right)
```

It deserves to be mentioned that, although the formula has 4 $\sum$s, the $\sum$ of $x,y,z$ are decoupled. So we can do the three `for loop`s respectively, which means the depth of `for loop` is not 4 but 2.


----
Reference

[^1]: A. N. Moskalev D. A. Varshalovich and V. K. Khersonskii, *Quantum theory of angular momentum*.

