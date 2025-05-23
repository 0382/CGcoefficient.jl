# Moshinsky

The Talmi-Moshinsky braket can be evaluated as following [^1]:
```math
\begin{aligned}
    &\langle n_3l_3,n_4l_4;\Lambda|n_1l_1,n_2l_2;\Lambda\rangle = \prod_{k=1}^4 i^{-l_k} \sqrt{\dfrac{n_k![2(n_k+l_k)+1]!!}{\sqrt{2^{l_k}}}} \\
    \times &\sum_{n_an_bn_cn_dl_a l_b l_c l_d} \Bigg[(-1)^{l_d} (\sin\beta)^{2n_a + l_a + 2n_d + l_d}(\cos\beta)^{2n_b+l_b + 2n_c + l_c} \begin{Bmatrix}l_a & l_b & l_1 \\ l_c & l_d & l_2 \\ l_3 & l_4 & \Lambda \end{Bmatrix} \\ 
    \times & C_{l_a0l_b0}^{l_10} C_{l_a0l_c0}^{l_30} C_{l_b0l_d0}^{l_40} C_{l_c0l_d0}^{l_20} \prod_{p = a}^{d} \dfrac{(-1)^{l_p}\sqrt{2^{l_p}}(2l_p+1)}{n_p![2(n_p+l_p) + 1]!!} \Bigg].
\end{aligned}
```
where
```math
\begin{aligned}
    2n_a + l_a + 2n_b + l_b = 2n_1 + l_1 \equiv e_1, \\
    2n_a + l_a + 2n_c + l_c = 2n_3 + l_3 \equiv e_3, \\
    2n_b + l_b + 2n_d + l_d = 2n_4 + l_4 \equiv e_4, \\
    2n_c + l_c + 2n_d + l_d = 2n_2 + l_2 \equiv e_2.
\end{aligned}
```

Here we define $\chi = e_1 + e_2 = e_3 + e_4$。

Similar as Wigner symbols, we try to rewrite it as binomials. First, we have:
```math
    (2k+1)!! = \dfrac{(2k+1)!}{(2k)!!} = \dfrac{(2k+1)!}{2^k k!}.
```
We can rewrite the following term:
```math
\begin{aligned}
    t_p &\equiv \dfrac{(-1)^{l_p}\sqrt{2^{l_p}}(2l_p+1)}{n_p![2(n_p+l_p) + 1]!!} \\
        &= \dfrac{(-1)^{l_p}\sqrt{2^{l_p}}(2l_p+1)2^{n_p+l_p}}{(2n_p+l_p+1)!}\begin{pmatrix}2n_p+l_p+1 \\ n_p\end{pmatrix} \begin{pmatrix}2(n_p+l_p)+1 \\ n_p + l_p\end{pmatrix}^{-1}.
\end{aligned}
```
Similiarly,
```math
\begin{aligned}
    r_k &\equiv \dfrac{n_k![2(n_k+l_k)+1]!!}{\sqrt{2^l_k}} \\
        &= \dfrac{(2n_k+l_k+2)!}{\sqrt{2^{l_k}}2^{n_k+l_k}(n_k+l_k+2)}\begin{pmatrix}2(n_k+l_k)+1 \\ n_k + l_k\end{pmatrix}\begin{pmatrix}2n_k+l_k+2 \\ n_k\end{pmatrix}^{-1}.
\end{aligned}
```
Next,
```math
\begin{aligned}
    t &\equiv \dfrac{\prod_{k=1}^4 \sqrt{(2n_k+l_k+2)!}}{\prod_{p=a}^d (2n_p+l_p + 1)!} \\
    &= \sqrt{\begin{pmatrix}2n_1+l_1+2 \\ 2n_a + l_a + 1\end{pmatrix}\begin{pmatrix}2n_2 + l_2 + 2 \\ 2n_c + l_c + 1\end{pmatrix} \begin{pmatrix}2n_3 + l_3 + 2 \\ 2n_a + l_a + 1\end{pmatrix}\begin{pmatrix}2n_4 + l_4 + 2 \\ 2n_b + l_b + 1\end{pmatrix}}
\end{aligned}
```

In the prod of $t_p$, there is
```math
    \prod_{p=a}^{d} \sqrt{2^{l_p}}2^{n_p+l_p} = \prod_{p=a}^d 2^{\frac{2n_p+l_p}{2}} 2^{l_p} = 2^{\chi/2} \prod_{p=a}^{d} 2^{l_p}.
```
In the pord of $r_k$, there is
```math
    \sqrt{\prod_{k=1}^4 \dfrac{1}{\sqrt{2^{l_k}}2^{n_k+l_k}}} = 2^{-\chi/2} \prod_{k=1}^4 2^{-l_k/2}
```

Finally, we have
```math
\begin{aligned}
    &\langle n_3l_3,n_4l_4;\Lambda|n_1l_1,n_2l_2;\Lambda\rangle\\
    = &\prod_{k=1}^4 i^{-l_k} \sqrt{\dfrac{1}{2^{l_k}(n_k+l_k+2)}\begin{pmatrix}2(n_k+l_k)+1 \\ n_k + l_k\end{pmatrix}\begin{pmatrix}2n_k+l_k+2 \\ n_k\end{pmatrix}^{-1}} \\
    \times & \sum_{n_an_bn_cn_dl_a l_b l_c l_d} \Bigg[(-1)^{l_d} (\sin\beta)^{2n_a + l_a + 2n_d + l_d}(\cos\beta)^{2n_b+l_b + 2n_c + l_c} \begin{Bmatrix}l_a & l_b & l_1 \\ l_c & l_d & l_2 \\ l_3 & l_4 & \Lambda \end{Bmatrix} \\
    \times & \sqrt{\begin{pmatrix}2n_1+l_1+2 \\ 2n_a + l_a + 1\end{pmatrix}\begin{pmatrix}2n_2 + l_2 + 2 \\ 2n_c + l_c + 1\end{pmatrix} \begin{pmatrix}2n_3 + l_3 + 2 \\ 2n_a + l_a + 1\end{pmatrix}\begin{pmatrix}2n_4 + l_4 + 2 \\ 2n_b + l_b + 1\end{pmatrix}} \\
    \times & C_{l_a0l_b0}^{l_10} C_{l_a0l_c0}^{l_30} C_{l_b0l_d0}^{l_40} C_{l_c0l_d0}^{l_20} \prod_{p=a}^d (-1)^{l_p}2^{l_p}(2l_p+1)\begin{pmatrix}2n_p+l_p+1 \\ n_p\end{pmatrix} \begin{pmatrix}2(n_p+l_p)+1 \\ n_p + l_p\end{pmatrix}^{-1} \Bigg]
\end{aligned}
```

Consider $m_1 = m_2,\;\omega_1 = \omega_2$, thus $\tan\beta = 1, \; \cos\beta = \sin\beta = 1/\sqrt{2}$, then
```math
\begin{aligned}
    (\sin\beta)^{2n_a + l_a + 2n_d + l_d}(\cos\beta)^{2n_b+l_b + 2n_c + l_c} &= 2^{-\frac{2n_a + l_a + 2n_b + l_b + 2n_c + l_c + 2n_d + l_d}{2}} \\
    &= 2^{-\frac{e_1+e_2}{2}} \\
    &= 2^{-\chi/2}
\end{aligned}
```

## Exact evaluation

Still start from the origin equation, define
```math
X_{abcd} = \left(\prod_{k=1}^4 i^{-l_k}\right)\left(\prod_{p = a}^d (-1)^{l_p}\right)\begin{Bmatrix}l_a & l_b & l_1 \\ l_c & l_d & l_2 \\ l_3 & l_4 & \Lambda \end{Bmatrix} C_{l_a0l_b0}^{l_10} C_{l_a0l_c0}^{l_30} C_{l_b0l_d0}^{l_40} C_{l_c0l_d0}^{l_20}.
```

According to Ref [^2], P251, Sec 8.5.2, Formula(32).
```math
C^{c0}_{a0b0} = \begin{cases}
0, & a + b + c = 2g + 1, \\
\dfrac{(-1)^{g-c}\hat{c}g!}{(g-a)!(g-b)!(g-c)!}\Delta(abc), & a +b +c = 2g.
\end{cases}
```
Define
```math
\Omega(abc) \equiv \dfrac{g!}{(g-a)!(g-b)!(g-c)!} = \frac{1}{(g+1)\Delta^2(\frac{a}{2}\frac{b}{2}\frac{c}{2})} = \binom{g}{c}\binom{c}{g-a}.
```
It is obvious a integer. Then
```math
X_{abcd} = (-1)^{f_x}\hat{l}_1\hat{l}_2\hat{l}_3\hat{l}_4\Omega(l_al_bl_1)\Delta(l_al_bl_1)\Omega(l_al_cl_3)\Delta(l_al_cl_3)\Omega(l_bl_dl_4)\Delta(l_bl_dl_4)\Omega(l_cl_dl_2)\Delta(l_cl_dl_2)\begin{Bmatrix}l_a & l_b & l_1 \\ l_c & l_d & l_2 \\ l_3 & l_4 & \Lambda \end{Bmatrix},
```
where
```math
\begin{aligned}
f_x &= -\frac{l_1+l_2+l_3+l_4}{2} + l_a+l_b+l_c+l_d + \dfrac{l_a+l_b-l_1}{2} + \dfrac{l_c+l_d-l_2}{2} + \dfrac{l_a+l_c-l_3}{2} + \dfrac{l_b+l_d-l_4}{2} \\
&= 2(l_a+l_b+l_c+l_d) - (l_1+l_2+l_3+l_4).
\end{aligned}
```
It is easy to prove that $l_1+l_2+l_3+l_4$ is even, so the phase factor $(-1)^{f_x}$ is $1$. In the Wigner symbols section, we have obtain
```math
\begin{Bmatrix}l_a & l_b & l_1 \\ l_c & l_d & l_2 \\ l_3 & l_4 & \Lambda \end{Bmatrix} = P_0\sum_t (-1)^{2t}(2t+1) \left(\sum_x A(t,x)\right) \left(\sum_y B(t,y)\right) \left(\sum_z C(t,z)\right).
```
Where, $A(t,x), B(t,y), C(t,z)$ are all integers, and
```math
P_0 = \Delta(l_al_bl_1)\Delta(l_cl_dl_2)\Delta(l_al_cl_3)\Delta(l_bl_dl_4)\Delta(l_1l_2\Lambda)\Delta(l_3l_4\Lambda).
```
So
```math
\begin{aligned}
X_{abcd} &= \hat{l}_1\hat{l}_2\hat{l}_3\hat{l}_4 \Omega(l_al_bl_1)\Omega(l_al_cl_3)\Omega(l_bl_dl_4)\Omega(l_cl_dl_2) \\
&\times \Delta^2(l_al_bl_1)\Delta^2(l_al_cl_3)\Delta^2(l_bl_dl_4)\Delta^2(l_cl_dl_2)\Delta(l_1l_2\Lambda)\Delta(l_3l_4\Lambda) \\
&\times\sum_t \hat{t}^2 \left(\sum_x A(t,x)\right) \left(\sum_y B(t,y)\right) \left(\sum_z C(t,z)\right).
\end{aligned}
```
In this formula, the only sqrt occurs in
```math
\hat{l}_1\hat{l}_2\hat{l}_3\hat{l}_4\Delta(l_1l_2\Lambda)\Delta(l_3l_4\Lambda),
```
which is only related to $l_1,l_2,l_3,l_4,\Lambda$, so can be taken out of the summation. If we only deal with $\tan\beta = 1$,
that's $(\sin\beta)^{2n_a + l_a + 2n_d + l_d}(\cos\beta)^{2n_b+l_b + 2n_c + l_c} = 2^{-\chi/2}$,
then the Moshinsky braket is a square root of a rational.

Define
```math
\begin{aligned}
M_{9j}(l_al_bl_cl_dl_1l_2l_3l_4\Lambda) &= \Omega(l_al_bl_1)\Omega(l_al_cl_3)\Omega(l_bl_dl_4)\Omega(l_cl_dl_2) \\
&\times \Delta^2(l_al_bl_1)\Delta^2(l_al_cl_3)\Delta^2(l_bl_dl_4)\Delta^2(l_cl_dl_2) \\
&\times \sum_t \hat{t}^2 \left(\sum_x A(t,x)\right) \left(\sum_y B(t,y)\right) \left(\sum_z C(t,z)\right).
\end{aligned}
```
With some simplification, we can write down the Moshinsky as follows:
```math
\begin{aligned}
&\langle n_3l_3,n_4l_4;\Lambda|n_1l_1,n_2l_2;\Lambda\rangle\\
= &\sqrt{\frac{\binom{\chi+2}{e_1+1}}{\binom{\chi+2}{e_3+1}\binom{l_1+l_2+\Lambda+1}{2\Lambda+1}\binom{2\Lambda}{\Lambda+l_1-l_2}\binom{l_3+l_4+\Lambda+1}{2\Lambda+1}\binom{2\Lambda}{\Lambda + l_3-l_4}}} \\
&\times\frac{1}{(e_1+2)(e_2+2)(2\Lambda+1)}\prod_{k=1}^4 \sqrt{\frac{(2l_k+1)\binom{e_k+l_k+1}{n_k+l_k}}{2^{l_k}\binom{e_k+1}{n_k}}} \\
&\times \sum_{n_an_bn_cn_dl_a l_b l_c l_d} \Bigg[(-1)^{l_d} (\sin\beta)^{e_a+e_d}(\cos\beta)^{e_b+e_c} M_{9j}(l_al_bl_cl_dl_1l_2l_3l_4\Lambda) \\
&\times \binom{e_1+2}{e_a+1}\binom{e_2+2}{e_c+1}  \prod_{p=a}^d 2^{l_p} (2l_p+1) \frac{\binom{e_p+1}{n_p}}{\binom{e_p+l_p+1}{n_p+l_p}}\Bigg]
\end{aligned}
```

Finally, if $\tan\beta = (m_1\omega_1/m_2\omega_2)^{1/2} = \sqrt{D}$ is a square root of rational, then
```math
\begin{aligned}
    (\sin\beta)^{e_a+e_d}(\cos\beta)^{e_b+e_c} &= (\sin\beta)^{e_4-e_1}(\cos\beta)^{e_1+e_3}(\tan\beta)^{2e_a} \\
    &= \frac{\sqrt{D}^{e_4-e_1}}{\sqrt{1+D}^{e_3+e_4}} D^{e_a},
\end{aligned}
```
where the $\frac{\sqrt{R}^{e_4-e_1}}{\sqrt{1+R}^{e_3+e_4}}$ can be taken aout from the summation. Espicially, when $D = 1$, we have
```math
    (\sin\beta)^{e_a+e_d}(\cos\beta)^{e_b+e_c} = 2^{-\frac{e_1+e_2}{2}}.
```

In summary, we prove that if $D = (m_1\omega_1)/(m_2\omega_2)$ is a rational, then the Moshinsky bracket can be expressed as a square root of a rational.

---

[^1]: Buck et al. Nuc. Phys. A 600 (1996) 387-402.
[^2]: A. N. Moskalev D. A. Varshalovich and V. K. Khersonskii, *Quantum theory of angular momentum*.