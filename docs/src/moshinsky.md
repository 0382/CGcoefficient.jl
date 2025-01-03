# Moshinsky

The Talmi-Moshinsky braket can be evaluated as following [^1]：
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
    2n_a + l_a + 2n_b + l_b = 2n_1 + l_1 \equiv f_1, \\
    2n_a + l_a + 2n_c + l_c = 2n_3 + l_3 \equiv f_3, \\
    2n_b + l_b + 2n_d + l_d = 2n_4 + l_4 \equiv f_4, \\
    2n_c + l_c + 2n_d + l_d = 2n_2 + l_2 \equiv f_2.
\end{aligned}
```

Here we define $\chi = f_1 + f_2 = f_3 + f_4$。

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
    &= 2^{-\frac{f_1+f_2}{2}} \\
    &= 2^{-\chi/2}
\end{aligned}
```

[^1]: Buck et al. Nuc. Phys. A 600 (1996) 387-402.