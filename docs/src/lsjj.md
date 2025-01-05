# LS-jj recoupling

By defination
```math
\begin{aligned}
\mathrm{lsjj}(l_1,j_1,l_2,j_2,L,S,J) &= \langle l_1l_2L,s_1s_2S,J|l_1s_1j_1,l_2s_2j_2,J\rangle \\
&= \begin{bmatrix}l_1 & s_1 & j_1 \\ l_2 & s_2 & j_2 \\ L & S & J\end{bmatrix} \\
&= \sqrt{(2j_1+1)(2j_2+1)(2L+1)(2S+1)}\begin{Bmatrix}l_1 & s_1 & j_1 \\ l_2 & s_2 & j_2 \\ L & S & J\end{Bmatrix}.
\end{aligned}
```
We only deal with $s_1 = s_2 = 1/2$ condition.

## $S = 0$

Ref [^1] P357, 10.9.1, (1).
```math
\begin{Bmatrix}l_1 & s_1 & j_1 \\ l_2 & s_2 & j_2 \\ J & 0 & J\end{Bmatrix} = \dfrac{(-1)^{l_1+s_1+j_2+J}}{\sqrt{2(2J+1)}} \begin{Bmatrix}l_1 & s_1 & j_1 \\ j_2 & J & l_2\end{Bmatrix}.
```
So
```math
\begin{bmatrix}l_1 & s_1 & j_1 \\ l_2 & s_2 & j_2 \\ J & 0 & J\end{bmatrix} = (-1)^{l_1+s_1+j_2+J}\sqrt{\dfrac{(2j_1+1)(2j_2+1)}{2}} \begin{Bmatrix}l_1 & s_1 & j_1 \\ j_2 & J & l_2\end{Bmatrix}.
```

Here we define
```math
H_J(l_1,l_2,j_1,j_2) \equiv (-1)^{l_1+\frac{1}{2}+j_2+J}\sqrt{\dfrac{(2j_1+1)(2j_2+1)}{2}} \begin{Bmatrix}l_1 & \frac{1}{2} & j_1 \\ j_2 & J & l_2\end{Bmatrix}.
```
It is useful in laster formula.

## $S = 1$

### $L = J$

Ref [^1] P358, 10.9.2, (6).
```math
\begin{Bmatrix}l_1 & s_1 & j_1 \\ l_2 & s_2 & j_2 \\ J & 1 & J\end{Bmatrix} = (-1)^{l_1+s_1+j_2+J}\frac{(j_1-j_2)(j_1+j_2+1) - (l_1-l_2)(l_1+l_2+1)}{\sqrt{6J(J+1)(2J+1)}}\begin{Bmatrix}l_1 & \frac{1}{2} & j_1 \\ j_2 & J & l_2\end{Bmatrix}.
```
So
```math
\begin{bmatrix}l_1 & s_1 & j_1 \\ l_2 & s_2 & j_2 \\ J & 1 & J\end{bmatrix} = \frac{(j_1-j_2)(j_1+j_2+1) - (l_1-l_2)(l_1+l_2+1)}{\sqrt{J(J+1)}} H_J(l_1,l_2,j_1,j_2).
```

### $J = L + 1$
Ref [^1] P358, 10.9.2, (7).
```math
\begin{aligned}
&(-1)^{l_1+s_1+j_2+J-1}\sqrt{6J(2J-1)(2J+1)}\begin{Bmatrix}l_1 & s_1 & j_1 \\ l_2 & s_2 & j_2 \\ J-1 & 1 & J\end{Bmatrix} \\
=& \sqrt{(j_1-j_2+J)(j_2-j_1+J)(j_1+j_2+J+1)(j_1+j_2-J+1)}\begin{Bmatrix}l_1 & \frac{1}{2} & j_1 \\ j_2 & J-1 & l_2\end{Bmatrix} \\
 &+ \sqrt{(l_1-l_2+J)(l_2-l_1+J)(l_1+l_2+J+1)(l_1+l_2-J+1)} \begin{Bmatrix}l_1 & \frac{1}{2} & j_1 \\ j_2 & J & l_2\end{Bmatrix}.
\end{aligned}
```
So
```math
\begin{aligned}
\begin{bmatrix}l_1 & s_1 & j_1 \\ l_2 & s_2 & j_2 \\ J-1 & 1 & J\end{bmatrix} &= \sqrt{\frac{(j_1-j_2+J)(j_2-j_1+J)(j_1+j_2+J+1)(j_1+j_2-J+1)}{J(2J+1)}}H_{J-1}(l_1,l_2,j_1,j_2) \\
&-\sqrt{\frac{(l_1-l_2+J)(l_2-l_1+J)(l_1+l_2+J+1)(l_1+l_2-J+1)}{J(2J+1)}} H_J(l_1,l_2,j_1,j_2).
\end{aligned}
```

### $J = L - 1$
Ref [^1] P358, 10.9.2, (7).
```math
\begin{aligned}
&(-1)^{l_1+s_1+j_2+J}\sqrt{6J(2J+1)(2J+3)}\begin{Bmatrix}l_1 & s_1 & j_1 \\ l_2 & s_2 & j_2 \\ J-1 & 1 & J\end{Bmatrix} \\
=&\sqrt{(j_1-j_2+L)(j_2-j_1+L)(j_1+j_2+L+1)(j_1+j_2-L+1)} \begin{Bmatrix}l_1 & \frac{1}{2} & j_1 \\ j_2 & J+1 & l_2\end{Bmatrix} \\
&+\sqrt{(l_1-l_2+L)(l_2-l_1+L)(l_1+l_2+L+1)(l_1+l_2-L+1)} \begin{Bmatrix}l_1 & \frac{1}{2} & j_1 \\ j_2 & J & l_2\end{Bmatrix}.
\end{aligned}
```
So
```math
\begin{aligned}
\begin{bmatrix}l_1 & s_1 & j_1 \\ l_2 & s_2 & j_2 \\ J+1 & 1 & J\end{bmatrix} 
&=\sqrt{\frac{(l_1-l_2+L)(l_2-l_1+L)(l_1+l_2+L+1)(l_1+l_2-L+1)}{(J+1)(2J+1)}} H_J(l_1,l_2,j_1,j_2) \\
&- \sqrt{\frac{(j_1-j_2+L)(j_2-j_1+L)(j_1+j_2+L+1)(j_1+j_2-L+1)}{(J+1)(2J+1)}}H_{J+1}(l_1,l_2,j_1,j_2).
\end{aligned}
```

## $H$ function

According to the symmetry of Wigner-6j symbol, one can find
```math
H_J(l_1,l_2,j_1,j_2) = (-1)^{l_1+l_2+j_1+j_2+1} H_{J}(l_2,l_1,j_2,j_1).
```

- if $j_1 = l_1 + \frac12,\; j_2 = l_2 + \frac12$

Ref [^1] P300, 9.5.2, (5).
```math
H_J(l_1,l_2,j_1,j_2) = \sqrt{\dfrac{(j_1+j_2+J+1)(j_1+j_2-J)}{8j_1j_2}}.
```

- if $j_1 = l_1 + \frac12,\; j_2 = l_2 - \frac12$

Ref [^1] P300, 9.5.2, (4).
```math
H_J(l_1,l_2,j_1,j_2) = \sqrt{\dfrac{(j_1-j_2+J)(j_2-j_1+J+1)}{8j_1(j_2+1)}}.
```

- if $j_1 = l_1 - \frac12,\; j_2 = l_2 + \frac12$
```math
H_J(l_1,l_2,j_1,j_2) = -\sqrt{\dfrac{(j_1-j_2+J+1)(j_2-j_1+J)}{8(j_1+1)j_2}}.
```

- if $j_1 = l_1 - \frac12,\; j_2 = l_2 - \frac12$
```math
H_J(l_1,l_2,j_1,j_2) = \sqrt{\dfrac{(j_1+j_2+J+2)(j_1+j_2-J+1)}{8(j_1+1)(j_2+1)}}.
```

## Result

Define
```math
p = l_1+l_2, \quad, m = l_1 - l_2
```

- if $j_1 = l_1 + \frac12,\; j_2 = l_2 + \frac12$
```math
\begin{aligned}
(S=0, L = J) &\to \sqrt{\frac{(J+p+2) (-J+p+1)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}, \\
(S=1, L = J) &\to \frac{m}{\sqrt{J(J+1)}}\sqrt{\frac{(J+p+2) (-J+p+1)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}, \\
(S=1, L = J-1)&\to \sqrt{\frac{J^2-m^2}{J (2 J+1)}} \sqrt{\frac{(J+p+1) (J+p+2)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}, \\
(S=1, L = J+1)&\to -\sqrt{\frac{(J+1)^2-m^2}{(J+1) (2 J+1)}} \sqrt{\frac{(p-J) (-J+p+1)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}.
\end{aligned}
```

- if $j_1 = l_1 + \frac12,\; j_2 = l_2 - \frac12$
```math
\begin{aligned}
(S=0, L = J) &\to \sqrt{\frac{(J+m+1) (J-m)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}, \\
(S=1, L = J) &\to \frac{p+1}{\sqrt{J(J+1)}}\sqrt{\frac{(J+m+1) (J-m)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}, \\
(S=1, L = J-1)&\to -\sqrt{\frac{(p+1)^2-J^2}{J (2 J+1)}} \sqrt{\frac{(J+m+1) (J+m)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}, \\
(S=1, L = J+1)&\to -\sqrt{\frac{(p+1)^2-(J+1)^2}{(J+1) (2 J+1)}} \sqrt{\frac{(J-m+1) (J-m)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}.
\end{aligned}
```

- if $j_1 = l_1 - \frac12,\; j_2 = l_2 + \frac12$
```math
\begin{aligned}
(S=0, L = J) &\to -\sqrt{\frac{(J+m) (J-m+1)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}, \\
(S=1, L = J) &\to \frac{p+1}{\sqrt{J(J+1)}}\sqrt{\frac{(J+m) (J-m+1)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}, \\
(S=1, L = J-1)&\to \sqrt{\frac{(p+1)^2-J^2}{J (2 J+1)}} \sqrt{\frac{(J-m) (J-m+1)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}, \\
(S=1, L = J+1)&\to \sqrt{\frac{(p+1)^2-(J+1)^2}{(J+1) (2 J+1)}} \sqrt{\frac{(J+m+1) (J+m)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}.
\end{aligned}
```

- if $j_1 = l_1 - \frac12,\; j_2 = l_2 - \frac12$
```math
\begin{aligned}
(S=0, L = J) &\to \sqrt{\frac{(J+p+1) (p-J)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}, \\
(S=1, L = J) &\to \frac{-m}{\sqrt{J(J+1)}}\sqrt{\frac{(J+p+1) (p-J)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}, \\
(S=1, L = J-1)&\to -\sqrt{\frac{J^2-m^2}{J (2 J+1)}} \sqrt{\frac{(p-J) (-J+p+1)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}, \\
(S=1, L = J+1)&\to \sqrt{\frac{(J+1)^2-m^2}{(J+1) (2 J+1)}} \sqrt{\frac{(J+p+1) (J+p+2)}{2 (2 \text{l1}+1) (2 \text{l2}+1)}}.
\end{aligned}
```

[^1]: A. N. Moskalev D. A. Varshalovich and V. K. Khersonskii, *Quantum theory of angular momentum*.