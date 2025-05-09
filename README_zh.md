# CGcoefficient.jl

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE)
[![CI](https://github.com/0382/CGcoefficient.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/0382/CGcoefficient.jl/actions/workflows/CI.yml)
[![codecov.io](https://codecov.io/gh/0382/CGcoefficient.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/0382/CGcoefficient.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://0382.github.io/CGcoefficient.jl/dev)

[[English](README.md)]

计算CG系数，Racah系数，Wigner 3j, 6j, 9j系数，Moshinsky系数等。

你可以使用定义的`SqrtRational`类型以保存准确结果，它内部使用`BigInt`计算以避免溢出。同时我们也提供浮点数运算的快速版本，比[GNU Scientific Library](https://www.gnu.org/software/gsl/)快好几倍。

我也用c++重写了浮点数版本，用于数值计算：[WignerSymbol](https://github.com/0382/WignerSymbol).

关于更多更多细节以及计算公式，请看文档[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://0382.github.io/CGcoefficient.jl/dev)。

### 安装

使用julia REPL安装即可
```julia-repl
julia> ]
pkg> add CGcoefficient
```

### 示例

```julia-repl
julia> CG(1,2,3,1,1,2)
√(2//3)

julia> nineJ(1,2,3,4,5,6,3,6,9)
1//1274√(3//5)

julia> f6j(6,6,6,6,6,6)
-0.07142857142857142
```

更多示例请看文档。

### 函数列表

这个包提供两类函数。

1. 精确函数。返回`SqrtRational`；这类函数主要是为了演示，内部使用`BigInt`进行计算，不需要缓存，计算速度相对较慢。
2. 浮点数函数。返回双精度浮点数；这类函数主要是为了数值计算，内部使用`Float64`进行计算。在使用之前，你需要先调用[`wigner_init_float`](https://0382.github.io/CGcoefficient.jl/stable/api/#CGcoefficient.wigner_init_float)来预先计算二项式系数表并缓存这个表，用于后面的计算。当角动量量子数非常大的时候，它们可能会由于浮点数计算产生一些误差。对于误差的分析详见：[wigner-benchmark](https://github.com/0382/wigner-benchmark)。

#### 精确函数

精确函数仅用于演示，其中一些可以是`Real`参数，比如`1, 3//2, 0.5`，但是`4//3, 0.6`会报错。

- `CG(j1, j2, j3, m1, m2, m3)`, CG系数，参数是均为`Real`。
- `CG0(j1, j2, j3)`, CG系数特殊情况`m1 = m2 = m3 = 0`，此时当然角动量是整数才有意义。
- `threeJ(j1, j2, j3, m1, m2, m3)`, Wigner 3j系数，参数均为`Real`。
- `sixJ(j1, j2, j3, j4, j5, j6)`, Wigner 6j系数，参数均为`Real`。
- `Racah(j1, j2, j3, j4, j5, j6)`, Racah系数，参数均为`Real`。
- `nineJ(j1, j2, j3, j4, j5, j6, j7, j8, j9)`, Wigner 9j系数，参数均为`Real`。
- `norm9J(j2, j3, j4, j5, j5, j6, j7, j8, j9)`, normalized 9j系数，参数均为`Real`。
- `lsjj(l1, l2, j1, j2, L, S, J)`, LS耦合到jj耦合的转换系数，它实际上等于一个normalized 9j系数，但更易于使用且更快。`j1, j2`是`Real`，其余必须是整数。
- `Moshinsky(N, L, n, l, n1, l1, n2, l2, Λ, D)`, Moshinsky括号，量子数都是整数，`D`是有理数，默认为`1`。

#### 浮点数函数

对于数值计算而言，为了效率我们避免使用分数类型，如果某个参数可能是半整数，函数的则使用两倍的角动量量子数作为参数，以此避免半整数。在形参命名中，如果某个参数接受两倍的角动量，那么会有一个`d`前缀。

- `fCG(dj1, dj2, dj3, dm1, dm2, dm3)`, CG系数
- `fCG0(j1, j2, j3)`, CG系数特殊情况`m1 = m2 = m3 = 0`。
- `fCGspin(ds1, ds2, S)`, 快速计算两个1/2自旋的CG系数。
- `fCG3spin(ds1, ds2, ds3, S12, dS)`, 快速计算 `<S12,M12|1/2,m1;1/2,m2><S,M|S12,M12;1/2,m3>`.
- `f3j(dj1, dj2, dj3, dm1, dm2, dm3)`, Wigner 3j系数。
- `f6j(dj1, dj2, dj3, dj4, dj5, dj6)`, Wigner 6j系数。
- `fRacah(dj1, dj2, dj3, dj4, dj5, dj6)`, Racah系数。
- `f9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9)`, Wigner 9j系数。
- `fnorm9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9)`, normalized 9j系数。
- `flsjj(l1, l2, dj1, dj2, L, S, J)`, LS耦合到jj耦合的转换系数。
- `fMoshinsky(N, L, n, l, n1, l1, n2, l2, Λ)`, oshinsky括号。
- `dfunc(dj, dm1, dm2, β)`, Wigner d 函数。


### 参考资料

- [https://github.com/ManyBodyPhysics/CENS](https://github.com/ManyBodyPhysics/CENS)
- D. A. Varshalovich, A. N. Moskalev and V. K. Khersonskii, *Quantum Theory of Angular Momentum*, (World Scientific, 1988).
- Buck et al. Nuc. Phys. A 600 (1996) 387-402.
- H. T. Johansson and C. Forssén, Fast and Accurate Evaluation of Wigner 3$j$, 6$j$, and 9$j$ Symbols Using Prime Factorization and Multiword Integer Arithmetic, SIAM J. Sci. Comput. 38, A376 (2016).
