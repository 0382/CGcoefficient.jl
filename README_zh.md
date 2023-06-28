# CGcoefficient.jl

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE)
[![CI](https://github.com/0382/CGcoefficient.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/0382/CGcoefficient.jl/actions/workflows/CI.yml)
[![codecov.io](http://codecov.io/github/0382/CGcoefficient.jl/coverage.svg?branch=master)](http://codecov.io/github/0382/CGcoefficient.jl?branch=master)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://0382.github.io/CGcoefficient.jl/dev)

[[English](README.md)]

计算CG系数，Racah系数和Wigner 3j, 6j, 9j系数。

你可以使用定义的`SqrtRational`类型以保存准确结果，它内部使用`BigInt`计算以避免溢出。同时我们也提供浮点数运算的快速版本，比[GNU Scientific Library](https://www.gnu.org/software/gsl/)快一倍。

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

### 参考

- [https://github.com/ManyBodyPhysics/CENS](https://github.com/ManyBodyPhysics/CENS)
- D. A. Varshalovich, A. N. Moskalev and V. K. Khersonskii, *Quantum Theory of Angular Momentum*, (World Scientific, 1988).