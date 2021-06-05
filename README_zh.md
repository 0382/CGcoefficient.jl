# CGcoefficient.jl

[[English](README.md)]

计算CG系数和Wigner 3j, 6j, 9j系数，并给出准确的结果。为了得到准确结果，所有的中间计算都是使用`BigInt`，并且定义了`SqrtRational`类型以保存结果。

关于更多更多细节以及计算公式，请看[文档](https://0382.github.io/CGcoefficient.jl-docs/)。

### 安装

使用julia REPL安装即可
```julia-repl
julia> ]
pkg> add CGcoefficient
```

### 示例

```julia-repl
julia> CG(1,2,3,1,1,2)
√2/√3
julia> nineJ(1,2,3,4,5,6,3,6,9)
√3/1274√5
```

更多示例请看[文档](https://0382.github.io/CGcoefficient.jl-docs/)。

### 参考

- [https://github.com/ManyBodyPhysics/CENS](https://github.com/ManyBodyPhysics/CENS)
- D. A. Varshalovich, A. N. Moskalev and V. K. Khersonskii, Quantum Theory of Angular Momentum, (World Scientific, 1988).