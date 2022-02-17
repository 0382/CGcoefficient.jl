# CGcoefficient.jl

[[English](README.md)]

计算CG系数，Racah系数和Wigner 3j, 6j, 9j系数，并使用定义的`SqrtRational`类型以保存准确结果。

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
4√(1//24)

julia> simplify(ans)
√(2//3)

julia> simplify(nineJ(1,2,3,4,5,6,3,6,9))
1//1274√(3//5)
```

更多示例请看[文档](https://0382.github.io/CGcoefficient.jl-docs/)。

### 参考

- [https://github.com/ManyBodyPhysics/CENS](https://github.com/ManyBodyPhysics/CENS)
- D. A. Varshalovich, A. N. Moskalev and V. K. Khersonskii, Quantum Theory of Angular Momentum, (World Scientific, 1988).