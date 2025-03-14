module kinetics
export log, Log, exp, Exp, ∨, ∧
using LinearAlgebra
# 使用定长向量与矩阵实现对 SO(3), SE(3) 元素的多重派发
using StaticArrays

include("rotation.jl")

export ⊕, unpack, normalize
include("movement.jl")

end
