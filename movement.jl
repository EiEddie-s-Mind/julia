"""
	R^3 × R^3 -> R^6

将轴角表示的旋转与平移粘合为 6 维向量, 也称运动旋量.
"""
⊕(ω::SVector{3}, v::SVector{3}) = [ω; v]

"""
	SO(3) × R^3 -> SE(3)

直接由旋转 R 与平移 ρ 生成刚体运动.
"""
⊕(R::SMatrix{3,3}, ρ::SVector{3}) = SMatrix{4,4}([R ρ; zeros(1,3) 1])

"""
	R^6 -> R^3 × R^3

将李代数对应的 6 维向量拆为 轴角表示的旋转 与 平移 两个 3 维向量.
"""
unpack(ξ::SVector{6}) = (SVector{3}(ξ[1:3]), SVector{3}(ξ[4:6]))

"""
	se(3) -> so(3) × R^3
	SE(3) -> SO(3) × R^3

将李代数 se(3) 或李群 SE(3) 拆为 李代数 so(3) 或李群 SO(3) 与 平移.
"""
unpack(X::SMatrix{4,4}) = (SMatrix{3,3}(X[1:3,1:3]), SVector{3}(X[1:3,4]))

import LinearAlgebra.normalize
normalize(ξ::SVector{6}) = ((ω,v)=unpack(ξ); l=norm(ω); l ≈ 0 ? ξ/norm(v) : ξ/l)


"""
	∧: R^6 -> se(3)
"""
∧(ξ::SVector{6}) = ((ω,v)=unpack(ξ); SMatrix{4,4}([∧(ω) v; zeros(1,4)]))

"""
	∨: se(3) -> R^6
"""
∨(X::SMatrix{4,4}) = ((Ω,v)=unpack(X); ∨(Ω) ⊕ v)

"""
	log: SE(3) -> se(3)
"""
function Base.log(T::SMatrix{4,4})
	(R, ρ) = unpack(T)
	Ω = log(R)
	θ = sqrt(-tr(Ω^2) / 2)
	# 注意 Ω 不是 *单位* 向量的外积矩阵
	# Ω = θ [ω̂]
	J′ = θ ≈ 0 ? I : I - 1/2*Ω + (1/θ^2-cot(θ/2)/2θ)*Ω^2
	return SMatrix{4,4}([Ω J′*ρ; zeros(1,4)])
end

"""
	exp: se(3) -> SE(3)
"""
function Base.exp(X::SMatrix{4,4})
	(Ω, v) = unpack(X)
	θ = sqrt(-tr(Ω^2)/2)
	R = exp(Ω)
	# 注意 Ω 不是 *单位* 向量的外积矩阵
	# Ω = θ [ω̂]
	Ω /= θ
	J = θ ≈ 0 ? I : I + (1-cos(θ))/θ*Ω + (1-sin(θ)/θ)*Ω^2
	return SMatrix{4,4}([R J*v; zeros(1,3) 1])
end

"""
	Log: SE(3) -> R^6
"""
Log(T::SMatrix{4,4}) = ∨(log(T))

"""
	Exp: R^6 -> SE(3)
"""
Exp(ξ::SVector{6}) = exp(∧(ξ))
