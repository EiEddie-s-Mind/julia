"""
	∨: so(3) -> R^3
"""
∨(Ω::SMatrix{3,3}) = @SVector [Ω[3, 2], Ω[1, 3], Ω[2, 1]]

"""
	∧: R^3 -> so(3)
"""
∧(ω::SVector{3}) = @SMatrix [0 -ω[3] ω[2]; ω[3] 0 -ω[1]; -ω[2] ω[1] 0]

"""
	log: SO(3) -> so(3)
"""
function Base.log(R::SMatrix{3,3})
	θ = acos((tr(R) - 1) / 2)
	# TODO: what is returned when R=I
	θ ≈ 0 && return @SMatrix zeros(3,3)
	# TODO: another better way
	# 可以证明, 当 θ = π 时, dim Ker(R-I) = 1
	θ ≈ pi && return pi * ∧(SVector{3}(nullspace(R-I)[:,1]))
	W = (R - R') / 2sin(θ)
	return θ * W
end

"""
	exp: so(3) -> SO(3)
"""
function Base.exp(Ω::SMatrix{3,3})
	θ = sqrt(-tr(Ω^2)/2)
	θ ≈ 0. && return SMatrix{3,3}(I)
	Ω = Ω / θ
	return I + sin(θ)*Ω + (1-cos(θ))*Ω^2
end

"""
	Log: SO(3) -> R^3
"""
Log(R::SMatrix{3,3}) = ∨(log(R))

rot(θ, ω::SVector{3}) = cos(θ)*I + (1 - cos(θ))*ω*ω' + sin(θ)*∧(ω)

"""
	Exp: R^3 -> SO(3)
"""
#Exp(v::Vector) = exp(∧(v))
Exp(v::SVector{3}) = norm(v) ≈ 0. ? SMatrix{3,3}(I) : rot(norm(v), normalize(v))

"""
	Ad: SO(3) -> R^3×3

将流形上的元素 R 处的切空间内的向量线性变换到另一元素处的切空间内.
即对于坐标系 ``{a}``, 其李群 ``R_a`` 切空间内表示的李代数的向量形式 ``ω_a``,
有另一坐标系 ``{b}`` 与李群 ``R_b``, 与坐标系变换 ``R_{ba}``;
那么 ``R_b`` 切空间内表示的向量 ``ω_b = Ad(R_{ba}) ω_a``.
"""
Ad(R::SMatrix{3,3}) = R
