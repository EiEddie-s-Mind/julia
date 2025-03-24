module series


_δ = 1e-2


"""
	depth(::Type{<:AbstractFloat}, δ::Float) -> Integer

计算浮点类型按精度 `δ` 求导的可行深度.
此值也是求导阶数 + 1.
"""
depth(::Type{Float}, δ) where {Float<:AbstractFloat} = floor(Int, log(eps(Float)) / log(δ)) - 1


mutable struct Series{Float<:AbstractFloat}
	const depth::Int
	coef::Vector{Float}
end

"""
	Series(seq::Vector, depth::Integer, δ::Float=_δ)

根据 `seq` 内的值, 按精度 `δ` 初始化求导列表.
若 `seq` 的长度大于 `depth`, 将只取后 `depth` 个值.
"""
function Series(seq::Vector{Float}, depth, δ::Float=_δ) where {Float<:AbstractFloat}
	depth = min(depth, length(seq))
	coef = seq[end-depth+1:end]
	for i in depth-1:-1:1
		for j in 1:i
			coef[j] = (coef[j+1] - coef[j]) / δ
		end
	end
	reverse!(coef)
	return Series(depth, coef)
end

"""
	Series(seq::Vector, δ::Float=_δ)

自动计算深度.
"""
Series(seq::Vector{Float}, δ::Float=_δ) where {Float<:AbstractFloat} = Series(seq, depth(Float,δ), δ)


"""
	update!(s::Series, x::Float, δ::Float=_δ) -> Vector

使用新的观测值 `x` 更新模型 `s`, 返回更新后模型的参数列表.
"""
function update!(s::Series{Float}, x::Float, δ::Float=_δ) where {Float<:AbstractFloat}
	it = x
	for i in 1:s.depth-1
		it, s.coef[i] = s.coef[i], it
		it = (s.coef[i] - it) / δ
	end
	s.coef[end] = it
	return s.coef
end


"""
	gen(coef::Vector) -> (Float -> Float)

根据参数列表 `coef` 返回对应的 `n` 阶泰勒展开.
可用于推算被预测函数未来的取值.
`n` 为模型的深度.

# Examples
```julia-repl
# exp 在 1 处各级导数都为 e
julia> coef = repeat([float(ℯ)], 6);
julia> f = gen(coef);

julia> f(0.1)
3.004166020116425

julia> exp(1.1)
3.0041660239464334
```
"""
function gen(coef::Vector{Float}) where {Float<:AbstractFloat}
	taylor = coef
	fact = 1.0
	for i in 1:length(coef)
		taylor[i] /= fact
		fact *= i
	end

	function f(x)
		sum = 0.0
		for (i, v) in enumerate(taylor)
			sum += v * x^(i - 1)
		end
		sum
	end
	return f
	# return x -> enumerate(taylor) .|> (((i, v),) -> v * x^(i - 1)) |> sum
end

"""
	gen(s::Series) -> (Float -> Float)

返回模型对应的预测函数. 此函数的变量为自此时刻之后的时间.

# Examples
```julia-repl
julia> s = Series(0:0.01:1 .|> exp);
julia> f = gen(s);

julia> f(0.1)
3.0026693046333848

julia> exp(1.1)
3.0041660239464334
```
"""
gen(s::Series) = gen(s.coef)


end
