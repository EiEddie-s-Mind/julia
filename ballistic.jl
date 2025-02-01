module ballistic
export env, default, θ, t


struct Env
    g::Float64
    μ::Float64
    θ_range::Tuple{Float64,Float64}
end

global env_::Union{Nothing,Env} = nothing


"""
    env(g::Float64, μ::Float64, θ_range::Tuple{Float64,Float64})

为后续计算设置一个物理环境.

`θ_range` 是射角的范围, 闭区间, 并且应当是 ``(-π/2, π/2)`` 的一个子集.
注意端点取值不能是 ``±π/2``. 可以取一个接近的值, 例如 ``π/2 - 10^{-2}`` 或 ``-π/2 + 10^{-2}``.

`μ` 是阻力系数与质量的比值. 模型假设空气阻力 ``f`` 满足 ``f = - k v``, 则 ``μ = k/m``.

`g` 是重力加速度.
"""
env(g, μ, θ_range) = (global env_ = Env(g, μ, θ_range); nothing)


"""
    default()

默认的环境参数.
- `g` = ``9.8 \\mathrm{m/s^2}``
- `μ` = ``0.1 \\mathrm{s^{-1}}``
- `θ_range` = ``(-π/6, π/3)``

参见 [`env`](@ref).
"""
default() = env(9.8, 0.1, (-pi / 6, pi / 3))


"""
    θ(r::Vector, v0, g, μ, θ_range::Tuple)

计算以初速率 `v0` 击中位置 `r` 的射角.
`r` 应当是个二维向量.

考虑了空气阻力, 假设阻力 ``f`` 与速度 ``v`` 成正比, 方向相反;
即满足 ``f = - k v``.

关于 `g`, `μ` 与 `θ_range`, 参见 [`env`](@ref).
"""
function θ(r, v0, g, μ, θ_range)
    ϵ = 1e-3
    θ_lower, θ_upper = θ_range
    rx, ry = r

    # 目标在后方
    rx <= 0 && return missing
    # 目标太高, 超过最大倾角
    atan(ry, rx) > θ_upper && return missing
    # 使用无阻力抛体模型 (抛物线) 的包络线测试
    # 在包络线之下的总能击中
    # 因实际模型有阻力, 这只是必要不充分条件
    ry > -g / (2v0^2) * rx^2 + v0^2 / (2g) && return missing
    # 保证接下来的计算, 根号中不会出现负数
    v0 > g / μ && rx > v0 / sqrt(μ^2 - (g / v0)^2) && return missing

    # 目标函数
    f(ψ) = ry - g * (μ * rx * sec(ψ) - v0 * log(v0 / (v0 - μ * rx * sec(ψ)))) / (v0 * μ^2) - rx * tan(ψ)
    # 目标函数的导数
    f′(ψ) = rx * sec(ψ)^2 * (g * rx * sin(ψ) / (v0^2 * cos(ψ) - μ * v0 * rx) - 1.0)
    # 使得目标函数取最小值的 θ
    θ_st_min = atan(v0^2, g * rx) - atan(rx * μ, sqrt(v0^2 + rx^2 * ((g / v0)^2 - μ^2)))

    # θ s.t. min 小于最小倾角, 则使用曲线右边的较大的 θ 值
    # θ s.t. min 大于最小倾角, 则使用曲线右边的较小的 θ 值 (默认)
    θ = θ_st_min < θ_lower ? θ_upper : θ_lower
    fθ = f(θ)
    fθ_m = f(θ_st_min)

    # 最佳初始值
    θ = (θ * fθ_m - θ_st_min * fθ) / (fθ_m - fθ)
    # 5 次牛顿迭代
    for _ in 1:5
        θ -= f(θ) / f′(θ)
    end

    !(θ_lower <= θ <= θ_upper) && return missing
    abs(f(θ)) > ϵ && return missing
    return θ
end

import StaticModules: @with
"""
    θ(r::Vector, v0)

使用设置的物理环境计算射角.

!!! warning
    必须先使用 [`env`](@ref) 或 [`default`](@ref) 初始化物理环境.
"""
function θ(r, v0)
    isnothing(env_) && throw(error("env has no inited"))
    @with env_ begin
        return θ(r, v0, g, μ, θ_range)
    end
end


"""
    t(r::Vector, v0, θ)

计算在设置的物理环境下, 以初速率 `v0`, 射角 `θ` 发射的弹丸到目标 `r` 的运动时间.
`r` 应当是二维向量.

!!! note
    计算的是从起始点到 ``x`` 坐标与目标相等时的时间. ``y`` 可能仍有偏差.

!!! warning
    必须先使用 [`env`](@ref) 或 [`default`](@ref) 初始化物理环境.
"""
function t(r, v0, θ)
    isnothing(env_) && throw(error("env has no inited"))
    @with env_ begin
        return log(v0 / (v0 - μ * r[1] * sec(θ))) / μ
    end
end

import LinearAlgebra: norm
"""
    t(r::Vector, v::Vector)

计算在设置的物理环境下, 以初速度 `v` 发射的弹丸到目标 `r` 的运动时间.
`r`, `v` 应当是二维向量.
"""
t(r, v) = t(r, norm(v), atan(v[2], v[1]))


end
