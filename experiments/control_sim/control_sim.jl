using CSV, DataFrames, LinearAlgebra, Plots
using DifferentialEquations
using Optim


function load_data(filepath="data_m1.csv")
    data = CSV.File(filepath) |> DataFrame

    flexed_forces = [PointR3(row.P1dir1, row.P1dir2, row.P1dir3) for row in eachrow(data) if !any(ismissing, [row.P1dir1, row.P1dir2, row.P1dir3])]
    extended_forces = [PointR3(row.P2dir1, row.P2dir2, row.P2dir3) for row in eachrow(data) if !any(ismissing, [row.P2dir1, row.P2dir2, row.P2dir3])]

    return flexed_forces, extended_foerces
end


struct PointR3
    x::Float64
    y::Float64
    z::Float64
end

#  Σ: ẍ = f(x, ẋ, u, t)

mutable struct ControlSystem
    m::Float64
    beta::Float64
    k::Float64
    sigma::Float64
end

function system_init()
    return ControlSystem(0.1, 0.5, 10, 0.1)
end

# ∂x/∂t = f(x, u, t) + σ(u)ω, ω ∼ 𝒩(0, 1)
function system_dynamics!(du, u, p,t)

    system, tgt, ctrl = p

    # position
    x,y,z = u[1:3] # q

    # veloicty
    dx,dy,dz = u[4:6] # q̇

    # gaussian noise
    # control dependent st. ξ(u) = σ|u|ω
    gaussian_noise = system.noise * randn(3) .* abs.(ctrl)

    # mẍ + βẋ + κ(x - x*) = u + ξ

    du[1:3] = u[4:6]

    du[4] = (-system.beta * dx - system.k * (x - tgt.x) + ctrl[1] + gaussian_noise[1]) / system.m
    du[5] = (-system.beta * dy - system.k * (y - tgt.y) + ctrl[2] + gaussian_noise[2]) / system.m
    du[6] = (-system.beta * dz - system.k * (z - tgt.z) + ctrl[3] + gaussian_noise[3]) / system.m

end

# cost func (J(x,u)), optim (argmin(J)), etc. to be implemented
