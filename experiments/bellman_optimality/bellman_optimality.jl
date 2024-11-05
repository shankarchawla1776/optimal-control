using Base: ExtensionId
using LinearAlgebra
using DataFrames
using CSV
using Plots
using StaticArrays

# Load data from CSV file
df = CSV.File("data_m1.csv") |> DataFrame

# Define a function to handle missing values or unequal data lengths
function safe_svector(x, y, z)
    if any(ismissing.([x, y, z]))
        return nothing  # Skip this vector if any value is missing
    else
        return SVector(x, y, z)
    end
end

# Define flexed and extended vectors, filtering out rows with incomplete data
flexed_vecs = [safe_svector(df.P1dir1[i], df.P1dir2[i], df.P1dir3[i]) for i in 1:size(df,1) if safe_svector(df.P1dir1[i], df.P1dir2[i], df.P1dir3[i]) !== nothing]
extended_vecs = [safe_svector(df.P2dir1[i], df.P2dir2[i], df.P2dir3[i]) for i in 1:size(df,1) if safe_svector(df.P2dir1[i], df.P2dir2[i], df.P2dir3[i]) !== nothing]

# Ensure flexed_vecs and extended_vecs have the same length
n = min(length(flexed_vecs), length(extended_vecs))
flexed_vecs = flexed_vecs[1:n]
extended_vecs = extended_vecs[1:n]

target = normalize(SVector(1,1,1))

# the cost function l(x,u,t)
# a state vector x has its deviation
# from a target vector
function l(x, tgt)
    # find the normed difference
    return norm(x - tgt)
end

flexed_costs = [l(vector, target) for vector in flexed_vecs]
extended_costs = [l(vector, target) for vector in extended_vecs]

# backwards time iteration
# to minimize the cost-to-go
# we compute V(t - Δt, x) from V(t, x)
function backwards_iteration(vecs, tgt, iters=10)
    # the value function V(t,x)
    V = [l(x, tgt) for x in vecs]

    for _ in 1:iters
        V_new = similar(V)
        for i in 1:length(vecs)
            # define an update rule
            V_new[i] = min(V[i], V[i] - 0.1 * (V[i] - l(vecs[i], tgt)))
        end
        # the .= syntax is analogous to Δ=
        V .= V_new
    end
    return V
end

# do value iteration on both postures

function H(x,tgt,lambda)
    # f(x,u) is shown as the dot
    # product of x and the target

    # initialize the instantaneous
    # cost l(x,u,t) -> l(x,u)
    cost = l(x,tgt)
    dynamics = dot(x,tgt)
    return cost + lambda * dynamics # this is the full hamiltonian form
end

opt_ctrl_flexed = [flexed_vecs[argmin([H(x, target, 1) for x in flexed_vecs])]]
opt_ctrl_extended = [extended_vecs[argmin([H(x, target, 1) for x in extended_vecs])]]

plot()
scatter!(flexed_vecs, label="Posture 1 State Vectors (x_P1)", legend=:topright)
scatter!(extended_vecs, label="Posture 2 State Vectors (x_P2)")
plot!([target], label="Target Vector", lw=3, marker=:circle)
title!("3D Force Vectors and Target Vector Alignment")
savefig("bellman_plot.png")
println("Initial costs for Posture 1 (V_P1):", flexed_costs)
println("Initial costs for Posture 2 (V_P2):", extended_costs)
println("Optimal controls for P1 under Hamiltonian minimization (u*):", opt_ctrl_flexed)
println("Optimal controls for P2 under Hamiltonian minimization (u*):", opt_ctrl_extended)
