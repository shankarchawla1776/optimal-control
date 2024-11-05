# Implementation of a hamiltonian system for numerical optimal control
# does not use empirical force vector data

using Plots
using DifferentialEquations

# define system parameters
a = 1
b = 1

# initial conditioin for the system state
x0 = 5

# time horizon & time step
T = 10
dt = 0.1

# control cost weight
R = 0.1

# the hamiltonian dynamics (for the state and costate)
function hamiltonian_dynamics!(du, u, t, p)
    # defint the state & costate
    x, p = u[1], u[2]

    # the optimal control law u*
    opt_ctr_law = -p * b / R

    # equations for the state and costaet
    # dx/dt
    du[1] = a * x + b * opt_ctr_law
    # dp/dt
    du[2] = -x - a * p
end


# initial and boundary conditions
u0 = [x0, 0]
t_int = (0, T) # the full time interval

# use Julia's ODE solver
# to solve the differential equation
ODE = ODEProblem(hamiltonian_dynamics!, u0, t_int)
println("here")
solution = solve(ODE, Tsit5(), reltol=1e-8, abstol=1e-8)

# find the optimal control trajectory
u_optimal = [-solution[i][2] * b / R for i in 1:length(solution.t)]
x_traj = [solution[i][1] for i in 1:length(solution.t)]
p_traj = [solution[i][2] for i in 1:length(solution.t)]


# Step 6: Plot results
plot(solution.t, x_traj, label="State x(t)", xlabel="Time (t)", ylabel="State / Costate", lw=2)
plot!(solution.t, x_traj, label="Costate p(t)", lw=2)
plot!(solution.t, u_optimal, label="Control u(t)", lw=2, ylabel="Control u(t)")
title!("Optimal Control of Linear System")

# Save the plot as an image file
savefig("optimal_control_plot.png")  # Saves the plot as 'optimal_control_plot.png' in the current directory


# find J's cost

terminal_J = 0.5 * x_traj[end]^2
integral = sum(0.5 * x^2 + 0.5 * R * u^2 for (x, u) in zip(x_traj, u_optimal)) * dt

cost = terminal_J * integral
println(cost)
