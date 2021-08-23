# Cart-pole control using Hermite-Simpson collocation
using JuMP, LinearAlgebra, Interpolations, Plots, JLD2
using Ipopt
# using NLopt
# include("cartPoleDynamics.jl")
# using .cartPoleDynamics

cd("/Users/blake/Dropbox/CMU/julia/controls_problems/pendulum_opt")

# Set initial conditions
const θ_s = 0.0
const θdot_s = 0
const u_s = 0
const t_s = 0.000
const t_F = 1.0

const θ_f = pi
const θdot_f = 0
const u_f = 0
const n = 50
const integration_rule = "trapezoidal"
const h_k = t_F / n

# Create JuMP model for NLopt
# model = Model(NLopt.Optimizer)
# set_optimizer_attribute(model, "algorithm", :LN_COBYLA)

# Create JuMP model for NLopt
model = Model(Ipopt.Optimizer)

@variables(model, begin
    -Inf ≤ θ[1:n] ≤ Inf
    -Inf ≤ θdot[1:n] ≤ Inf
    -10 ≤ u[1:n] ≤ 10
end)

# Fix initial conditions
fix(θ[1], θ_s, force=true)
fix(θdot[1], θdot_s, force=true)

# Fix final conditions
fix(θ[n], θ_f, force=true)
fix(θdot[n], θdot_f, force=true)

# Initial guess: linear interpolation between boundary conditions
x_s = [θ_s, θdot_s, u_s]
x_f = [θ_f, θdot_f, u_s]
interp_linear = Interpolations.LinearInterpolation([1, n], [x_s, x_f])
initial_guess = mapreduce(transpose, vcat, interp_linear.(1:n))
set_start_value.(all_variables(model), vec(initial_guess))

# # Get an initial guess
# JLD2.@load "data/guessSolution.jld2" solOutput
# initial_guess = solOutput
# set_start_value.(all_variables(model), vec(initial_guess))

# Define helper functions for discrete cartpole dynamics
# This is fₖ
@NLexpression(model, δθ[j = 1:n], θdot[j])
@NLexpression(model, δθdot[j = 1:n], u[j]-sin(θ[j]))
@NLexpression(model, δθdotmid[j = 2:n], u[j]-sin(θ[j]))

# This is f_{k+1}
# @NLexpression(model, δxck[j = 2:n], xcdot[j])
# @NLexpression(model, δθk[j = 2:n], θdot[j])
# @NLexpression(model, δxcdotk[j = 2:n], (1/(1+sin(θ[j])^2)*(u[j] + sin(θ[j])*(θdot[j]^2 + cos(θ[j])))))
# @NLexpression(model, δθdotk[j = 2:n], (1/(1+sin(θ[j])^2)*(-u[j]*cos(θ[j]) - θdot[j]^2 * cos(θ[j])*sin(θ[j]) - 2*sin(θ[j]))))

@NLexpression(model, θmid[j = 2:n-1], 0.5*(θ[j-1] - θ[j]) + (h_k/8)*(δθ[j] - δθ[j+1]))
@NLexpression(model, θdotmid[j = 2:n-1], 0.5*(θdot[j-1] - θdot[j]) + (h_k/8)*(δθdotmid[j] - δθdotmid[j+1]))

# Fix initial conditions
@NLconstraint(model, θ[1] - θ_s == 0)
@NLconstraint(model, θdot[1] - θdot_s == 0)

# Fix final conditions
@NLconstraint(model, θ[n] - θ_f == 0)
@NLconstraint(model, θdot[n] - θdot_f == 0)

# System dynamics 
for j in 2:n-1
    i = j - 1
    # Trapezoidal integration
    @NLconstraint(model, θ[j] == θ[i] + 0.5 * h_k * (δθ[j] + δθ[i]))
    @NLconstraint(model, θdot[j] == θdot[i] + 0.5 * h_k * (δθdot[j] + δθdot[i]))

    # Constraints for Hermite-Simpson collocation
    # @NLconstraint(model, θ[j] - θ[i] - (h_k/6)*(δθ[i] + 4*θmid[j] + δθ[j]) == 0)
    # @NLconstraint(model, θdot[j] - θdot[i] - (h_k/6)*(δθdot[i] + 4*θdotmid[j] + δθdot[j]) == 0)
end

# Objective is to minimize control effort
# @NLexpression(model, wk[1:n], u[1:n])
# @NLexpression(model, wkmid[j=1:n-1], u[j])
# @objective(model, Min, (h_k/6)*(u[1] + 4*wkmid[2] + u[2]) + sum((h_k/6)*(u[j-1] + 4*wkmid[j] + u[j]) for j in 3:n-1))
@objective(model, Min, dot(u,u))

# set_silent(model)
JuMP.optimize!(model)
# @assert termination_status(model) == MOI.LOCALLY_SOLVED

ts = [h_k*i for i in 1:n]
# plt_cartPos = plot(ts, value.(xc), title="Cart position")
plt_pendPos = plot(ts, value.(θ), title="Pendulum rotational position")
# plt_cartVel = plot(ts, value.(xcdot), title="Cart velocity")
plt_pendVel = plot(ts, value.(θdot), title="Pendulum rotational velocity")
plt_control = plot(ts, value.(u), title="Control input")
# plt = plot(plt_pendPos, plt_pendVel, plt_control, layout=grid(2, 2), color=:blue, linewidth=2, size=(700,700))
