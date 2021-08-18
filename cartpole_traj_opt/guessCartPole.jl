# Get an initial guess
using DifferentialEquations, LinearAlgebra, Ipopt, Interpolations, Plots, JLD2
include("cartPoleDynamics.jl")
using .cartPoleDynamics

cd("/Users/blake/Dropbox/CMU/julia/controlsProblems/cartPoleTrajOpt")

# Get an initial guess
alg = ExplicitRK(tableau=constructDormandPrince())
deltaT = 5.0/201
Tspan = (0.0, 5.0-deltaT)

x0 = [0, pi/2, 0, 0]
S = 0; R = 0; B = 0;
p = S,R,B;
guessProb = ODEProblem(cartPoleDyn, x0, Tspan, p);
guessSol = solve(guessProb, alg, dt = deltaT, saveat=deltaT);
controlSol = [0*rand() for i in 1:length(guessSol.t)]

solOutput = convert(Array,guessSol)
solOutput = [solOutput' controlSol]
@save "data/guessSolution.jld2" solOutput