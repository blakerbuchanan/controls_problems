# LQR implementation of the cart-pole system
using DifferentialEquations, LinearAlgebra, Dierckx
include("cartPoleDynamics.jl")
using .cartPoleDynamics
using Plots
using MatrixEquations
using Debugger, LaTeXStrings

function riccati(dS,S,p,t)
    A = p.A; B = p.B; Q = p.Q; R = p.R; tspan = p.Tspan;
    # At = Spline1D([0.0, 15.0], A'; k=1)(t)
    # Bt = Spline1D(tspan, B'; k=1)(t)

    dS .= (Q .- S*B*inv(R)*transpose(B)*S .+ S*A .+ transpose(A)*S);
end

# Linearized cart-pole dynamics with parameters equal to 1
A = [0 0 1 0; 0 0 0 1; 0 1 0 0; 0 2 0 0];
B = [0; 0; 1; 1];

Q = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
Q[1,1] = 1;
Q[3,3] = 1;
R = 1;

S0 = 10.0 .*Q;
T = 10.0;
Tspan = (0.0, T)
deltaT = 0.001;
tspan = range(0.0,stop=T,step=deltaT)

priccati = cartPoleLQRParams(A,B,Q,R,tspan)
alg = ExplicitRK(tableau=constructDormandPrince())
riccatiProb = ODEProblem(riccati, S0, Tspan, priccati);
riccatiSol = solve(riccatiProb, alg, dt = deltaT, saveat=deltaT);
pcart = riccatiSol, R, B

x0 = [0; pi+pi/12 ;0;0];
prob = ODEProblem(cartPoleDyn,x0,Tspan,pcart);
sol = solve(prob,alg,dt = deltaT, saveat=0.01);
plot(sol, xaxis = (L"t"), linewidth=1.5, label = [L"x" L"$\theta$" L"$\dot{x}$" L"$\dot{\theta}$"],legendfontsize=12, legend=:right)
