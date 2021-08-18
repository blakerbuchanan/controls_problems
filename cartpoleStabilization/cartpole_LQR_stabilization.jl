# Use LQR to stabilize the cart-pole system within the basin of attraction of the unstable fixed point
using DifferentialEquations, LinearAlgebra, Dierckx, Plots, MatrixEquations, LaTeXStrings

include("cartPoleDynamics.jl")
using .cartPoleDynamics

logocolors = Colors.JULIA_LOGO_COLORS # Use the Julia Logo colors for plotting

function riccati(dS,S,p,t)
    A = p.A; B = p.B; Q = p.Q; R = p.R; tspan = p.Tspan;
    dS .= (Q .- S*B*inv(R)*transpose(B)*S .+ S*A .+ transpose(A)*S);
end

# Linearized cart-pole dynamics with parameters equal to 1
A = [0 0 1 0; 0 0 0 1; 0 1 0 0; 0 2 0 0]
B = [0; 0; 1; 1]

# Define cost matrix Q
Q = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
Q[1,1] = 1
Q[3,3] = 1
R = 1

# Solve the Riccati equation
S0 = 10.0 .*Q
T = 30.0
Tspan = (0.0, T)
deltaT = 0.001;
tspan = range(0.0,stop=T,step=deltaT)

priccati = cartPoleLQRParams(A,B,Q,R,tspan)
alg = ExplicitRK(tableau=constructDormandPrince())
riccatiProb = ODEProblem(riccati, S0, Tspan, priccati);
riccatiSol = solve(riccatiProb, alg, dt = deltaT, saveat=deltaT);

# Define parameters and initial conditions
pcart = riccatiSol, R, B
x0 = [0, pi+pi/12, 0, 0];

# Generate an ODEProblem and solve
alg = Tsit5()
prob = ODEProblem(cartPoleDynLQR, x0, Tspan, pcart);
sol = solve(prob, alg, dt = deltaT, saveat=0.01);

# Visualize results
plot(sol, xaxis = (L"t"), linewidth=1.5, label = [L"x" L"$\theta$" L"$\dot{x}$" L"$\dot{\theta}$"],legendfontsize=12, legend=:bottomright)

thetaData = [sol(t)[2] for t in range(0,stop=T,length=2000)]
thetaDotData = [sol(t)[4] for t in range(0,stop=T,length=2000)]

plot(thetaData,thetaDotData,linewidth=1.5, color=logocolors.blue,xguidefontsize=18,yguidefontsize=18,xaxis = (L"\theta"),yaxis=(L"\dot{\theta}"),legend=false,grid=true)
plot!(circleShape(-0.01,0.0,0.05),color=logocolors.green,linewidth=1.5,grid=true)
plot!(circleShape(pi,0.0,0.05),color=logocolors.red,linewidth=1.5,aspect_ratio=1.0,grid=true)
savefig("phaseportrait.png")
