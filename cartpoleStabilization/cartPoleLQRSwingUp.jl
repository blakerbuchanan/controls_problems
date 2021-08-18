# Energy-based swing up with LQR stamilization implementation of the cart-pole system
using DifferentialEquations, LinearAlgebra, Dierckx
include("cartPoleDynamics.jl")
using .cartPoleDynamics
using Plots
using MatrixEquations
using Debugger, LaTeXStrings
logocolors = Colors.JULIA_LOGO_COLORS

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
T = 30.0;
Tspan = (0.0, T)
deltaT = 0.001;
tspan = range(0.0,stop=T,step=deltaT)

priccati = cartPoleLQRParams(A,B,Q,R,tspan)
alg = ExplicitRK(tableau=constructDormandPrince())
riccatiProb = ODEProblem(riccati, S0, Tspan, priccati);
riccatiSol = solve(riccatiProb, alg, dt = deltaT, saveat=deltaT);

# Numerically estimate basin of attraction for cartpole system
thmin = -pi/6;
thmax = pi/6;
deltath = 0.1;

thbase = thmin:deltath:thmax;
thdbase = thmin:deltath:thmax;
basin = [0,pi,0,0];

thd = 0;
th = 0.01;
# E0 = (1/2)*thd^2 - cos(th);

pcart = riccatiSol, R, B

# for i in thbase
#     for j in thdbase
#
#         x0 = [0, pi + i, 0, j];
#         prob = ODEProblem(cartPoleDynLQR, x0, Tspan, pcart);
#         sol = solve(prob, alg, dt = deltaT, saveat=0.01);
#
#         if norm(sol(T) - [0, pi, 0, 0]) < 0.01
#             basin = [basin x0];
#         end
#     end
# end

kE = 8;
kp = 0.5;
kd = 0.5;

Kswing = [kE, kp, kd];

pcart = riccatiSol, R, B, Kswing, basin
x0 = [0, -0.01, 0, 0];

alg = Tsit5()
prob = ODEProblem(cartPoleDynSwingUp,x0,Tspan,pcart);
solSwing = solve(prob, alg, dt = deltaT, saveat=0.01);
plot(solSwing, xaxis = (L"t"), linewidth=1.5, label = [L"x" L"$\theta$" L"$\dot{x}$" L"$\dot{\theta}$"],legendfontsize=12, legend=:bottomright)

thetaData = [solSwing(t)[2] for t in range(0,stop=T,length=2000)]
thetaDotData = [solSwing(t)[4] for t in range(0,stop=T,length=2000)]

plot(thetaData,thetaDotData,linewidth=1.5, color=logocolors.blue,xguidefontsize=18,yguidefontsize=18,xaxis = (L"\theta"),yaxis=(L"\dot{\theta}"),legend=false,grid=true)
plot!(circleShape(-0.01,0.0,0.05),color=logocolors.green,linewidth=1.5,grid=true)
plot!(circleShape(pi,0.0,0.05),color=logocolors.red,linewidth=1.5,aspect_ratio=1.0,grid=true)
savefig("phaseportrait.png")
