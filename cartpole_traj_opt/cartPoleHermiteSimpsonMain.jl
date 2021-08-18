# Cart-pole control using Hermite-Simpson collocation
using NLopt, LinearAlgebra
include("cartPoleDynamics.jl")
using .cartPoleDynamics

# model = Model()

# Initial conditions
N_st = 4;
N = 50;
dt = 0.1;
t = 0:dt:N*dt;
# tc = range(t[1]+dt/2, stop=t[end]-dt/2, length=N);

x = zeros(N*N_st);
u = zeros(N);
Xopt = [u; x];
X0 = rand(N*N_st + N)
# testcon = constraintsFunc(X0,N,N_st,dt)
# testobj = costFunction(X0)

lbc = -5*ones(N); # lower bound on control
ubc = 5*ones(N); # upper bound on control
lbs = -Inf*ones(N*N_st); # lower bound on states
ubs = Inf*ones(N*N_st); # upper bound on states
lb = [lbc; lbs]; # lower bounds
ub = [ubc; ubs]; # lower bounds

# LD_SLSQP
# LN_COBYLA
opt = Opt(:LD_SLSQP, N_st*N+N)
opt.lower_bounds = lb;
opt.upper_bounds = ub;
# opt.xtol_rel = 0.00000001
opt.ftol_rel = 0.00000001
opt.maxeval = 100000
# opt.initial_step = 0.001*ones(N*N_st + N)

opt.min_objective = costFunction
# inequality_constraint!(opt, (x) -> constraintsFunc(x,u), 1e-8)
equality_constraint!(opt, (X) -> constraintsFunc)
X0 = zeros(N*N_st + N)
(minf,minx,ret) = optimize(opt, X0)
numevals = opt.numevals # the number of function evaluations
println("got $minf at $minx after $numevals iterations (returned $ret)")

function constraintsFunc(X)
    N_st = 4;
    N = 50;
    dt = 0.1;
    x = X[N+1:end];
    u = X[1:N];
    x = reshape(x,(N,N_st));

    xdot = cartPoleDynv2(x,u);
    xll = x[1:end-1,:];
    xrr = x[2:end,:];
    xdotll = xdot[1:end-1,:];
    xdotrr = xdot[2:end,:];
    ull = u[1:end-1,:];
    urr = u[2:end,:];
    xc = .5*(xll .+ xrr) .+ dt/8*(xdotll .- xdotrr);
    uc = (ull .+ urr)./2;
    xdotc = cartPoleDynv2(xc,uc);

    ceq = (xll .- xrr) .+ dt/6*(xdotll .+ 4*xdotc .+ xdotrr);
    # print(x[end,:])
    xT = x[end,:] - [0,pi,0,0];
    # xT = xT';
    x1 = x[1,:];
    println(size(ceq))
    println(size(xT))
    println(size(x1))

    # ceq = reshape(ceq,(N*N_st))
    ceq = vec(ceq)
    ceq = [ceq; xT; x1];

    # print(size(ceq))
    # c = [];
    return ceq
end

# function cartPoleDynExpr(model,q,u)
#     # How do I incoporate the control "u"?
#     # dq = zeros(size(q));
#     dq1 = @NLexpression(model,q[3]);
#     dq2 = @NLexpression(model,q[4]);
#     dq3 = @NLexpression(model, (1/(1+sin(q[2])^2))*(u + sin(q[2])*(q[4]^2 + cos(q[2]))) );
#     dq4 = @NLexpression(model, (1/(1+sin(q[2])^2))*(-u*cos(q[2]) - q[4]^2 * cos(q[2])*sin(q[2]) - 2*sin(q[2])));
#     return dq1,dq2,dq3,dq4
# end

function costFunction(u)
    N = 50;
    J = dot(u,u);
    return J
end

function f = inBasin(q,basin)
    conv = q - basin;
    norms = sqrt(sum(conv.^2,2));
    min(norms)
    if min(norms) < 0.1
        f = 1;
    else
        f = 0;
    end
end
