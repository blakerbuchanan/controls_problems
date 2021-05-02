module cartPoleDynamics

cd("/Users/blake/Dropbox/CMU/julia/controls/cartpoleStabilization")

using LinearAlgebra, PolynomialRoots, DifferentialEquations
# using Debugger

export cartPoleDyn, cartPoleLQRParams

function cartPoleDyn(dq,q,p,t)
    S,R,B = p
    u = -inv(R)*transpose(B)*S(t)*(q - [0;pi;0;0]);
    u = 0

    dq[1] = q[3];
    dq[2] = q[4];
    dq[3] = (1/(1+sin(q[2])^2))*(u + sin(q[2])*(q[4]^2 + cos(q[2])));
    dq[4] = (1/(1+sin(q[2])^2))*(-u*cos(q[2]) - q[4]^2 * cos(q[2])*sin(q[2]) - 2*sin(q[2]));
end

mutable struct cartPoleLQRParams
    A
    B
    Q
    R
    Tspan
end

end
