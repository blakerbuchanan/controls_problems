module cartPoleDynamics

cd("/Users/blake/Dropbox/CMU/julia/controls/cartpoleStabilization")

using LinearAlgebra, PolynomialRoots, DifferentialEquations
# using Debugger

export cartPoleDyn, cartPoleDynLQR, cartPoleDynSwingUp, cartPoleLQRParams, inBasin

function cartPoleDyn(dq,q,p,t)
    S,R,B = p
    u = 0

    dq[1] = q[3];
    dq[2] = q[4];
    dq[3] = (1/(1+sin(q[2])^2))*(u + sin(q[2])*(q[4]^2 + cos(q[2])));
    dq[4] = (1/(1+sin(q[2])^2))*(-u*cos(q[2]) - q[4]^2 * cos(q[2])*sin(q[2]) - 2*sin(q[2]));
end

function cartPoleDynLQR(dq,q,p,t)
    S,R,B = p
    u = -inv(R)*transpose(B)*S(t)*(q - [0;pi;0;0]);
    # u = 0

    # E = 0.5*q[4]^2 - cos(q[2]) - 1;
    # u = Kswing[1]*q[4]*cos(q[2])*E - Kswing[2]*q[1] - Kswing[3]*q[3];
    # flag = inBasin(q, basin);
    # if flag == true
    #     u = -inv(R)*transpose(B)*S(t)*(q - [0;pi;0;0]);
    # end

    dq[1] = q[3];
    dq[2] = q[4];
    dq[3] = (1/(1+sin(q[2])^2))*(u + sin(q[2])*(q[4]^2 + cos(q[2])));
    dq[4] = (1/(1+sin(q[2])^2))*(-u*cos(q[2]) - q[4]^2 * cos(q[2])*sin(q[2]) - 2*sin(q[2]));
end

function cartPoleDynSwingUp(dq,q,p,t)
    S,R,B,K,basin = p
    # q[2] = mod2pi(q[2])
    qcheck = [q[1],mod2pi(q[2]),q[3],q[4]]

    E = 0.5*q[4]^2 - cos(qcheck[2]) - 1;
    u = K[1]*q[4]*cos(qcheck[2])*E - K[2]*q[1] - K[3]*q[3];
    flag = inBasin(qcheck, basin);

    if flag == true
        u = -inv(R)*transpose(B)*S(t)*(qcheck - [0, pi, 0, 0]);
        print("LQR on")
    end

    dq[1] = q[3];
    dq[2] = q[4];
    dq[3] = (1/(1+sin(q[2])^2))*(u + sin(q[2])*(q[4]^2 + cos(q[2])));
    dq[4] = (1/(1+sin(q[2])^2))*(-u*cos(q[2]) - q[4]^2 * cos(q[2])*sin(q[2]) - 2*sin(q[2]));
end

function inBasin(q,basin)
    # conv = q .- basin;
    # norms = sqrt.(sum(conv.^2,dims=2));

    # if minimum(norms) < 0.01
    #     return true;
    # else
    #     return false;
    # end

    # if norm(q - [0, pi, 0, 0]) < 0.01
    #     return true
    # else
    #     return false
    # end

    if q[2] > pi-pi/6 && q[2] < pi+pi/6
        return true
    else
        return false
    end
end

mutable struct cartPoleLQRParams
    A
    B
    Q
    R
    Tspan
end

end
