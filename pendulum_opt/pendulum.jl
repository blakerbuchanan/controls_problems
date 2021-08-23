# The formulation for dynamics follows the format provided by the RobotDynamics tutorial from the Robotic Exploration lab at CMU
# Using this as a basis for learning a bit more about Julia best practices
using RobotDynamics
using StaticArrays

# Define the model struct with parameters
struct Pendulum{T} <: AbstractModel
    mp::T
    b::T
    l::T
    g::T
end

Pendulum() = Pendulum(0.2, 0.0, 0.5, 9.81)

# Define the continuous dynamics
function RobotDynamics.dynamics(model::Pendulum, x, u)
    mp = model.mp   # mass of the pole (point mass at the end) in kg
    l = model.l   # length of the pole in m
    g = model.g  # gravity m/s^2

    q = x[ @SVector [1] ] # What is an SVector?
    qd = x[ @SVector [2] ]

    s = sin(q[1])
    c = cos(q[1])

    H = @SMatrix [mp*l^2]
    C = @SMatrix [b]
    G = @SVector [mp*g*l*s]
    B = @SVector [1]

    qdd = -H\(C*qd + G - B*u[1])
    return [qd; qdd]
end

# Specify the state and control dimensions
RobotDynamics.state_dim(::Pendulum) = 2
RobotDynamics.control_dim(::Pendulum) = 1

# Create the model
model = Pendulum()
n,m = size(model)

# Generate random state and control vector
x,u = rand(model)
dt = 0.1  # time step (s)
z = KnotPoint(x,u,dt)

# Evaluate the continuous dynamics and Jacobian
ẋ = dynamics(model, x, u)
∇f = RobotDynamics.DynamicsJacobian(model)
jacobian!(∇f, model, z)

# Evaluate the discrete dynamics and Jacobian
x′ = discrete_dynamics(RK3, model, z)
discrete_jacobian!(RK3, ∇f, model, z)
