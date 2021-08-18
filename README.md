## Implementations of control strategies in Julia

This respository is home to implementations of various nonlinear control strategies on simple underactuated mechanical systems. *Currently only the cartpoleStabilization directory is in a working state.* Note that this repository undergoes continous updates.

## Downloading and executing

First download the latest version of the Julia programming language at https://julialang.org/. 

Next, download the controls_problems directory.

Navigate to the downloaded directory and start a new Julia REPL by running ```julia``` in the terminal. Execute the following in the REPL.

``` using Pkg``` 

To install all of the required packages, execute ```Pkg.add("DifferentialEquations")```, ```Pkg.add("LinearAlgebra")```, ```Pkg.add("Plots")```, ```Pkg.add("MatrixEquations")```, ```Pkg.add("LaTeXStrings")```, and ```Pkg.add("Dierckx")```.

You should now be able to execute julia scripts like ```cartpole_LQR_stabilization.jl``` from the Julia REPL. 

