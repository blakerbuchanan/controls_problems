## Implementations of control strategies in Julia

This respository is home to implementations of various nonlinear control strategies on simple underactuated mechanical systems. *Currently only the cartpole_stabilization directory is in a working state.* Note that this repository undergoes continous updates. If you are an employer interested in assessing  these implementations for the purposes of reviewing me as a potential candidate, firstly, thanks for taking a look. Secondly, if you are unfamiliar with the Julia programming language, please read on. Julia is incredibly readable as a programming language, and is much like Python in terms of its readability. My choice to implement these problems in Julia rather than Python is due to the fast rate at which Julia is growing, and I imagine its usage will be on par with Python in the coming years. It is also reported to be considerably faster than Python on several benchmarks.

## Downloading and executing

First download the latest version of the Julia programming language at https://julialang.org/. 

Next, download the controls_problems directory.

Navigate to the downloaded directory and start a new Julia REPL by running ```julia``` in the terminal. Execute the following in the REPL.

``` using Pkg``` 

To install all of the required packages, execute ```Pkg.add("DifferentialEquations")```, ```Pkg.add("LinearAlgebra")```, ```Pkg.add("Plots")```, ```Pkg.add("MatrixEquations")```, ```Pkg.add("LaTeXStrings")```, and ```Pkg.add("Dierckx")```.

You should now be able to execute julia scripts like ```cartpole_LQR_stabilization.jl``` from the Julia REPL. 

