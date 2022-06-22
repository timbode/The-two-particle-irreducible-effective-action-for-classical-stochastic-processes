# The Two-Particle Irreducible Effective Action for Classical Stochastic Processes

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to enable reproduction of the numerical results of the research article 

[The Two-Particle Irreducible Effective Action for Classical Stochastic Processes](https://doi.org/10.1088/1751-8121/ac73c6).

It is authored by Tim Bode.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
