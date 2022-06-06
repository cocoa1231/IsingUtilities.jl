# IsingUtilities.jl
A package containing algorithms and utilites I'm using for my undergraduate thesis regarding the Ising Model

**NOTE** -- This is NOT meant to be a package available for public use. This is only to make my life easier.
The code in this package is not well written and may not be easily accessibly and is definitely not well tested.
You are better off using this as a reference to write your own code at best. If you use this, expect things to
break.

---

## Basic Usage

This package exports a main type `IsingLattice` which is a subtype of `AbstractIsingLattice` on which, at the moment
these functions can act -

```julia
julia> methodswith(IsingUtilities.AbstractIsingLattice)
[1] energy(lattice::IsingUtilities.AbstractIsingLattice) in IsingUtilities at /home/cocoafedora/.julia/packages/IsingUtilities/njdYe/src/IsingAlgorithms.jl:3
[2] fill_internalenergy_history!(lattice::IsingUtilities.AbstractIsingLattice) in IsingUtilities at /home/cocoafedora/.julia/packages/IsingUtilities/njdYe/src/IsingUtilities.jl:32
[3] fill_magnetization_history!(lattice::IsingUtilities.AbstractIsingLattice) in IsingUtilities at /home/cocoafedora/.julia/packages/IsingUtilities/njdYe/src/IsingUtilities.jl:15
[4] metropolis!(lattice::IsingUtilities.AbstractIsingLattice, steps::Integer, β::Float64) in IsingUtilities at /home/cocoafedora/.julia/packages/IsingUtilities/njdYe/src/IsingAlgorithms.jl:21
```

The first three are utility functions, providing the ability to fill the magnetization history based on single-spin-flip history.
In the future I will be adding cluster algorithms too. But at the moment one can evolve the lattice a given number of Monte Carlo
steps using the `metropolis!` mutating function which will fill in the spin flip history. Then one can fill in the magnetization
and internal energy history of the lattice and take measurements on that. Basic usage goes something like this -

```julia
julia> using IsingUtilities

julia> l = IsingLattice(generate_lattice(50, :infty))
IsingLattice([1 -1 … 1 1; 1 1 … 1 -1; … ; 1 -1 … -1 -1; -1 1 … 1 -1], [1 -1 … 1 1; 1 1 … 1 -1; … ; 1 -1 … -1 -1; -1 1 … 1 -1], Tuple{Int64, Int64}[], Float64[], Float64[])

julia> metropolis!(l, Int(1e5), 0.9)
Progress: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████| Time: 0:00:00
IsingLattice([1 -1 … 1 1; 1 1 … 1 -1; … ; 1 -1 … -1 -1; -1 1 … 1 -1], [1 1 … 1 1; 1 1 … 1 1; … ; 1 1 … 1 1; 1 1 … 1 1], [(-1, -1), (-1, -1), (43, 29), (-1, -1), (39, 1), (-1, -1), (-1, -1), (-1, -1), (29, 39), (49, 39)  …  (-1, -1), (-1, -1), (-1, -1), (-1, -1), (-1, -1), (-1, -1), (-1, -1), (-1, -1), (-1, -1), (-1, -1)], Float64[], Float64[])

julia> fill_internalenergy_history!(l);

julia> fill_magnetization_history!(l);

julia> # Take your measurements on the lattice now!
```