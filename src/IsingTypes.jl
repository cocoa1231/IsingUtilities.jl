abstract type AbstractIsingLattice end

@Base.kwdef mutable struct IsingLattice <: AbstractIsingLattice
    initial_state
    final_state
    spin_flip_history
    internal_energy_history
    magnetization_history
end

IsingLattice(lattice) = IsingLattice(copy(lattice), copy(lattice), Vector{Tuple{Int64, Int64}}(), Float64[], Float64[])
Base.size(l::IsingLattice) = size(l.initial_state)