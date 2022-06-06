module IsingUtilities

using CircularArrays
using SciPy:signal as signal

export IsingLattice # Main type used
export energy
export ΔE_at_site
export fill_magnetization_history!, fill_internalenergy_history!
export generate_lattice
export metropolis!

include("IsingTypes.jl")

function fill_magnetization_history!(lattice::IsingLattice)
    l = copy(lattice.initial_state)
    lattice.magnetization_history = Float64[sum(l)]
    for (idx, s_k) in enumerate(lattice.spin_flip_history)
        if s_k == (-1, -1)
            push!(lattice.magnetization_history, lattice.magnetization_history[end])
        else
            # Update the lattice for the next calculation
            l[s_k...] *= -1

            # Calculate Mν and append it
            push!(lattice.magnetization_history, lattice.magnetization_history[end] + 2*l[s_k...])
        end
    end
    return lattice.magnetization_history
end

function fill_internalenergy_history!(lattice::IsingLattice)
    l = copy(lattice.initial_state)
    lattice.internal_energy_history = Float64[energy(l)]
    for (idx, s_k) in enumerate(lattice.spin_flip_history)
        if s_k == (-1, -1)
            push!(lattice.internal_energy_history, lattice.internal_energy_history[end])
        else
            # Calculate ΔE at the given site and append it
            push!(lattice.internal_energy_history, lattice.internal_energy_history[end] + ΔE_at_site(l, s_k))

            # Update the lattice for the next calculation
            l[s_k...] *= -1
        end
    end
    return lattice.internal_energy_history 
end

function generate_lattice(N::Integer, T::Symbol)
    if T == :zero
        lattice = -1 .* ones(N, N) |> CircularArray .|> Integer
    elseif T == :infty
        lattice = rand([-1, 1], N, N) |> CircularArray .|> Integer
    else
        @error "Valid initial temperature is :zero and :infty"
    end

    return lattice
end

include("IsingAlgorithms.jl")

end # module
