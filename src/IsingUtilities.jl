module IsingUtilities

using CircularArrays
using SciPy:signal as signal

export IsingLattice # Main type used
export energy
export ΔE_at_site
export fill_magnetization_history!, fill_internalenergy_history!
export generate_lattice
export metropolis!
export estimate_eq_time
export generate_autocorr_data

include("IsingTypes.jl")

function fill_magnetization_history!(lattice::AbstractIsingLattice)
    l = copy(lattice.initial_state)
    lattice.magnetization_history = Float64[sum(l)]
    for s_k in lattice.spin_flip_history
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

function fill_internalenergy_history!(lattice::AbstractIsingLattice)
    l = copy(lattice.initial_state)
    lattice.internal_energy_history = Float64[energy(l)]
    for s_k in lattice.spin_flip_history
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

function estimate_eq_time(lattice::IsingLattice; deviation_tolerance = 1e-5, estimator = :magnetization_history,
    step_sweep_size = 10, iterlim = Int(1e5))

    if isempty(getfield(lattice, estimator))
        
        if estimator == :magnetization_history
            fill_magnetization_history!(lattice)
        elseif estimator == :internal_energy_history
            fill_internalenergy_history!(lattice) end 
        # If estimator is still missing, exit
        if isempty(getfield(lattice, estimator)[2:end])
            @error "Unable to fill estimator data! Exiting"
            throw(EOFError())
        end
    end
    
    
    m = getfield(lattice, estimator) ./ prod(size(lattice))
    mc_perc  = 0:1/length(m):1
    
    
    N, N = size(lattice.initial_state)
    step_size = floor(Int, step_sweep_size * N^2)
    
    # slope = (y₂ - y₁)/(x₂ - x₁)
    slope(p1, p2) = (p2[2] - p1[2]) / (p2[1] - p1[1])
    slope_min, slope_max = -deviation_tolerance, deviation_tolerance
    
    # Pick the starting point of my slope calculation
    τ_eq = 1
    p2 = (mc_perc[end], m[end])
    p1 = (mc_perc[τ_eq], m[τ_eq])
    
    niters = 0
    while !(slope_min <= slope(p1, p2) <= slope_max)
        
        
        τ_eq += step_size
        if τ_eq > length(m)
            @error "Unable to converge!"
            return -1
        end
        
        if niters > iterlim
            @warn "Unable to converge to thermalization time within given iteration limit! Ensure that the system has thermalized or increase iteration limit."
            return -1
        end
        
        
        p1 = (mc_perc[τ_eq], m[τ_eq])
        
        niters += 1
    end
    
    return τ_eq ./ N^2
    
end

# Naieve autocorrelation function
function slow_autocorr(magnetization_history, t)
    # t_max is the number of monte carlo steps we've taken in total
    t_max = length(magnetization_history)
    
    # Magnetization history m, assuming I've run fill_magnetization_history! already on the lattice
    m = magnetization_history
    
    # Normalization constant
    norm = 1/(t_max - t)
    
    mean_m² = norm * sum(m[1:t_max-t]) * norm * sum(m[t+1:t_max])
    
    return norm * sum(m[1:t_max - t] .* m[t+1:t_max]) - mean_m²
end

function generate_autocorr_data(lattice::IsingLattice, eq_time; n_sweeps = 100, show_progress = true)
    # Check if magnetization data is available
    if isempty(lattice.magnetization_history)
        @info "Filling magnetization data"
        fill_magnetization_history!(lattice)
    end
    
    N = size(lattice.initial_state)[1]
    m = lattice.magnetization_history[ceil(Int, eq_time*N^2):end]
    npts = floor(Int, n_sweeps*N^2)
    
    # Create an array of corresponding χ values
    χ = zeros(npts)
    
    if show_progress
        P = Progress(npts)
        Threads.@threads for t in 1:npts
            χ[t] = slow_χ(m, t)
            next!(P)
        end
    else
        Threads.@threads for t in 1:npts
            χ[t] = slow_χ(m, t)
        end
    end
    
    return χ
end

include("IsingAlgorithms.jl")

end # module
