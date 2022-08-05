using ProgressMeter

function energy(lattice::CircularArray)
    kern = [
        0 1 0;
        1 0 1;
        0 1 0
    ]
    atom_energy = signal.convolve2d(lattice, kern, mode = "same", boundary = "wrap") .* lattice
    return -Int(sum(atom_energy)/2)
end

function ΔE_at_site(lattice::CircularArray, k)
    x, y = k
    sk = lattice[x, y]
    nn_sum = sum([lattice[x+1,y], lattice[x-1, y], lattice[x, y+1], lattice[x, y-1]])
    return 2*sk*nn_sum
end

# Metropolis Algorithm for evolving the lattice
function metropolis!(lattice::AbstractIsingLattice, steps::Integer, β::Float64; showprogress = false)
    
    possible_EChange = [-8:2:8;]
    exponentials = Dict(possible_EChange .* -β .=> exp.(possible_EChange .* -β))
    N = size(lattice.initial_state)[1]
    
    p = Progress(steps)
    for step in 1:steps
        
        # Pick a random lattice point
        x, y = rand(1:N, 2)
        
        # Calculate dE
        dE = ΔE_at_site(lattice.final_state, (x, y))
        
        # If new energy is lower accept
        if dE < 0
            lattice.final_state[x, y] *= -1
            push!(lattice.spin_flip_history, (x, y))
        # Otherwise probabilistically accept the new state with probability exp(-β*dE)
        elseif dE >= 0
            u = rand()
            try
                if u < exponentials[-β*dE]
                    lattice.final_state[x, y] *= -1
                    push!(lattice.spin_flip_history, (x, y))
                else
                    push!(lattice.spin_flip_history, (-1, -1))
                end
            catch(e)
                @warn e
                if isa(e, KeyError)
                    exponentials[-β*dE] = exp(-β*dE)
                    if u < exponentials[-β*dE]
                        lattice.final_state[x, y] *= -1
                        push!(lattice.spin_flip_history, (x, y))
                    end
                end
            end
        end
        if showprogress == true
            next!(p)
        end
    end
    
    return lattice

end
