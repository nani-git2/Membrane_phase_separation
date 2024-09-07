#=
Lipids in membrane phase separating on the introduction of bacterial contacts.

Model:
2D Ising model NxN, with patches (of size pxp) introduced in one half.
Inside the patches we have a non-zero h.
=#

using StatsBase, Plots, JLD
"""
Returns list of neighbours of a point
"""
function neighbours(ind::Vector{Int64})::Vector{Vector{Int64}}
	neigh = Vector{Vector{Int64}}(undef, 4)
	i,j = ind

	# PBC for x direction exists
	neigh[1] = [mod(i,N)+1, j]
	neigh[2] = [mod(i-2,N)+1, j]
	# PBC for y direction exists
	neigh[3] = [i, mod(j,N)+1]
	neigh[4] = [i, mod(j-2,N)+1]	
	return neigh
end

"""
Calculates the effective local field Hi, at site i
Hi = -J*sum(nj)-hi
"""
function field(k::Vector{Int64}, lattice::Matrix{Int8}, J::Float64, patch::Vector{Vector{Int64}})::Float64
	Hi = 0.0
	neigh = neighbours(k)
	#contributions from the J and Js terms
	for n in neigh
		if (([k[1],k[2]] in patch) && ([n[1],n[2]] in patch)) 
			Hi += -J*Js*lattice[n[1], n[2]]
		else
			Hi += -J*lattice[n[1], n[2]]
		end
	end
	#contribution from the hi term
	if ([k[1],k[2]] in patch)
        Hi += -h
    end
	return Hi
end

"""
Calculates the total energy of the system
bmark:
"""
function calculate_energy(lattice::Matrix{Int8}, J::Float64, patch::Vector{Vector{Int64}})
    Emat = Matrix{Float64}(undef, N, N)
    for n in 1:N
        for m in 1:N
            Emat[m,n] = (field([m,n], lattice, J, patch))*(lattice[m,n])
        end
    end
    
    E = 0.5*sum(Emat)                   #since each pair counted twice
    return E
end

"""
Does a particle "displacement". Global kawasaki exchange
"""
function displace(lattice::Matrix{Int8}, J::Float64, patch::Vector{Vector{Int64}})
	k = rand(1:N, 2)         # choose random sites for global exchange            
	q = rand(1:N, 2)
	dE = 0.0
	# making the step
	if !(lattice[k[1], k[2]] == lattice[q[1], q[2]])
		Hk = field([k[1], k[2]], lattice, J, patch)
		Hq = field([q[1], q[2]], lattice, J, patch)
		dE = -2*(Hk-Hq)*(lattice[k[1], k[2]])
		dE = round(dE, digits=6)					# change in energy
		if rand() < minimum([1, exp(-dE)])	# acceptance prob.
			lattice[k[1], k[2]] = Int8(-1*lattice[k[1], k[2]])
			lattice[q[1], q[2]] = Int8(-1*lattice[q[1], q[2]])
			return lattice, dE
		else								# if Monte Carlo step is rejected
			return lattice, 0.0
	    end
	else									# if both sites have the same value
		return lattice, 0.0
	end
end



""" 
Main loop 
"""
function main(J::Float64, patch::Vector{Vector{Int64}})
    lattice = sample([Int8(-1), Int8(1)], Weights([1-phi, phi]), (N,N)) 
	energy = zeros(101)
	dE = 0.0	
	dat = Array{Int8}(undef, (N,N,101))
	dat[:,:,1] = lattice
	nsave = div(Nsteps, 100)
	energy[1] = calculate_energy(lattice, J, patch)
	for n in 2:Nsteps
		lattice, dE = displace(lattice, J, patch)
		dE = 0.0
		if n%nsave == 0 
			dat[:,:,div(n,nsave)+1] = lattice
			energy[div(n,nsave)+1] = calculate_energy(lattice, J, patch)
		end
	end
	return energy, dat
end



# Parameters
const N = 100                    # lattice size
const J = 0.3                   # 2D Ising Jc~0.441
const Nsteps = 5*10^8
const phi = 0.5                 # proportion of +1 to total

# for patches
const ps = 2                    # patch size
const h = 5.0                   # field within patch
const Np = 0                  # number of patches

#creating the patches
tempp = [Int64.([rand(1:N-ps), rand(1:(N/2-ps))]) for n in 1:Np]
patches = []
for p in tempp
    for m in p[1]:1:p[1]+ps-1
        for n in p[2]:1:p[2]+ps-1
            global patches = append!(patches, [[m,n]])
        end
    end
end
const patch = patches
