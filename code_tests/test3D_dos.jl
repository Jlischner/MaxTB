using DelimitedFiles

pathw = "../src/"
include(pathw*"libs/cheb.jl")
include(pathw*"libs/dist.jl")
include(pathw*"libs/dos.jl")
include(pathw*"libs/hcg.jl")
include(pathw*"libs/potential.jl")
include(pathw*"libs/shape_lib.jl")
include(pathw*"libs/slater_koster.jl")
include(pathw*"libs/sumcheb.jl")


function test(L,N, mater)

   # L: number of crystal planes for the PBC
   # N: total number of Chebyshev polynomials
   # mater: material

   # Properties of the nanoparticle (see documentation for options)
    shape = "periodic"

   # Construct the TB model corresponding to the selected material
    onsite, first_neighbor, second_neighbor, A, B, fermi_Ha, a0, diel = tightbinding(mater)

    rad = (L-1+0.01)*a0/2  # NP radius [nm]

   # Generate list of atomic positions
    Elist, Edict, R = generate_shape_FCC(shape, a0, rad)
    @printf("\nNumber of atoms: %d\n", length(R))

   # Use list of atomic positions to determine the Hamiltonian H [KPM units]
    H, v = slater_koster_FCC(Elist, Edict, onsite, first_neighbor, second_neighbor, A, B, L1=L, L2=L, L3=L, periodic=true)

   # Build the Chebyshev matrix
    flag = 2    # use random vectors for the calculation
    Nk   = 10   # number of random vectors
    mu   = doscompute_mu!(H, N, Nk, flag)

   # Get the DOS from the Chebyshev matrix
    perc = 100  # percentage of polynomials to keep
    NE = 1000  # number of energy points

   # Set large limits to trigger the automatic DOS rescale inside get_dos
    minE_Ha = -10  # smallest energy point
    maxE_Ha = +10  # largest energy point

   # Resum the Chebyshev vector into the DOS
    dos = get_dos(mu, perc, minE_Ha, maxE_Ha, NE, A, B)

    return dos
end


println(length(ARGS))
if length(ARGS) > 0 && ARGS[1] == "run"

   # Properties of the nanoparticle (see documentation for options)
    L = 8
   # L = 128
    N = 200
    mater = "3DTB"
    dos = test(L, N, mater)

   # Save to file
    open("dos.txt", "w") do io
        for i in 1:size(dos,1)
            for j in 1:size(dos,2)
                @printf(io,"%12.8f ",dos[i,j])
            end
            Base.write(io,"\n")
        end
    end

#    open("dos.txt", "w") do io
#        writedlm(io, dos)
#    end

end
