using DelimitedFiles

pathw = "../"
include(pathw*"libs/cheb.jl")
include(pathw*"libs/dist.jl")
include(pathw*"libs/dos.jl")
include(pathw*"libs/hcg.jl")
include(pathw*"libs/potential.jl")
include(pathw*"libs/shape_lib.jl")
include(pathw*"libs/slater_koster.jl")
include(pathw*"libs/sumcheb.jl")


# Script with full control over every parameter of the DOS calculation
function lowlevel_dos()

   # Properties of the nanoparticle (see documentation for options)
    mater = "gold"
    shape = "cube"
    rad = 5.1  # NP radius [nm]
    outname = "dos.dat"

    N  = 300  # total number of Chebyshev polynomials
    Nk = 2    # number of random vectors

    perc = 100    # percentage of Chebyshev moments to keep
    minE_Ha = -1  # smallest energy point
    maxE_Ha = +1  # largest energy point
    NE = 1000     # number of energy points
    flag = 2      # f=2 -> random vectors


   # Construct the TB model corresponding to the selected material
    onsite, first_neighbor, second_neighbor, A, B, fermi_Ha, a0, diel = tightbinding(mater)

   # Generate list of atomic positions
    Elist, Edict, R = generate_shape_FCC(shape, a0, rad)
    @printf("\nNumber of atoms: %d\n", length(R))

   # Use list of atomic positions to determine the Hamiltonian H [KPM units]
    H, v = slater_koster_FCC(Elist, Edict, onsite, first_neighbor, second_neighbor, A, B)

   # Build the Chebyshev vector
    mu = doscompute_mu!(H, N, Nk, flag)

   # Resum the Chebyshev vector into the DOS
    dos = get_dos(mu, perc, minE_Ha, maxE_Ha, NE, A, B)


   # Save to file
    open(outname, "w") do io
        for i in 1:size(dos,1)
            for j in 1:size(dos,2)
                @printf(io,"%12.8f ",dos[i,j])
            end
            Base.write(io,"\n")
        end
    end
end

lowlevel_dos()
