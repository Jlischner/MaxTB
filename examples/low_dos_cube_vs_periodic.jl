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


function dos()

    N  = 300  # total number of Chebyshev polynomials
    Nk = 2    # number of random vectors

    perc = 100    # percentage of Chebyshev moments to keep
    minE_Ha = -1  # smallest energy point
    maxE_Ha = +1  # largest energy point
    NE = 1000     # number of energy points
    flag = 2      # f=2 -> random vectors

   # Select material
    mater = "gold"

   # Construct the TB model corresponding to the selected material
    onsite, first_neighbor, second_neighbor, A, B, fermi_Ha, a0, diel = tightbinding(mater)

    shape = "periodic"
    L = 8  # number of crystal planes for the PBC
    rad = (L-1+0.01)*a0/2  # NP radius [nm]

   # Generate list of atomic positions
    Elist, Edict, R = generate_shape_FCC(shape, a0, rad)
    @printf("\nNumber of atoms: %d\n", length(R))

   # Use list of atomic positions to determine the Hamiltonian H [KPM units]
    H, v = slater_koster_FCC(Elist, Edict, onsite, first_neighbor, second_neighbor, A, B, L1=L, L2=L, L3=L, periodic=true)

   # Build the Chebyshev vector
    mu = doscompute_mu!(H, N, Nk, flag)

   # Resum the Chebyshev vector into the DOS
    dos_periodic = get_dos(mu, perc, minE_Ha, maxE_Ha, NE, A, B)





    shape = "cube"
    rad = 4.1  # size of the cube

   # Generate list of atomic positions
    Elist, Edict, R = generate_shape_FCC(shape, a0, rad)
    @printf("\nNumber of atoms: %d\n", length(R))

   # Use list of atomic positions to determine the Hamiltonian H [KPM units]
    H, v = slater_koster_FCC(Elist, Edict, onsite, first_neighbor, second_neighbor, A, B)

   # Build the Chebyshev vector
    mu = doscompute_mu!(H, N, Nk, flag)

   # Resum the Chebyshev vector into the DOS
    dos_cube = get_dos(mu, perc, minE_Ha, maxE_Ha, NE, A, B)


   # Save to file
    open("dos_periodic.dat", "w") do io
        for i in 1:size(dos_periodic,1)
            for j in 1:size(dos_periodic,2)
                @printf(io,"%12.8f ",dos_periodic[i,j])
            end
            Base.write(io,"\n")
        end
    end

    open("dos_cube.dat", "w") do io
        for i in 1:size(dos_cube,1)
            for j in 1:size(dos_cube,2)
                @printf(io,"%12.8f ",dos_cube[i,j])
            end
            Base.write(io,"\n")
        end
    end

#    open("dos_periodic.dat", "w") do io
#        writedlm(io, dos_periodic)
#    end
#    open("dos_cube.dat", "w") do io
#        writedlm(io, dos_cube)
#    end

end

dos()
