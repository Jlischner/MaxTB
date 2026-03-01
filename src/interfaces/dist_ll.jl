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


# Script with full control over every parameter of the distribution calculation
function lowlevel_dist()

   # Properties of the nanoparticle (see documentation for options)
    shape = "cube"
    mater = "gold"
    rad   = 1.5  # NP radius [nm]

    freq  = 2.4  # optical potential frequency [eV]
    eps_m = 1.0  # relative dielectric constant of the environment

    tempr = 298  # temperature [K]

    outname = "dist.dat"

   # Relaxation [eV]
    function Gamma(E,Ep)
        return 0.1
    end

    N = 200  # total number of Chebyshev polynomials

   # Number of Chebyshev polynomials per block
    NTx = 10
    NTy = 10

   # Number of blocks
    NNL = N÷NTx + 1
    NNR = N÷NTy + 1
    Nk  = 5  # number of random vectors

   # Number of integration points in the energy discretization
    Nx = N*2
    Ny = Nx

   # Percentage of Chebyshev polynomials to keep
    perc1 = 100
    perc2 = perc1

    NE  = N*2+1  # number of energy points in the integration
    NEp = N*2+1  # number of output energy points

   # Temperature
    kbT  = tempr*boltzmann  # [eV]
    beta = 1.0/kbT

    write  = false  # write debug information
    Gammaf = 1  # prefactor for the relaxation function


   # Construct the TB model corresponding to the selected material
    onsite, first_neighbor, second_neighbor, A, B, fermi_Ha, a0, diel = tightbinding(mater)

   # Get the dielectric constant from the material and the optical frequency
    eps = diel(freq)

   # Generate list of atomic positions
    Elist, Edict, R = generate_shape_FCC(shape, a0, rad)

   # Use list of atomic positions to determine the Hamiltonian H [KPM units]
    H, v = slater_koster_FCC(Elist, Edict, onsite, first_neighbor, second_neighbor, A, B)

   # Get the potential Phi [eV]
    if shape == "sphere"
        Phi = potential_sphere(R, eps, eps_m)
    elseif shape == "import"
        filename = "pot.dat"
        Phi = comsol_read(filename)
    else
        Phi = potential_sphere(R, eps, eps_m)  # calls sphere by default if no COMSOL file
    end

   # Build the Chebyshev matrix
    mumn = compute_mumn!(H, Phi, NNL, NNR, NTx, NTy, Nk)

   # Transform the Chebyshev matrix into the energy-resolved optical matrix
   # Ei,Ej [eV], opt [dimensionless]
    Ei,Ej,opt = resum_mu_dd(mumn, Nx, Ny, A, B, perc1, perc2)

   # Get the distribution
    dist = sum_dist(opt, Ei, Ej, NE, NEp, fermi_Ha, freq, beta, Gammaf, Gamma, write)


   # Write the data
    col1 = collect(dist[1])  # energy list
    col2 = collect(dist[2])  # population
    dist2 = hcat(col1,col2)
    open(outname, "w") do io
        for i in 1:size(dist2,1)
            for j in 1:size(dist2,2)
                @printf(io,"%12.8f ",dist2[i,j])
            end
        Base.write(io,"\n")
        end
    end
end

lowlevel_dist()
