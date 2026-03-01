using DelimitedFiles

pathw = "../src/"
include(pathw*"libs/cheb.jl")
include(pathw*"libs/dist.jl")
include(pathw*"libs/hcg.jl")
include(pathw*"libs/potential.jl")
include(pathw*"libs/shape_lib.jl")
include(pathw*"libs/slater_koster.jl")
include(pathw*"libs/sumcheb.jl")


# Script with full control over every parameter of the distribution calculation
function raw()

   # Properties of the nanoparticle (see documentation for options)
    shape = "cube"
    mater = "3DTB"
    rad   = 1.1  # NP radius [nm]

    freq  = 2.0  # optical potential frequency [eV]
    eps_m = 1.0  # relative dielectric constant of the environment

    tempr = 298  # temperature [K]

    outname = "dist.dat"

   # Relaxation [eV]
    function Gamma(E,Ep)
        return 0.1
    end

    N = 400  # total number of Chebyshev polynomials

   # Hot carrier generation parameters
    kappa = 0.0  # parameter in the broadening function
    gph   = 0.06  # broadening [eV]

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

   # Get the hot carrier generation rate
    hcg  = sum_hcg( opt, Ei, Ej, NE, NEp, fermi_Ha, freq, beta, kappa, gph, write)

   # columns = hcg_conv(mumn, A, B, NE, fermi, freq, beta, kappa, gph, write)
    print_conv = true
    if print_conv
        hcg_conv(mumn, A, B, NE, fermi_Ha, freq, beta, kappa, gph, print_conv)
    end

   # Write the data
    col1 = collect(hcg[1])   # energy list
    col2 = collect(hcg[2])   # hot electron
    col3 = collect(hcg[3])   # hot hole
    col4 = collect(dist[2])  # population
    columns = hcat(col1,col2,col3,col4)
    open(outname, "w") do io
        for i in 1:size(columns,1)
            for j in 1:size(columns,2)
                @printf(io,"%12.8f ",columns[i,j])
            end
        Base.write(io,"\n")
        end
    end

#    columns = Array{Float64}(undef, NE, 4)
#    columns[:,1] = hcg[1]  # energy list
#    columns[:,2] = hcg[2]  # hot electron
#    columns[:,3] = hcg[3]  # hot hole
#    columns[:,4] = dist[2] # population

#    writedlm(outname, columns)

end

raw()
