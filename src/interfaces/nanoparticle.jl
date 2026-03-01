using DelimitedFiles
using Printf

pathw = "../"
include(pathw*"libs/aux.jl")
include(pathw*"libs/cheb.jl")
include(pathw*"libs/dist.jl")
include(pathw*"libs/dos.jl")
include(pathw*"libs/hcg.jl")
include(pathw*"libs/potential.jl")
include(pathw*"libs/shape_lib.jl")
include(pathw*"libs/slater_koster.jl")
include(pathw*"libs/sumcheb.jl")


# BUILDER HAMILTONIAN PREDEFINED
function hamiltonian_builder_predefined(shape::AbstractString, rad_nm::Float64, material::AbstractString)
    #
   # Check if the material selected is implemented
    materials_implemented = ["aluminum", "copper", "gold", "iridium", "palladium", "platinum", "rhodium", "silver"]
    if !(material in materials_implemented)
        @printf("\nMaterial %s not implemented yet. EXITING\n\n", material) ; flush(stdout)
        exit()
    end
    #
   # Check if the shape selected is implemented
    shapes_implemented = ["cube", "dodecahedron", "octahedron", "sphere"]
    if !(shape in shapes_implemented)
        @printf("\nShape %s not implemented yet. EXITING\n\n", shape) ; flush(stdout)
        exit()
    end
    #
   # Check if the radius makes sense
    if rad_nm <= 0
        @printf("\nNP radius must be positive. EXITING\n\n") ; flush(stdout)
    end
    #
    description = "predefined"
    #
    return [description, shape, rad_nm, material]
end


# BUILDER POTENTIAL SPHERE
function potential_builder_sphere(freq_eV::Float64; eps = "default", eps_m::Float64 = 1.0)
    #
    if freq_eV <= 0.0
        @printf("\nThe frequency must be positive. EXITING\n\n") ; flush(stdout)
        exit()
    end
    #
    description = "sphere"
    #
    return [description, freq_eV, eps, eps_m]
end


# BUILDER POTENTIAL IMPORTED
function potential_builder_import(freq_eV::Float64, filename::AbstractString)
    #
    if freq_eV <= 0.0
        @printf("\nThe frequency must be positive. EXITING\n\n") ; flush(stdout)
        exit()
    end
    #
    description = "import"
    #
    return [description, freq_eV, filename]
end


# BUILDER DOS
function calculation_builder_dos(;NE::Int64=1000, minE_Ha::Float64=-2, maxE_Ha::Float64=2, N::Int64=200, NK::Int=5, name::AbstractString="dos.dat")
    #
    description = "dos"
    #
    return [description, NE, minE_Ha, maxE_Ha, N, NK, name]
end


# BUILDER BANDS
function calculation_builder_bands(;name::AbstractString="bands.dat")
    #
    description = "band structure"
    #
    return [description, name]
end


# BUILDER HCG
function calculation_builder_hcg(;temp::Float64=298.0, gph::Float64=0.01, N::Int64=200, NR::Int64=5, NB::Int64=10, NE::Int64=1000, conv::Bool=false, name::AbstractString="hcg.dat")
    #
   # Check if the quantities make sense
    quantities = [temp, gph, N, NR, NB, NE]
    varnames = ["Temperature", "Lifetime", "Number of Chebyshev polynomials", "Number of random vectors", "Number of blocks", "Number of output energies"]
    #
    for i in 1:6
        q = quantities[i]
        v = varnames[i]
        if q<0
            @printf("\n%s cannot be negative. EXITING\n\n", v) ; flush(stdout)
            exit()
        end
    end
    #
    description = "hcg"
    #
    return [description, temp, gph, N, NR, NB, NE, conv, name]
end


# BUILDER DISTRIBUTION
function calculation_builder_dist(;temp::Float64=298, Gamma=default, N::Int64=200, NR::Int64=5, NB::Int64=10, NE::Int64=1000, conv::Bool=false, name::AbstractString="dist.dat")
    #
   # Check if the quantities make sense
    quantities = [temp, N, NR, NB, NE]
    varnames = ["Temperature", "Number of Chebyshev polynomials", "Number of random vectors", "Number of blocks", "Number of output energies"]
    #
    for i in 1:5
        q = quantities[i]
        v = varnames[i]
        if q<0
            @printf("\n%s cannot be negative. EXITING\n\n", v) ; flush(stdout)
            exit()
        end
    end
    #
    description = "dist"
    #
    return [description, temp, Gamma, N, NR, NB, NE, conv, name]
end


# DEFAULT
function default(E,Ep)
    return 1.0
end




# Perform the complete workflow:  hcg + dist + dos + band
function nanoparticle_run(;ham_builder=[], pot_builder=[], hcg_builder=[], dist_builder=[], dos_builder=[], band_builder=[])
    #
    @printf("\n\n >>>  MAXTB PROGRAM STARTED  --  %s  <<<\n\n",Dates.format(now(), "HH:MM:SS / yyyy-mm-dd")) ; flush(stdout)
    t0 = time()
    #
   # Check which builders have been specified
    ham_b = length(ham_builder)  > 0
    pot_b = length(pot_builder)  > 0
    hcg_b = length(hcg_builder)  > 0
    dis_b = length(dist_builder) > 0
    dos_b = length(dos_builder)  > 0
    ban_b = length(band_builder) > 0
    #
   # Check what needs to be done
    need_mat = hcg_b || dis_b || dos_b || ban_b
    need_at  = hcg_b || dis_b || dos_b 
    need_H   = hcg_b || dis_b || dos_b 
    need_V   = hcg_b || dis_b
    #
    if need_H && !ham_b
        @printf("\nHamiltonian required for computation but no Hamiltonian builder specified. EXITING\n\n") ; flush(stdout)
        exit()
    end
    #
    if need_V && !pot_b
        @printf("\nPotential required for computation but no potential builder specified. EXITING\n\n") ; flush(stdout)
        exit()
    end
    #
    if need_mat
        mater = ham_builder[4]
       # Fetch the tight-binding model parameters
        onsite, first_neighbor, second_neighbor, A, B, fermi_Ha, a0, diel = tightbinding(mater)
    end
    #
    if need_at
        shape  = ham_builder[2]
        rad_nm = ham_builder[3]
       # Generate list of atomic positions
        Elist, Edict, R = generate_shape_FCC(shape, a0, rad_nm)
        @printf("\nNumber of atoms:  %d\n", length(R)) ; flush(stdout)
    end
    #
   # Get Hamiltonian from the predefined set
    if need_H
        ham_desc = ham_builder[1]
        #
        if ham_desc == "predefined"
           # Use list of atomic positions to determine the Hamiltonian and velocity operator
            H, v = slater_koster_FCC(Elist, Edict, onsite, first_neighbor, second_neighbor, A, B)
        #
       # Get Hamiltonian sparse matrix from file
        elseif ham_desc == "custom"
            ham_fname = ham_builder[2]
            # H, A, B, fermi = get_hami(ham_fname)
        #
        else
            @printf("\nRequested Hamiltonian builder is not supported. EXITING\n\n") ; flush(stdout)
            exit()
        end
    end
    #
   # Get the potential
    if need_V
       # Description and frequency of the optical potential (always required)
        pot_desc = pot_builder[1]
        freq     = pot_builder[2]
        #
       # Dielectric sphere:  pot_builder = [pot_desc, freq, eps, eps_m]
        if pot_desc == "sphere"
            eps = pot_builder[3]
            println("Input:  ", pot_builder) ; flush(stdout)
            if eps == "default"
                eps = diel(freq)
                @printf("Diel. constant:  %.3f + %.3fi\n", real(eps), imag(eps)) ; flush(stdout)
            end
            eps_m = pot_builder[4]  # dielectric constant of the environment
            Phi   = potential_sphere(R, eps, eps_m)
        #
       # Reading from file:  pot_builder = [pot_desc, freq, filename]
        elseif pot_desc == "import"
            pot_fname = pot_builder[3]
            Phi = comsol_read(pot_fname)
        #
        else
            @printf("Requested potential builder is not supported. EXITING\n\n") ; flush(stdout)
            exit()
        end
    end



# Computation of the hot-carrier generation rate
    if hcg_b
        @printf("\n\n><  Computing the hot carrier generation rate  ><\n") ; flush(stdout)
        hcg_desc = hcg_builder[1]
       # description, temperature, lifetime, number of Chebyshev polynomials,
       # number of random vectors, number of blocks, convergence, name
        desc, tempr, gph, N, NR, NB, NE, conv, name = hcg_builder
        #
       # Number of Chebyshev polynomials per block
        NTx = NB
        NTy = NB
        #
       # Number of blocks
        NNL = N÷NTx + 1
        NNR = N÷NTy + 1
        #
       # Build the Chebyshev matrix
        mumn = compute_mumn!(H, Phi, NNL, NNR, NTx, NTy, NR)
        #
       # Get inverse temperature
        kbT   = tempr*boltzmann  # [eV]
        beta  = 1.0/kbT
        kappa = 0.0  # parameter in the broadening function
        #
        if conv
            print_conv = true
            hcg = hcg_conv(mumn, A, B, NE, fermi_Ha, freq, beta, kappa, gph, print_conv)
        else
            hcg = hcg_auto(mumn, A, B, NE, fermi_Ha, freq, beta, kappa, gph)
        end
        #
       # Write output hcg data
        open(name, "w") do io
            for i in 1:size(hcg,1)
                for j in 1:size(hcg,2)
                    @printf(io,"%12.8f ",hcg[i,j])
                end
                Base.write(io,"\n")
            end
        end
        @printf("    |CPU time: %.0f s|\n", time() - t0) ; flush(stdout)
    end


# Computation of the hot-carrier distribution
    if dis_b
        @printf("\n\n><  Computing the population distribution  ><\n") ; flush(stdout)
        dist_desc = dist_builder[1]
       # description, temperature, lifetime, number of Chebyshev polynomials,
       # number of random vectors, number of blocks, convergence, name
        desc, tempr, Gamma, N, NR, NB, NE, conv, name = dist_builder
        #
       # Number of Chebyshev polynomials per block
        NTx = NB
        NTy = NB
        #
       # Number of blocks
        NNL = N÷NTx + 1
        NNR = N÷NTy + 1
        #
       # Build the Chebyshev matrix
        mumn = compute_mumn!(H, Phi, NNL, NNR, NTx, NTy, NR)
        #
       # Get inverse temperature
        kbT  = tempr*boltzmann  # [eV]
        beta = 1.0/kbT
        #
       # Number of integration points in the energy discretization
        Nx = N*2 + 1
        Ny = Nx
        perc1 = 100
        perc2 = 100
        Ei,Ej,opt = resum_mu_dd(mumn, Nx, Ny, A, B, perc1, perc2)
        NEp = Nx
        Gammaf = 1
        write=false
        #
        dist = sum_dist(opt, Ei, Ej, NE, NEp, fermi_Ha, freq, beta, Gammaf, Gamma, write)
        #
       # Write output dist data
        col1 = collect(dist[1])
        col2 = collect(dist[2])
        dist2 = hcat(col1,col2)
        open(name, "w") do io
            for i in 1:size(dist2,1)
                for j in 1:size(dist2,2)
                    @printf(io,"%12.8f ",dist2[i,j])
                end
            Base.write(io,"\n")
            end
        end
        @printf("    |CPU time: %.0f s|\n", time() - t0) ; flush(stdout)
    end


# Computation of the density of states
    if dos_b
        @printf("\n\n><  Computing the density of states  ><\n") ; flush(stdout)
       # Description, number of output energies, minimum energy point [Ha], maximum energy point [Ha]
       # number of polynomials, number of random vectors, name of output file
        description, NE, minE_Ha, maxE_Ha, N, Nk, outname = dos_builder
        #
        flag = 2  # f=2 -> random vectors are used
        perc = 100
        #
       # Build the Chebyshev vector
        mu = doscompute_mu!(H, N, Nk, flag)
        #
       # Resum the Chebyshev vector into the DOS
        dos = get_dos(mu, perc, minE_Ha, maxE_Ha, NE, A, B)
        #
       # Write output DOS data
        open(outname, "w") do io
            for i in 1:size(dos,1)
                for j in 1:size(dos,2)
                    @printf(io,"%12.8f ",dos[i,j])
                end
                Base.write(io,"\n")
            end
        end
        @printf("    |CPU time: %.0f s|\n", time() - t0) ; flush(stdout)
    end


# Computation of the band structure
    if ban_b
        @printf("\n\n><  Computing the band structure  ><\n") ; flush(stdout)
       # Description and name of output file
        description, outname = band_builder
        #
       # Calculate the band structure along a predefined path
        bands = get_bands(onsite, first_neighbor, second_neighbor)
        #
       # Write output band data
        open(outname, "w") do io
            for i in 1:size(bands,1)
                for j in 1:size(bands,2)
                    if j <= 3
                        @printf(io,"%8.6f ",bands[i,j])
                    else
                        @printf(io,"%12.8f ",bands[i,j]*Ha2eV)  # [eV]
                    end
                end
                Base.write(io,"\n")
            end
        end
        @printf("    |CPU time: %.0f s|\n", time() - t0) ; flush(stdout)
    end
    #
    @printf("\n\n >>>>  MAXTB PROGRAM ENDED  --  %s  <<<<\n\n\n",Dates.format(now(), "HH:MM:SS / yyyy-mm-dd")) ; flush(stdout)
    #
end
