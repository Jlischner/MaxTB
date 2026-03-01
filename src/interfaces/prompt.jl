include("nanoparticle.jl")


# Script with limited control over the simulation parameters:  hcg + dist + dos + band
function run_full()

   # Calculations to perform
    dos_run  = true
    band_run = true
    hcg_run  = true
    dist_run = true


   # Properties of the nanoparticle (see documentation for options)
    shape  = "cube"
    mater  = "gold"
    rad_nm = 1.5  # NP radius [nm]

    hb = hamiltonian_builder_predefined(shape, rad_nm, mater)


   # Properties of the external perturbation
    freq_eV = 2.4  # optical potential frequency [eV]
    eps_m   = 2.0  # relative dielectric constant of the environment (default: 1)
    eps     = -200+3im  # relative dielectric constant of the nanoparticle (default: prebuild)
    if shape == "sphere"
        pb = potential_builder_sphere(freq_eV, eps=eps, eps_m=eps_m)
    elseif shape == "import"
        pb = potential_builder_import(freq_eV, filename)
    else
        pb = potential_builder_sphere(freq_eV, eps=eps, eps_m=eps_m)
    end


   # Select the calculation:  density of states
    N  = 200        # number of Chebyshev polynomials (default: 200)
    NK = 4          # number of random vectors (default: 5)
    NE = 1000       # number of energy points (default: 1000)
    minE_Ha = -1.0  # minimum energy point [Ha] (default: -2) (dos printed in eV)
    maxE_Ha = +1.0  # maximum energy point [Ha] (default: +2) (dos printed in eV)
    dosname = "dos.dat"  # output file name (default: "dos.dat")

    if (dos_run)
        dosb = calculation_builder_dos(N=N, NE=NE, minE_Ha=minE_Ha, maxE_Ha=maxE_Ha, NK=NK, name=dosname)
    else
        dosb = []
    end


   # Select the calculation:  band structure
    bandname = "bands.dat"  # output file name (default: "bands.dat")

    if (band_run)
        bb = calculation_builder_bands(name=bandname)
    else
        bb = []
    end


   # Select the calculation:  hot carrier generation
    T     = 298.0  # temperature [K] (default: 298 K)
    gph   = 0.010  # lifetime [eV] (default: 10meV)
    conv  = true   # include convergence data in the output (default: false)
    N     = 200    # number of Chebyshev polynomials (default: 200)
    NR    = 5      # number of random vectors (default: 5)
    NB    = 10     # number of Chebyshev polynomial blocks (default: 10)
    NE    = 1000   # number of output energies (default: 1000)
    hname = "hcg.dat"  # hot carrier generation output filename (default: hcg.dat)

    if (hcg_run)
        gb = calculation_builder_hcg(temp=T, gph=gph, N=N, NR=NR, NB=NB, NE=NE, conv=conv, name=hname)
    else
        gb = []
    end


   # Select the calculation:  hot carrier distribution
    T     = 298.0  # temperature [K] (default: 298 K)
    conv  = true   # include convergence data in the output (default: false)
    N     = 200    # number of Chebyshev polynomials (default: 200)
    NR    = 5      # number of random vectors (default: 5)
    NB    = 10     # number of Chebyshev polynomial blocks (default: 10)
    NE    = 1000   # number of output energies (default: 1000)
    dname = "dist.dat"  # hot carrier generation output filename (default: dist.dat)


   # Relaxation rates [eV]
    function Gamma(E,Ep)
        return 0.100
    end

    if (dist_run)
        db = calculation_builder_dist(temp=T, Gamma=Gamma, N=N, NR=NR, NB=NB, NE=NE, conv=conv, name=dname)
    else
        db = []
    end


   # Run all
    nanoparticle_run(ham_builder=hb, pot_builder=pb, hcg_builder=gb, dist_builder=db, dos_builder=dosb, band_builder=bb)
end

run_full()
