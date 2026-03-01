include("../src/interfaces/nanoparticle.jl")


function run()

   # Properties of the nanoparticle (see documentation for options)
    shape  = "octahedron"
    mater  = "palladium"
    rad_nm = 2.5  # NP radius [nm]

    hb = hamiltonian_builder_predefined(shape, rad_nm, mater)

   # Select the calculation:  density of states
    N  = 200        # number of Chebyshev polynomials (default: 200)
    minE_Ha = -0.0  # minimum energy point [Ha] (default: -2) (dos printed in eV)
    maxE_Ha = +1.0  # maximum energy point [Ha] (default: +2) (dos printed in eV)
    dosname = "dos.dat"  # output file name (default: "dos.dat")

   # Builder for DOS
    dosb = calculation_builder_dos(N=N, minE_Ha=minE_Ha, maxE_Ha=maxE_Ha, name=dosname)

   # Run for DOS
    nanoparticle_run(ham_builder=hb, dos_builder=dosb)
end

run()
