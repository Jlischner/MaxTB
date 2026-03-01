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


# Script with full control over every parameter of the shape construction
function lowlevel_shp()

   # Properties of the nanoparticle (see documentation for options)
    shape = "cube"
    mater = "gold"
    outname = "positions.dat"

   # NP geometry [nm]
    length = 1.0
    height = 1.5
    width  = 2.0


   # Construct the TB model corresponding to the selected material
    onsite, first_neighbor, second_neighbor, A, B, fermi_Ha, a0, diel = tightbinding(mater)

   # Generate list of atomic positions
    Elist, Edict, R = generate_shape_FCC(shape, a0, length, height, width)
    print_positions(R, outname)
end

lowlevel_shp()
