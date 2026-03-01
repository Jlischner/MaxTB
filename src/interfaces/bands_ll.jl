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


# Script with full control over every parameter of the band structure calculation
function lowlevel_bs()

   # Properties of the nanoparticle (see documentation for options)
    mater = "gold"
    outname = "bands.dat"


   # Construct the TB model corresponding to the selected material
    onsite, first_neighbor, second_neighbor, A, B, fermi_Ha, a0, diel = tightbinding(mater)

   # Calculate the band structure along a predefined path
    bands = get_bands(onsite, first_neighbor, second_neighbor)


   # Save to file
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
end

lowlevel_bs()
