# Select materials
include("../materials/aluminum.jl")
include("../materials/copper.jl")
include("../materials/gold.jl")
include("../materials/iridium.jl")
include("../materials/palladium.jl")
include("../materials/platinum.jl")
include("../materials/rhodium.jl")
include("../materials/silver.jl")
include("../materials/cube3D.jl")


function tightbinding(material)
    if material == "aluminum"
        return aluminum()
    elseif material == "copper"
        return copper()
    elseif material == "gold"
        return gold()
    elseif material == "iridium"
        return iridium()
    elseif material == "palladium"
        return palladium()
    elseif material == "platinum"
        return platinum()
    elseif material == "rhodium"
        return rhodium()
    elseif material == "silver"
        return silver()
    elseif material == "3DTB"
        return cube()
    else
        @printf("\nMaterial not supported\n\n")
        return -1
    end
end
