using DataStructures
using Dates
using Interpolations
using LinearAlgebra
using MKL
using MKLSparse
using NearestNeighbors
using Printf
using Random
using SparseArrays

include("../materials/materials.jl")


function Deque_to_vec(Elist,Edict,vararg...;Transform=nothing)
    NN=length(Elist)
    R=zeros(Float64,NN,3)
    if isnothing(Transform)
        for Ri in Elist
            idx=Edict[Ri]
            R[idx,:]=Ri
        end
        return R
    else
        for Ri in Elist
            idx=Edict[Ri]
            R[idx,:]=Transform(Ri,vararg)
        end
    return R
    end
end


function getVAA(R, hops)
    V=Array{Float64,2}(undef,9,9)

   # Direction cosines
    l,m,n = R./norm(R)

   # basis ordering: s, dxy, dzx, dyz, dx2-y2, d3z2-r2, px, py, pz

    sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del = hops

   # s-s interaction
    V[1,1] = ss_sig  # s -> s

   # p-p interaction
    V[7,8] = l*m*pp_sig -     l*m*pp_pi  # px -> py
    V[7,9] = l*n*pp_sig -     l*n*pp_pi  # px -> pz
    V[8,9] = m*n*pp_sig -     m*n*pp_pi  # py -> pz
    V[7,7] = l^2*pp_sig + (1-l^2)*pp_pi  # px -> px
    V[8,8] = m^2*pp_sig + (1-m^2)*pp_pi  # py -> py
    V[9,9] = n^2*pp_sig + (1-n^2)*pp_pi  # pz -> pz
    V[8,7] = V[7,8]                      # py -> px
    V[9,7] = V[7,9]                      # pz -> px
    V[9,8] = V[8,9]                      # pz -> py

   # d1-d1 interaction
    V[2,3] = 3*l^2*m*n*dd_sig +            m*n*(1-4*l^2)*dd_pi +   m*n*(l^2-1)*dd_del  # dxy -> dzx
    V[3,4] = 3*n^2*m*l*dd_sig +            m*l*(1-4*n^2)*dd_pi +   l*m*(n^2-1)*dd_del  # dzx -> dyz
    V[2,4] = 3*l*m^2*n*dd_sig +            l*n*(1-4*m^2)*dd_pi +   l*n*(m^2-1)*dd_del  # dxy -> dyz
    V[2,2] = 3*l^2*m^2*dd_sig + ( l^2 + m^2 - 4*l^2*m^2)*dd_pi + (n^2+l^2*m^2)*dd_del  # dxy -> dxy
    V[3,3] = 3*l^2*n^2*dd_sig + ( l^2 + n^2 - 4*l^2*n^2)*dd_pi + (m^2+l^2*n^2)*dd_del  # dzx -> dzx
    V[4,4] = 3*m^2*n^2*dd_sig + ( m^2 + n^2 - 4*m^2*n^2)*dd_pi + (l^2+m^2*n^2)*dd_del  # dyz -> dyz
    V[3,2] = V[2,3]                                                                    # dzx -> dxy
    V[4,2] = V[2,4]                                                                    # dyz -> dxy
    V[4,3] = V[3,4]                                                                    # dyz -> dzx

   # d1-d2 interaction
    V[2,5] =                   (1.5*l*m*(l^2-m^2)*dd_sig +        2*l*m*(m^2-l^2)*dd_pi +         .5*l*m*(l^2-m^2)*dd_del)  # dxy -> dx2-y2
    V[2,6] =         (3^.5*l*m*(n^2-.5*(l^2+m^2))*dd_sig -         2*3^.5*l*m*n^2*dd_pi +      .5*3^.5*l*m*(1+n^2)*dd_del)  # dxy -> d3z2-r2
    V[3,5] =                   (1.5*n*l*(l^2-m^2)*dd_sig +    n*l*(1-2*(l^2-m^2))*dd_pi -     n*l*(1-.5*(l^2-m^2))*dd_del)  # dzx -> dx2-y2
    V[3,6] =         (3^.5*l*n*(n^2-.5*(l^2+m^2))*dd_sig + 3^.5*l*n*(l^2+m^2-n^2)*dd_pi -    .5*3^.5*l*n*(l^2+m^2)*dd_del)  # dzx -> d3z2-r2
    V[4,5] =                   (1.5*m*n*(l^2-m^2)*dd_sig -    m*n*(1+2*(l^2-m^2))*dd_pi +      m*n*(1+(l^2-m^2)/2)*dd_del)  # dyz -> dx2-y2
    V[4,6] =         (3^.5*m*n*(n^2-.5*(l^2+m^2))*dd_sig + 3^.5*m*n*(l^2+m^2-n^2)*dd_pi -    .5*3^.5*m*n*(l^2+m^2)*dd_del)  # dyz -> d3z2-r2
    V[5:6,2:4] = transpose(V[2:4,5:6])                                                                                      # transpose

   # d2-d2 interaction
    V[5,6] = .5*3^.5*(l^2-m^2)*(n^2-.5*(l^2+m^2))*dd_sig +     3^.5*n^2*(m^2-l^2)*dd_pi + 3^.5*(1+n^2)*(l^2-m^2)/4*dd_del  # dx2-y2 -> d3z2-r2
    V[5,5] =                      3/4*(l^2-m^2)^2*dd_sig +  (l^2+m^2-(l^2-m^2)^2)*dd_pi +      (n^2+(l^2-m^2)^2/4)*dd_del  # dx2-y2 -> dx2-y2
    V[6,6] =                 (n^2-.5*(l^2+m^2))^2*dd_sig +        3*n^2*(l^2+m^2)*dd_pi +          3/4*(l^2+m^2)^2*dd_del  # d3z2-r2 -> d3z2-r2
    V[6,5] = V[5,6]                                                                                                        # d3z2-r2 -> dx2-y2

   # s-d interaction
    V[1,2] =           3^.5*l*m*sd_sig  # s -> dxy
    V[1,3] =           3^.5*l*n*sd_sig  # s -> dzx
    V[1,4] =           3^.5*n*m*sd_sig  # s -> dyz
    V[1,5] =   3^.5/2*(l^2-m^2)*sd_sig  # s -> dx2-y2
    V[1,6] = (n^2-.5*(l^2+m^2))*sd_sig  # s -> d3z2-r2
    V[2:6,1] = transpose(V[1,2:6])      # transpose

   # s-p interaction
    V[1,7] = l*sp_sig                # s -> px
    V[1,8] = m*sp_sig                # s -> py
    V[1,9] = n*sp_sig                # s -> pz
    V[7:9,1] = -transpose(V[1,7:9])  # transpose

   # p-d interaction
    V[7,2] =         3^0.5*l^2*m*pd_sig +         m*(1-2*l^2)*pd_pi  # px -> dxy
    V[8,4] =         3^0.5*m^2*n*pd_sig +         n*(1-2*m^2)*pd_pi  # py -> dyz
    V[9,3] =         3^0.5*n^2*l*pd_sig +         l*(1-2*n^2)*pd_pi  # pz -> dzx

    V[7,4] =         3^0.5*l*m*n*pd_sig -             2*l*m*n*pd_pi  # px -> dyz
    V[8,3] =         3^0.5*l*m*n*pd_sig -             2*l*m*n*pd_pi  # py -> dzx
    V[9,2] =         3^0.5*l*m*n*pd_sig -             2*l*m*n*pd_pi  # pz -> dxy

    V[7,3] =         3^0.5*l^2*n*pd_sig +         n*(1-2*l^2)*pd_pi  # px -> dzx
    V[8,2] =         3^0.5*m^2*l*pd_sig +         l*(1-2*m^2)*pd_pi  # py -> dxy
    V[9,4] =         3^0.5*n^2*m*pd_sig +         m*(1-2*n^2)*pd_pi  # pz -> dyz

    V[7,5] = 3^0.5/2*l*(l^2-m^2)*pd_sig +       l*(1-l^2+m^2)*pd_pi  # px -> dx2-y2
    V[8,5] = 3^0.5/2*m*(l^2-m^2)*pd_sig -       m*(1+l^2-m^2)*pd_pi  # py -> dx2-y2
    V[9,5] = 3^0.5/2*n*(l^2-m^2)*pd_sig -         n*(l^2-m^2)*pd_pi  # pz -> dx2-y2
    V[7,6] = l*(n^2-(l^2+m^2)/2)*pd_sig -         3^0.5*l*n^2*pd_pi  # px -> d3z2-r2
    V[8,6] = m*(n^2-(l^2+m^2)/2)*pd_sig -         3^0.5*m*n^2*pd_pi  # py -> d3z2-r2
    V[9,6] = n*(n^2-(l^2+m^2)/2)*pd_sig +   3^0.5*n*(l^2+m^2)*pd_pi  # pz -> d3z2-r2
    V[2:6,7:9] = -transpose(V[7:9,2:6])                              # transpose

    return V
end


function get_Hk_FCC(kpoint, onsite, nearest, second)
   # 18 positions of the nearest (and next nearest) neighbors in cartesian coordinates?
    nnlist=[[1, 1, 0];;[-1, 1, 0];;[1,-1, 0];;[-1,-1, 0];;
            [1, 0, 1];;[-1, 0, 1];;[1, 0,-1];;[-1, 0,-1];;
            [0, 1, 1];;[ 0,-1, 1];;[0, 1,-1];;[ 0,-1,-1];;
            [2, 0, 0];;[-2, 0, 0];;[0, 2, 0];;[ 0,-2, 0];;[0,0,2];;[0,0,-2]]

    # Hk = Array{Float64,2}(undef,9,9)
    Hk = zeros(ComplexF64, 9, 9)

    for j=1:18
        dR    = nnlist[:,j]
        dist  = norm(dR)

        ph = dot(dR,kpoint)
        phase = exp(im*ph*pi)

        if dist < 1.99
            # println("dist<1.99")
            Hk += getVAA(dR, nearest).*phase
        else
            # println("dist>1.99")
            vv = getVAA(dR, second )
            # println("norm", norm(vv - vv'))

            Hk = Hk + vv.*phase
            # println(kpoint," ", j, " ", ph, " ", phase, " ", dR, " ", vv[1,1])
            # println("Diag of matrix:")
            # println(diag(Hk))
        end
    end
    # println("after orbitals:")
    # println(diag(Hk))

    for i=1:9
        Hk[i,i] += onsite[i]
    end

    return Hk
end


function get_bands(onsite, nearest, second)

   # Define the path in k-space
    G = [0.000, 0.000, 0.000]
    X = [0.000, 1.000, 0.000]
    W = [0.500, 1.000, 0.000]
    K = [0.750, 0.750, 0.000]
    L = [0.500, 0.500, 0.500]
    U = [0.250, 1.000, 0.250]

    # path = [G, X, W, L, G, K]
    path = [G, X, W, K, G, L, U, W, L]
    Npaths = length(path) - 1
    kpp = 30  # number of kpts per path
    NK = (kpp-1)*Npaths + 1
    klist = zeros(Float64, NK,  3)
    line  = zeros(Float64, kpp, 3)

   # Create k-paths
    for i in 1:Npaths
        Pi = path[i]
        Pf = path[i+1]
        # println(Pi, Pf)

        for d in 1:3
            line[:,d] = LinRange(Pi[d], Pf[d], kpp)
        end

        # println(line)

        n = (i-1)*(kpp-1) + 1
        klist[n:n+kpp-1,:] = line

    end
    klist[end,:] = path[end]
    # println(klist)

   # Calculate eigenenergies for each path
    bands = zeros(Float64, NK, 9+4)
    for k in 1:NK
        kvec = klist[k,:]
        hk = Hermitian(get_Hk_FCC(kvec, onsite, nearest, second))
        # bands[k,4:end-1] = diag(hk)
        bands[k,4:end-1] = eigvals(hk)
        bands[k,end] = norm(hk)
        # println("---", bands[k,:])
    end

    bands[:,1:3] = klist
    return bands
end



############LEONARDO#################


function read_xyz_positions(filename::AbstractString)
    # Read coordinates from a standard XYZ file and return an Nat x 3 matrix.
    lines = readlines(filename)
    if length(lines) < 2
        error("XYZ file $filename must contain at least two header lines")
    end

    Nat = try
        parse(Int, strip(lines[1]))
    catch
        error("First line of XYZ file $filename must be the number of atoms")
    end

    R = zeros(Float64, Nat, 3)
    parsed = 0

    for (offset, line) in enumerate(lines[3:end])
        line_number = offset + 2
        parsed == Nat && break

        l = strip(line)
        isempty(l) && continue
        fields = split(l)

        # Accept either "Au x y z" or "x y z".
        first_as_float = tryparse(Float64, fields[1])
        coord_start = isnothing(first_as_float) ? 2 : 1
        if length(fields) < coord_start + 2
            error("Malformed XYZ line $line_number in $filename: expected 3 coordinates")
        end

        x = try
            parse(Float64, fields[coord_start])
        catch
            error("Malformed x coordinate at XYZ line $line_number in $filename")
        end
        y = try
            parse(Float64, fields[coord_start + 1])
        catch
            error("Malformed y coordinate at XYZ line $line_number in $filename")
        end
        z = try
            parse(Float64, fields[coord_start + 2])
        catch
            error("Malformed z coordinate at XYZ line $line_number in $filename")
        end

        parsed += 1
        R[parsed, 1] = x
        R[parsed, 2] = y
        R[parsed, 3] = z
    end

    if parsed != Nat
        error("XYZ file $filename declares $Nat atoms but only $parsed were parsed")
    end

    return R
end


function infer_neighbor_shells(R::Array{Float64,2}; k_probe::Int64=24, shell_split::Float64=0.08)
    # Infer first and second neighbor distances from local nearest-neighbor statistics.
    Nat = size(R, 1)
    if Nat < 3
        error("Need at least 3 atoms to infer first and second neighbor shells")
    end

    points = permutedims(R) # 3 x Nat
    tree = KDTree(points)
    k = min(k_probe, Nat)
    _, dists = knn(tree, points, k, true)

    first_local = Float64[]
    second_local = Float64[]
    epsd = 1e-12

    for i in 1:Nat
        di = dists[i]
        d1 = 0.0
        d2 = 0.0
        for d in di
            if d <= epsd
                continue
            end
            if d1 <= epsd
                d1 = d
            elseif d > d1 * (1.0 + shell_split)
                d2 = d
                break
            end
        end

        if d1 > epsd
            push!(first_local, d1)
        end
        if d2 > epsd
            push!(second_local, d2)
        end
    end

    if isempty(first_local) || isempty(second_local)
        error("Could not infer first and second neighbor shells from coordinates")
    end

    sort!(first_local)
    sort!(second_local)
    mid1 = length(first_local) ÷ 2 + 1
    mid2 = length(second_local) ÷ 2 + 1

    return first_local[mid1], second_local[mid2]
end


function slater_koster_general(R::Array{Float64,2}, onsite, first_neighbour, second_neighbour, A, B;
                               first_shell=nothing, second_shell=nothing, shell_tol::Float64=0.12,
                               Norbitals::Int64=9, position_scale::Float64=1.0)
    # Slater-Koster for arbitrary coordinates (non-FCC):
    # neighbor shells are identified by distance and missing neighbors are simply omitted.
    Nat, ndim = size(R)
    if ndim != 3
        error("R must have size Nat x 3")
    end
    if Nat == 0
        error("R cannot be empty")
    end

    d1 = isnothing(first_shell) ? nothing : Float64(first_shell)
    d2 = isnothing(second_shell) ? nothing : Float64(second_shell)
    if isnothing(d1) || isnothing(d2)
        d1i, d2i = infer_neighbor_shells(R)
        d1 = isnothing(d1) ? d1i : d1
        d2 = isnothing(d2) ? d2i : d2
    end
    if d1 <= 0 || d2 <= 0 || d2 <= d1
        error("Invalid shell distances: first_shell=$d1 second_shell=$d2")
    end

    low1 = d1 * (1.0 - shell_tol)
    high1 = d1 * (1.0 + shell_tol)
    low2 = d2 * (1.0 - shell_tol)
    high2 = d2 * (1.0 + shell_tol)
    if high1 >= low2
        error("Neighbor shell windows overlap, lower shell_tol or set explicit shell distances")
    end

    onsiteKPM = (onsite .- A) ./ B

    iidx = Int64[]
    jidx = Int64[]
    value = ComplexF64[]

    iidx_r = Int64[]
    value_r = ComplexF64[]

    function addV!(iidx, jidx, value, V, iptr, jptr)
        local a, b = size(V)
        for i = 1:a
            for j = 1:b
                vv = V[i, j]
                if abs(vv) > 1e-10
                    push!(iidx, iptr + i)
                    push!(jidx, jptr + j)
                    push!(value, vv)
                end
            end
        end
    end

    # Onsite terms and position operator, mirroring the FCC workflow.
    for Atomi = 1:Nat
        iptr = (Atomi - 1) * Norbitals
        for k = 1:Norbitals
            push!(iidx, k + iptr)
            push!(jidx, k + iptr)
            push!(value, onsiteKPM[k])

            push!(iidx_r, k + iptr)
            push!(value_r, R[Atomi] * position_scale)
        end
    end

    points = permutedims(R) # 3 x Nat
    tree = KDTree(points)

    # Build hoppings by distance shells; surface atoms naturally have fewer neighbors.
    for Atomi = 1:Nat
        neigh = inrange(tree, points[:, Atomi], high2, true)
        iptr = (Atomi - 1) * Norbitals
        for Atomj in neigh
            Atomj <= Atomi && continue

            Rij = vec(points[:, Atomj] .- points[:, Atomi])
            dij = norm(Rij)

            hops = nothing
            if low1 <= dij <= high1
                hops = first_neighbour
            elseif low2 <= dij <= high2
                hops = second_neighbour
            else
                continue
            end

            jptr = (Atomj - 1) * Norbitals
            Vij = getVAA(Rij, hops) ./ B
            Vji = getVAA(-Rij, hops) ./ B
            addV!(iidx, jidx, value, Vij, iptr, jptr)
            addV!(iidx, jidx, value, Vji, jptr, iptr)
        end
    end

    N = Nat * Norbitals
    H = sparse(iidx, jidx, value, N, N)
    pos_operator = sparse(iidx_r, iidx_r, value_r, N, N)
    v = (H * pos_operator - pos_operator * H) .* im
    return H, v
end


function slater_koster_general(xyz_file::AbstractString, onsite, first_neighbour, second_neighbour, A, B; kwargs...)
    R = read_xyz_positions(xyz_file)
    return slater_koster_general(R, onsite, first_neighbour, second_neighbour, A, B; kwargs...)
end


function slater_koster_FCC(Elist, Edict, onsite, first_neighbour, second_neighbour, A, B; L1=999999, L2=999999, L3=999999, periodic=false)
    # Use the Slater-Koster parametrization for a set of atomic positions defined in Elist and Edict, 
    # using the parameters onsite, first_neighbour, second_neighbour
    # This is only valid for the FCC lattice

   # Get the TB parameters
   # onsite, first_neighbor, second_neighbor, A, B, fermi, diel_fun = tightbinding(material)
   # println(onsite, first_neighbor, second_neighbor)

   # 18 positions of the nearest (and next nearest) neighbors in cartesian coordinates?
    nnlist=[[1, 1, 0];;[-1, 1, 0];;[1,-1, 0];;[-1,-1, 0];;
            [1, 0, 1];;[-1, 0, 1];;[1, 0,-1];;[-1, 0,-1];;
            [0, 1, 1];;[ 0,-1, 1];;[0, 1,-1];;[ 0,-1,-1];;
            [2, 0, 0];;[-2, 0, 0];;[0, 2, 0];;[ 0,-2, 0];;[0,0,2];;[0,0,-2]]


    function getHamiltonian(Elist::Deque{Vector{Int64}}, Edict::Dict{Vector{Int64},Int64},
                            nnlist::Array{Float64,2}, nnlist_for_table::Array{Int64,2},
                            getHopping::Function;
                            Norbitals::Int64=9, A::Float64=0.0, B::Float64=1.0, L1::Int64=999999, L2::Int64=999999, L3::Int64=999999, periodic=false)

       # Elist  - list of all the positions of atoms inside the nanoparticle [(x1,y1,z1), (x2,y2,z2), ...]
       # Edict  - dictionary with all the sites inside the nanoparticle (x,y,z): index
       # nnlist - list of vectors to nearest all nearest neighbors (half of list below)
       # nnlist_for_table - list of vectors to nearest all nearest neighbors
       # A - KPM shift
       # B - KPM scale

        NNA = length(Elist)  # number of atoms
        onsiteKPM = (onsite.-A)./B  # turning Ha into KPM units

        function addV!(iidx,jidx,value,V,iptr,jptr)
           # Add the values of the V 9x9 matrix into a sparse matrix format
           # iidx (jidx) is the set of (lines,columns) of the sparse matrix
           # iptr, jptr index the atoms
            local i = 0
            local j = 0
            local vv = 0.0
            local a,b = size(V)
            for i = 1:a
                for j = 1:b
                    vv = V[i,j]
                    if abs(vv) > 1e-10
                        push!(iidx, iptr+i)
                        push!(jidx, jptr+j)
                        push!(value, vv)
                    end
                end
            end
        end


       # Build a neighbor table for every particle
       # [[n1, n2, ..., n9],     <- indices of neighbors of atom 1
       # [[n1, -1, ..., n9],     <- indices of neighbors of atom 2. -1 means it was not found
       # ...
        global table=Array{Int64,2}(undef, NNA, 18)
        for Ri in Elist
            Idx = Edict[Ri]
            for j = 1:18
                Rij = nnlist_for_table[:,j]
                Rj  = Ri + Rij
                # println(Rj,L1,L2,L3)

               # if it is periodic, then there are no negative positions, by construction
                if periodic
                    Rj[1] = (Rj[1]+L1) % L1
                    Rj[2] = (Rj[2]+L2) % L2
                    Rj[3] = (Rj[3]+L3) % L3
                end

                # println(Rj)
                Idj = get(Edict, Rj, -1)
                table[Idx,j] = Idj
            end
        end

        # println("table\n")
        # display(table)

        R = Deque_to_vec(Elist,Edict)

       # build the sparse Hamiltonian matrix from the neighbor table and the hoppings
        global iidx  = Array{Int64,1}([])
        global jidx  = Array{Int64,1}([])
        global value = Array{ComplexF64,1}([])  # Hamiltonian

        global iidx_r  = Array{Int64,1}([])
        global value_r = Array{ComplexF64,1}([])  # position operator
        for Atomi=1:NNA
            iptr = (Atomi-1)*Norbitals

           # onsite terms
            for k=1:Norbitals
                push!(iidx, k+iptr)
                push!(jidx, k+iptr)
                push!(value, onsiteKPM[k])

                push!(iidx_r, k+iptr)
                push!(value_r, R[Atomi])
            end


            for j=1:18
                Atomj = table[Atomi,j]
                if Atomj >= 0
                   # Neighboring atom is atomA
                    jptr = (Atomj-1)*Norbitals
                    Rij  = nnlist[:,j]
                    V    = getHopping(Rij)./B
                    addV!(iidx, jidx, value, V, iptr, jptr)
                end
            end
        end

        # mat = zeros(Float64,length(iidx), 3)
        # mat[:,1] = iidx
        # mat[:,2] = jidx
        # mat[:,3] = value
        # display(mat)
        H = sparse(iidx, jidx, value)
        pos_operator = sparse(iidx_r, iidx_r, value_r)
        v = (H*pos_operator - pos_operator*H).*im
        return H,v
    end


    function getVAA2(R)
        V=Array{Float64,2}(undef,9,9)
        function getHoppingParameters(d)
            if d<0.9999999
                return first_neighbor
            else
                return second_neighbor
            end
        end

        hops = getHoppingParameters(norm(R))
        return getVAA(R, hops)
    end

    H,v = getHamiltonian(Elist, Edict, nnlist./2, nnlist, getVAA2, A=A, B=B, L1=L1, L2=L2, L3=L3, periodic=periodic)
    return H,v
end
