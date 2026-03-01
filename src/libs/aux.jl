boltzmann = 8.617333e-5   # Boltzmann constant [eV/K]
Ha2eV = 27.211386245988   # Energy conversion ratio [eV/Ha]


function meshgrid(x,y)
    local Nx = length(x)
    local Ny = length(y)
    xx = zeros(Float64, Nx, Ny)
    yy = zeros(Float64, Nx, Ny)
    for i=1:Nx
        xx[i,:].= x[i]
    end
    for j=1:Ny
        yy[:,j].= y[j]
    end 
    return xx, yy
end     


# Jackson kernel
function jackson_kernel(n,N)
    a = pi/(N+1)
    b = (N+1-n)*cos(a*n)
    c = sin(a*n)/tan(a)
    return (b+c)/(N+1)
end


# Build array of N Jackson kernel coefficients
function jackson_coefs(N)
    J = zeros(N,1)
    for i in 1:N
        J[i] = jackson_kernel(i-1, N)
    end
    J[1] = J[1]*0.5
    return J
end


# Fermi dirac distribution
function fermi_dirac(x,beta,mu)
    arg = beta*(x-mu)
    return 1.0/(1+exp(arg))
end


# # Read the atomic positions
# function read_positions(filename)
#    fr = open(filename,"r")
#   # Ignore the first 8 lines (header) to output file
#    for i=1:9
#        line = readline(fr)
#    end
#
#   # Iterate over the positions and save to array
#    flines=readlines(fr)
#    for i in flines
#        j = split(i)
#        x = parse(Float64,j[1])
#        y = parse(Float64,j[2])
#        z = parse(Float64,j[3])
#    end
#    close(fr)
#    return R
# end


# Exact diagonalization
function exact(H, Phi, Fermi, beta)
    println("Exact diagonalization")

    H_dens = Matrix(H)
    Phi_dens = Matrix(Phi)
    a,b = size(H_dens)
    vals, vecs = eigen(H_dens)
   # Testing vecs
    # n = 2
    # no = norm(Hdens*vecs[:,n] - vals[n]*vecs[:,n])
    # println("norm", no)
    println(typeof(Phi_dens), typeof(vecs))
    phi_eigen = vecs'*Phi_dens*vecs
    Gamma = zeros(a,a)

    # function fermi(x)
    #     return 1.0/(1.0 + exp(beta*(x-Fermi)))
    # end

    for i in 1:a
        Ei = vals[i]
        for j in 1:a
            Ej = vals[j]
            Gamma[i,j] = abs(phi_eigen[j,i])^2 #*fermi(Ei)*(1-fermi(Ej))
        end
    end
    # open("matrix.txt", "w") do file
    #     write(file, Gamma)
    # end
    writedlm("Phi_nm_jl.dat", Gamma)
    writedlm("eigen_energies_jl.dat", vals)
end
