using Interpolations


function copper()
   # Cu Slater-Koster 2-center parameters (approximation)
   # -> taken from pp.182 of Papaconstantopoulos' Handbook of the Band Structure of Elemental Solids, 2nd ed. (2015)
   # All quantities are in Rydberg, before being converted to Hartree
    Ry2eV = 13.605704
   # Local energies at each orbital [Ry]
    Es  = 0.79466
    Ep  = 1.35351
    Ed1 = 0.37307
    Ed2 = 0.37180
    onsite = [Es,Ed1,Ed1,Ed1,Ed2,Ed2,Ep,Ep,Ep]

   # First neighbor [Ry]
    ss_sig = -0.07518
    pp_sig =  0.19669
    pp_pi  =  0.01940
    dd_sig = -0.02566
    dd_pi  =  0.01800
    dd_del = -0.00408
    sp_sig =  0.11571
    sd_sig = -0.03107
    pd_sig = -0.03289
    pd_pi  =  0.01753
    first_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Second neighbor [Ry]
    ss_sig = -0.00092
    pp_sig =  0.05389
    pp_pi  =  0.00846
    dd_sig = -0.00451
    dd_pi  =  0.00241
    dd_del = -0.00029
    sp_sig =  0.01221
    sd_sig = -0.00852
    pd_sig = -0.00536
    pd_pi  =  0.00321
    second_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Convert the SK parameters from Ry to Ha
    onsite = onsite./2
    first_neighbor = first_neighbor./2
    second_neighbor = second_neighbor./2



   # Cu optical properties (experimental)
   # -> taken from pp.2203 (12-128) of Hayne's CRC Handbook of Chemistry and Physics, 97th ed. (2016)
   # Dielectric function calculated from the index of refraction n and the extinction coefficient k

   # Frequency Cu [eV]
    frequencies = [  0.10,     0.50,    1.00,    1.50,    1.70,    1.75,    1.80,    1.85,    1.90,    2.00,
                     2.10,     2.20,    2.30,    2.40,    2.60,    2.80,    3.00,    3.20,    3.40,    3.60,
                     3.80,     4.00,    4.20,    4.40,    4.60,    4.80,    5.00,    5.20,    5.40,    5.60,
                     5.80,     6.00,    6.50,    7.00,    7.50,    8.00,    8.50,    9.00,    9.50,   10.00,
                    11.00,    12.00,   13.00,   14.00,   14.50,   15.00,   15.50,   16.00,   17.00,   18.00,
                    19.00,    20.00,   21.00,   22.00,   23.00,   24.00,   25.00,   26.00,   27.00,   28.00,
                    29.00,    30.00,   31.00,   32.00,   33.00,   34.00,   35.00,   36.00,   37.00,   38.00,
                    39.00,    40.00,   41.00,   42.00,   43.00,   44.00,   45.00,   46.00,   47.00,   48.00,
                    49.00,    50.00,   51.00,   52.00,   53.00,   54.00,   55.00,   56.00,   57.00,   58.00,
                    59.00,    60.00,   61.00,   62.00,   63.00,   64.00,   65.00,   66.00,   67.00,   68.00,
                    69.00,    70.00,   75.00,   80.00,   85.00,   90.00]

   # Real part of dielectric constant of bulk Cu:  eps_real = n^2 - k^2
    eps_real = [-4240.769, -307.893, -71.717, -27.600, -19.577, -18.018, -16.278, -14.774, -13.425, -10.425,
                   -7.675,   -6.071,  -5.627,  -5.506,  -4.928,  -4.201,  -3.492,  -2.772,  -2.190,  -1.781,
                   -1.481,   -1.163,  -0.673,  -0.470,  -0.479,  -0.583,  -1.008,  -1.336,  -1.530,  -1.635,
                   -1.579,   -1.447,  -0.955,  -0.499,  -0.188,   0.000,   0.101,   0.215,   0.304,   0.409,
                    0.582,    0.655,   0.648,   0.605,   0.543,   0.516,   0.484,   0.454,   0.444,   0.479,
                    0.514,    0.572,   0.642,   0.702,   0.747,   0.785,   0.762,   0.686,   0.630,   0.617,
                    0.633,    0.672,   0.717,   0.744,   0.766,   0.788,   0.806,   0.810,   0.810,   0.833,
                    0.836,    0.836,   0.858,   0.858,   0.861,   0.880,   0.880,   0.880,   0.883,   0.883,
                    0.883,    0.886,   0.886,   0.886,   0.907,   0.907,   0.907,   0.910,   0.910,   0.910,
                    0.929,    0.929,   0.929,   0.929,   0.912,   0.912,   0.931,   0.931,   0.933,   0.933,
                    0.933,    0.933,   0.952,   0.952,   0.933,   0.915]

   # Imaginary part of dielectric constant of bulk Cu:  eps_imag = 2nk
    eps_imag = [ 4249.827,   60.295,   7.462,   2.735,   1.949,   1.785,   1.697,   1.694,   1.541,   1.750,
                    2.641,    4.316,   5.387,   5.824,   5.750,   5.522,   5.216,   5.092,   4.953,   4.899,
                    4.851,    4.610,   4.658,   4.887,   5.077,   5.233,   5.233,   4.968,   4.557,   4.106,
                    3.674,    3.307,   2.630,   2.328,   2.180,   2.122,   2.019,   1.895,   1.792,   1.706,
                    1.605,    1.591,   1.555,   1.526,   1.483,   1.434,   1.352,   1.273,   1.128,   0.997,
                    0.898,    0.792,   0.738,   0.699,   0.696,   0.710,   0.768,   0.736,   0.669,   0.602,
                    0.510,    0.447,   0.422,   0.392,   0.378,   0.364,   0.368,   0.350,   0.350,   0.335,
                    0.316,    0.316,   0.301,   0.301,   0.282,   0.285,   0.285,   0.285,   0.266,   0.266,
                    0.266,    0.247,   0.247,   0.247,   0.230,   0.230,   0.230,   0.211,   0.211,   0.211,
                    0.213,    0.213,   0.213,   0.213,   0.192,   0.192,   0.194,   0.194,   0.175,   0.175,
                    0.175,    0.175,   0.176,   0.176,   0.175,   0.154]

   # Interpolation of the optical frequencies
    eps_real_interp = LinearInterpolation(frequencies, eps_real)
    eps_imag_interp = LinearInterpolation(frequencies, eps_imag)
   #
    function dielectric(hw)
        return round(eps_real_interp(hw), digits=3) + 1im*round(eps_imag_interp(hw), digits=3)
    end


   # FCC lattice constant [nm] - 6.822 bohr
    a0 = 0.3610

   # KPM shift and scale
    min_bandval = -9.435800./Ry2eV./2  # min, no fermi [Ha]
    max_bandval = 25.260141./Ry2eV./2  # max, no fermi [Ha]
    fermi = 0.5805./2  # Fermi energy [Ha]

    A = (max_bandval-min_bandval)./2 + min_bandval + fermi
    B = (max_bandval-min_bandval)./2 * 1.05
    println("A = ", round(A, digits=3))
    println("B = ", round(B, digits=3))


    return [onsite, first_neighbor, second_neighbor, A, B, fermi, a0, dielectric]
end
