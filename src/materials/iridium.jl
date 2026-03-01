using Interpolations


function iridium()
   # Ir Slater-Koster 2-center parameters (approximation)
   # -> taken from pp.288 of Papaconstantopoulos' Handbook of the Band Structure of Elemental Solids, 2nd ed. (2015)
   # All quantities are in Rydberg, before being converted to Hartree
    Ry2eV = 13.605704
   # Local energies at each orbital [Ry]
    Es  = 0.97100
    Ep  = 1.62396
    Ed1 = 0.59951
    Ed2 = 0.57185
    onsite = [Es,Ed1,Ed1,Ed1,Ed2,Ed2,Ep,Ep,Ep]

   # First neighbor [Ry]
    ss_sig = -0.08547
    pp_sig =  0.18728
    pp_pi  = -0.03000
    dd_sig = -0.08408
    dd_pi  =  0.04211
    dd_del = -0.00673
    sp_sig =  0.12129
    sd_sig = -0.07255
    pd_sig = -0.10110
    pd_pi  =  0.02813
    first_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Second neighbor [Ry]
    ss_sig = -0.00084
    pp_sig =  0.02919
    pp_pi  = -0.00647
    dd_sig = -0.00245
    dd_pi  =  0.00180
    dd_del =  0.00007
    sp_sig = -0.00120
    sd_sig = -0.00756
    pd_sig = -0.01703
    pd_pi  =  0.00559
    second_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Convert the SK parameters from Ry to Ha
    onsite = onsite./2
    first_neighbor = first_neighbor./2
    second_neighbor = second_neighbor./2



   # Ir optical properties (experimental)
   # -> taken from pp.2206-2207 (12-131--12-132) of Hayne's CRC Handbook of Chemistry and Physics, 97th ed. (2016)
   # Dielectric function calculated from the index of refraction n and the extinction coefficient k

   # Frequency Ir [eV]
    frequencies = [  0.10,      0.15,      0.20,     0.25,     0.30,     0.35,     0.40,     0.45,     0.50,     0.60,
                     0.70,      0.80,      0.90,     1.00,     1.10,     1.20,     1.30,     1.40,     1.50,     1.60,
                     1.70,      1.80,      1.90,     2.00,     2.10,     2.20,     2.30,     2.40,     2.50,     2.60,
                     2.70,      2.80,      2.90,     3.00,     3.20,     3.40,     3.60,     3.80,     4.00,     4.20,
                     4.40,      4.60,      4.80,     5.00,     5.20,     5.40,     5.60,     5.80,     6.00,     6.20,
                     6.40,      6.60,      6.80,     7.00,     7.20,     7.40,     7.60,     7.80,     8.00,     8.20,
                     8.40,      8.60,      8.80,     9.00,     9.20,     9.40,     9.60,     9.80,    10.00,    10.20,
                    10.40,     10.60,     10.80,    11.00,    11.20,    11.40,    11.60,    11.80,    12.00,    12.40,
                    12.80,     13.20,     13.60,    14.00,    14.40,    14.80,    15.20,    15.60,    16.00,    16.40,
                    16.80,     17.20,     17.60,    18.00,    18.40,    18.80,    19.20,    19.60,    20.00,    20.50,
                    21.00,     21.50,     22.00,    22.50,    23.00,    23.50,    24.00,    24.50,    25.00,    25.50,
                    26.00,     26.50,     27.00,    27.50,    28.00,    28.50,    29.00,    29.50,    30.00,    32.00,
                    34.00,     36.00,     38.00,    40.00]

   # Real part of dielectric constant of bulk Ir:  eps_real = n^2 - k^2
    eps_real = [-2863.104, -1803.820, -1155.020, -784.686, -561.437, -415.332, -314.467, -240.970, -188.803, -126.312,
                  -87.064,   -64.273,   -51.918,  -43.514,  -37.544,  -32.327,  -28.722,  -25.549,  -22.030,  -18.624,
                  -16.970,   -16.167,   -15.298,  -14.635,  -14.310,  -13.940,  -13.395,  -12.855,  -12.080,  -11.252,
                  -10.490,    -9.756,    -9.187,   -8.772,   -8.003,   -6.962,   -5.586,   -4.644,   -4.493,   -4.848,
                   -5.080,    -5.044,    -4.808,   -4.313,   -3.847,   -3.368,   -2.960,   -2.538,   -2.266,   -1.939,
                   -1.645,    -1.366,    -1.114,   -0.855,   -0.633,   -0.400,   -0.239,    0.043,    0.216,    0.392,
                    0.528,     0.666,     0.762,    0.885,    0.947,    1.030,    1.076,    1.094,    1.082,    1.021,
                    0.929,     0.857,     0.734,    0.628,    0.496,    0.439,    0.384,    0.353,    0.371,    0.362,
                    0.396,     0.432,     0.466,    0.518,    0.595,    0.636,    0.711,    0.751,    0.775,    0.841,
                    0.890,     0.950,     0.933,    0.825,    0.672,    0.538,    0.379,    0.220,    0.086,   -0.021,
                   -0.102,    -0.157,    -0.208,   -0.275,   -0.298,   -0.269,   -0.224,   -0.199,   -0.148,   -0.115,
                   -0.069,    -0.027,     0.000,    0.039,    0.064,    0.088,    0.098,    0.107,    0.129,    0.191,
                    0.287,     0.403,     0.475,    0.529]

   # Imaginary part of dielectric constant of bulk Ir:  eps_imag = 2nk
    eps_imag = [ 3454.128,  1383.396,   684.889,  395.685,  250.260,  170.894,  123.530,   96.502,   83.798,   64.616,
                   57.311,    54.071,    50.274,   46.053,   41.587,   37.947,   34.599,   31.226,   28.567,   27.229,
                   26.470,    25.397,    24.055,   22.850,   21.504,   20.060,   18.574,   17.140,   15.840,   14.745,
                   13.801,    13.068,    12.425,   11.868,   10.562,    9.333,    8.542,    8.662,    8.790,    8.564,
                    7.772,     6.812,     5.876,    5.170,    4.618,    4.180,    3.881,    3.571,    3.382,    3.158,
                    2.989,     2.820,     2.698,    2.600,    2.515,    2.448,    2.348,    2.290,    2.328,    2.360,
                    2.391,     2.419,     2.451,    2.500,    2.584,    2.641,    2.755,    2.851,    2.929,    3.016,
                    3.082,     3.117,     3.158,    3.119,    3.055,    2.961,    2.867,    2.750,    2.678,    2.541,
                    2.404,     2.313,     2.223,    2.111,    2.059,    2.053,    1.999,    1.992,    2.009,    2.017,
                    2.050,     2.125,     2.262,    2.418,    2.464,    2.480,    2.472,    2.415,    2.332,    2.184,
                    2.059,     1.918,     1.780,    1.663,    1.517,    1.398,    1.270,    1.162,    1.090,    1.034,
                    0.965,     0.925,     0.871,    0.832,    0.805,    0.779,    0.741,    0.704,    0.678,    0.546,
                    0.448,     0.373,     0.350,    0.334]

   # Interpolation of the optical frequencies
    eps_real_interp = LinearInterpolation(frequencies, eps_real)
    eps_imag_interp = LinearInterpolation(frequencies, eps_imag)
   #
    function dielectric(hw)
        return round(eps_real_interp(hw), digits=3) + 1im*round(eps_imag_interp(hw), digits=3)
    end


   # FCC lattice constant [nm] - 7.257 bohr
    a0 = 0.3840

   # KPM shift and scale
    min_bandval = -11.179502./Ry2eV./2  # min, no fermi [Ha]
    max_bandval =  19.096639./Ry2eV./2  # max, no fermi [Ha]
    fermi = 0.7620./2  # Fermi energy [Ha]

    A = (max_bandval-min_bandval)./2 + min_bandval + fermi
    B = (max_bandval-min_bandval)./2 * 1.05
    println("A = ", round(A, digits=3))
    println("B = ", round(B, digits=3))


    return [onsite, first_neighbor, second_neighbor, A, B, fermi, a0, dielectric]
end
