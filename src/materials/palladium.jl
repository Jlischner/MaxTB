using Interpolations


function palladium()
   # Pd Slater-Koster 2-center parameters (approximation)
   # -> taken from pp.237 of Papaconstantopoulos' Handbook of the Band Structure of Elemental Solids, 2nd ed. (2015)
   # All quantities are in Rydberg, before being converted to Hartree
    Ry2eV = 13.605704
   # Local energies at each orbital [Ry]
    Es  = 0.94261
    Ep  = 1.36110
    Ed1 = 0.37285
    Ed2 = 0.36265
    onsite = [Es,Ed1,Ed1,Ed1,Ed2,Ed2,Ep,Ep,Ep]

   # First neighbor [Ry]
    ss_sig = -0.07962
    pp_sig =  0.17119
    pp_pi  = -0.00540
    dd_sig = -0.05216
    dd_pi  =  0.02878
    dd_del = -0.00533
    sp_sig =  0.11332
    sd_sig = -0.04885
    pd_sig = -0.06563
    pd_pi  =  0.02124
    first_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Second neighbor [Ry]
    ss_sig = -0.00105
    pp_sig =  0.04282
    pp_pi  = -0.00044
    dd_sig = -0.00385
    dd_pi  =  0.00212
    dd_del = -0.00026
    sp_sig =  0.01048
    sd_sig = -0.00837
    pd_sig = -0.00738
    pd_pi  =  0.00351
    second_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Convert the SK parameters from Ry to Ha
    onsite = onsite./2
    first_neighbor = first_neighbor./2
    second_neighbor = second_neighbor./2



   # Pd optical properties (experimental)
   # -> taken from pp.2212 (12-137) of Hayne's CRC Handbook of Chemistry and Physics, 97th ed. (2016)
   # Dielectric function calculated from the index of refraction n and the extinction coefficient k

   # Frequency Pd [eV]
    frequencies = [  0.10,      0.15,     0.20,     0.26,     0.30,     0.36,     0.40,     0.46,     0.50,    0.56,
                     0.60,      0.72,     0.80,     1.00,     1.10,     1.20,     1.30,     1.40,     1.50,    1.60,
                     1.70,      1.80,     1.90,     2.00,     2.10,     2.20,     2.30,     2.40,     2.50,    2.60,
                     2.70,      2.80,     2.90,     3.00,     3.10,     3.20,     3.30,     3.40,     3.50,    3.60,
                     3.70,      3.80,     3.90,     4.00,     4.20,     4.40,     4.60,     4.80,     5.00,    5.20,
                     5.40,      5.60,     5.80,     6.00,     6.20,     6.40,     6.60,     6.80,     7.00,    7.20,
                     7.40,      7.60,     7.80,     8.00,     8.20,     8.40,     8.60,     8.80,     9.00,    9.50,
                    10.00,     10.50,    11.00,    11.50,    12.00,    12.50,    13.00,    13.50,    14.00,   14.50,
                    15.00,     15.50,    16.00,    16.50,    17.00,    17.50,    18.00,    18.50,    19.00,   19.50,
                    20.00,     20.50,    21.00,    21.50,    22.00,    22.50,    23.00,    23.50,    24.00,   25.00,
                    26.40,     27.80,    29.20]

   # Real part of dielectric constant of bulk Pd:  eps_real = n^2 - k^2
    eps_real = [-2915.166, -1273.276, -697.603, -396.350, -285.579, -191.808, -157.860, -128.419, -114.064, -94.674,
                  -84.762,   -63.370,  -53.741,  -38.532,  -33.836,  -30.188,  -27.158,  -24.774,  -22.540, -20.176,
                  -18.278,   -16.925,  -15.610,  -14.410,  -13.452,  -12.494,  -11.722,  -10.871,  -10.122,  -9.413,
                   -8.820,    -8.133,   -7.593,   -7.131,   -6.683,   -6.304,   -5.883,   -5.506,   -5.140,  -4.836,
                   -4.520,    -4.213,   -3.960,   -3.735,   -3.287,   -2.979,   -2.703,   -2.590,   -2.538,  -2.394,
                   -2.168,    -1.968,   -1.763,   -1.525,   -1.329,   -1.131,   -0.946,   -0.744,   -0.570,  -0.398,
                   -0.235,    -0.098,    0.080,    0.242,    0.394,    0.432,    0.578,    0.659,    0.735,   0.832,
                    0.877,     0.923,    0.983,    0.994,    1.004,    0.967,    0.944,    0.944,    0.920,   0.860,
                    0.801,     0.748,    0.731,    0.727,    0.773,    0.752,    0.797,    0.797,    0.818,   0.794,
                    0.722,     0.612,    0.531,    0.467,    0.419,    0.390,    0.392,    0.409,    0.414,   0.396,
                    0.455,     0.512,    0.550]

   # Imaginary part of dielectric constant of bulk Pd:  eps_imag = 2nk
    eps_imag = [  447.279,   224.233,  163.263,  125.333,  122.962,  114.704,  113.326,  103.419,   93.808,  82.242,
                   75.696,    61.074,   54.002,   41.202,   36.305,   32.330,   28.900,   25.740,   22.655,  20.592,
                   18.880,    17.434,   15.834,   14.630,   13.460,   12.416,   11.475,   10.613,    9.814,   9.206,
                    8.580,     8.075,    7.636,    7.232,    6.840,    6.482,    6.110,    5.824,    5.544,   5.292,
                    5.093,     4.897,    4.725,    4.511,    4.347,    4.141,    3.996,    3.838,    3.571,   3.222,
                    2.890,     2.624,    2.402,    2.204,    2.028,    1.883,    1.742,    1.650,    1.533,   1.470,
                    1.401,     1.343,    1.295,    1.285,    1.316,    1.344,    1.300,    1.352,    1.370,   1.456,
                    1.482,     1.508,    1.510,    1.547,    1.584,    1.595,    1.581,    1.581,    1.568,   1.564,
                    1.559,     1.496,    1.426,    1.336,    1.305,    1.293,    1.263,    1.263,    1.274,   1.318,
                    1.391,     1.380,    1.327,    1.254,    1.165,    1.091,    1.015,    0.952,    0.907,   0.826,
                    0.688,     0.616,    0.574]

   # Interpolation of the optical frequencies
    eps_real_interp = LinearInterpolation(frequencies, eps_real)
    eps_imag_interp = LinearInterpolation(frequencies, eps_imag)
   #
    function dielectric(hw)
        return round(eps_real_interp(hw), digits=3) + 1im*round(eps_imag_interp(hw), digits=3)
    end


   # FCC lattice constant [nm] - 7.351 bohr
    a0 = 0.3890

   # KPM shift and scale
    min_bandval = -7.321616./Ry2eV./2  # min, no fermi [Ha]
    max_bandval = 21.327423./Ry2eV./2  # max, no fermi [Ha]
    fermi = 0.5190./2  # Fermi energy [Ha]

    A = (max_bandval-min_bandval)./2 + min_bandval + fermi
    B = (max_bandval-min_bandval)./2 * 1.05
    println("A = ", round(A, digits=3))
    println("B = ", round(B, digits=3))


    return [onsite, first_neighbor, second_neighbor, A, B, fermi, a0, dielectric]
end
