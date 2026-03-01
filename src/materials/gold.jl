using Interpolations


function gold()
   # Au Slater-Koster 2-center parameters (approximation)
   # -> taken from pp.298 of Papaconstantopoulos' Handbook of the Band Structure of Elemental Solids, 2nd ed. (2015)
   # All quantities are in Rydberg, before being converted to Hartree
    Ry2eV = 13.605704
   # Local energies at each orbital [Ry]
    Es  = 0.56220
    Ep  = 1.27897
    Ed1 = 0.26097
    Ed2 = 0.25309
    onsite = [Es,Ed1,Ed1,Ed1,Ed2,Ed2,Ep,Ep,Ep]

   # First neighbor [Ry]
    ss_sig = -0.06680
    pp_sig =  0.17866
    pp_pi  = -0.01645
    dd_sig = -0.04971
    dd_pi  =  0.02624
    dd_del = -0.00457
    sp_sig =  0.09721
    sd_sig = -0.04722
    pd_sig = -0.06399
    pd_pi  =  0.01896
    first_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Second neighbor [Ry]
    ss_sig =  0.00277
    pp_sig =  0.03707
    pp_pi  = -0.01025
    dd_sig = -0.00305
    dd_pi  =  0.00240
    dd_del = -0.00057
    sp_sig =  0.00261
    sd_sig = -0.00784
    pd_sig = -0.00762
    pd_pi  =  0.00470
    second_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Convert the SK parameters from Ry to Ha
    onsite = onsite./2
    first_neighbor = first_neighbor./2
    second_neighbor = second_neighbor./2



   # Electropolished Au(110) optical properties (experimental)
   # -> taken from pp.2204-2205 (12-129--12-130) of Hayne's CRC Handbook of Chemistry and Physics, 97th ed. (2016)
   # Dielectric function calculated from the index of refraction n and the extinction coefficient k

   # Frequency Au(110) [eV]
    frequencies = [  0.10,      0.20,     0.30,     0.40,     0.50,     0.60,     0.70,     0.80,    0.90,    1.00,
                     1.20,      1.40,     1.60,     1.80,     2.00,     2.10,     2.20,     2.40,    2.50,    2.60,
                     2.70,      2.80,     2.90,     3.00,     3.10,     3.20,     3.30,     3.40,    3.50,    3.60,
                     3.70,      3.80,     3.90,     4.00,     4.10,     4.20,     4.30,     4.40,    4.50,    4.60,
                     4.70,      4.80,     4.90,     5.00,     5.20,     5.40,     5.60,     5.80,    6.00,    6.20,
                     6.40,      6.60,     6.80,     7.00,     7.20,     7.40,     7.60,     7.80,    8.00,    8.20,
                     8.40,      8.60,     8.80,     9.00,     9.20,     9.40,     9.60,     9.80,   10.00,   10.20,
                    10.40,     10.60,    10.80,    11.00,    11.20,    11.40,    11.60,    11.80,   12.00,   12.40,
                    12.80,     13.20,    13.60,    14.00,    14.40,    14.80,    15.20,    15.60,   16.00,   16.40,
                    16.80,     17.20,    17.60,    18.00,    18.40,    18.80,    19.20,    19.60,   20.00,   20.40,
                    20.80,     21.20,    21.60,    22.00,    22.40,    22.80,    23.20,    23.60,   24.00,   24.40,
                    24.80,     25.20,    25.60,    26.00,    26.40,    26.80,    27.20,    27.60,   28.00,   28.40,
                    28.80,     29.20,    29.60,    30.00]

   # Real part of dielectric constant of Au(110):  eps_real = n^2 - k^2
    eps_real = [-6794.060, -1736.856, -772.972, -433.541, -275.740, -189.810, -138.014, -104.212, -81.158, -64.464,
                  -42.762,   -29.587,  -20.787,  -14.584,   -9.969,   -8.033,   -6.394,   -3.210,  -1.856,  -0.834,
                   -0.914,    -1.001,   -0.954,   -0.868,   -0.905,   -0.868,   -0.766,   -0.664,  -0.497,  -0.369,
                   -0.373,    -0.547,   -0.748,   -0.874,   -0.924,   -0.978,   -1.030,   -1.110,  -1.205,  -1.166,
                   -1.077,    -0.966,   -0.859,   -0.732,   -0.496,   -0.305,   -0.149,    0.024,   0.189,   0.350,
                    0.460,     0.593,    0.749,    0.893,    0.947,    0.983,    0.944,    0.842,   0.795,   0.844,
                    0.898,     0.942,    0.977,    1.001,    1.060,    1.161,    1.241,    1.253,   1.237,   1.210,
                    1.183,     1.172,    1.203,    1.218,    1.248,    1.290,    1.331,    1.400,   1.428,   1.541,
                    1.478,     1.311,    1.137,    1.029,    0.925,    0.882,    0.849,    0.832,   0.840,   0.832,
                    0.839,     0.854,    0.869,    0.869,    0.883,    0.892,    0.887,    0.824,   0.704,   0.577,
                    0.453,     0.328,    0.226,    0.144,    0.103,    0.098,    0.110,    0.150,   0.204,   0.256,
                    0.304,     0.359,    0.397,    0.435,    0.462,    0.473,    0.500,    0.500,   0.527,   0.544,
                    0.544,     0.544,    0.527,    0.509]

   # Imaginary part of dielectric constant of Au(110):  eps_imag = 2nk
    eps_imag = [ 1353.442,   177.770,   55.084,   24.579,   12.956,    7.717,    5.170,    3.676,   2.703,   2.088,
                    1.308,     0.870,    0.730,    0.688,    0.822,    1.022,    1.219,    1.860,   2.608,   3.819,
                    4.919,     5.168,    5.370,    5.544,    5.575,    5.544,    5.518,    5.491,   5.467,   5.605,
                    5.740,     5.835,    5.756,    5.611,    5.406,    5.269,    5.133,    4.963,   4.698,   4.394,
                    4.166,     3.975,    3.788,    3.636,    3.388,    3.219,    3.073,    2.904,   2.782,   2.703,
                    2.625,     2.565,    2.522,    2.546,    2.584,    2.650,    2.705,    2.673,   2.515,   2.392,
                    2.314,     2.306,    2.253,    2.158,    2.122,    2.075,    2.122,    2.165,   2.192,   2.176,
                    2.160,     2.117,    2.064,    2.037,    1.983,    1.971,    1.958,    1.960,   1.974,   2.102,
                    2.291,     2.386,    2.356,    2.288,    2.219,    2.117,    2.058,    1.976,   1.912,   1.872,
                    1.809,     1.785,    1.761,    1.761,    1.737,    1.776,    1.839,    1.936,   1.959,   1.938,
                    1.914,     1.848,    1.760,    1.617,    1.477,    1.343,    1.230,    1.120,   1.056,   0.992,
                    0.928,     0.918,    0.896,    0.874,    0.867,    0.850,    0.843,    0.843,   0.835,   0.845,
                    0.845,     0.845,    0.835,    0.826]

   # Interpolation of the optical frequencies
    eps_real_interp = LinearInterpolation(frequencies, eps_real)
    eps_imag_interp = LinearInterpolation(frequencies, eps_imag)
   #
    function dielectric(hw)
        return round(eps_real_interp(hw), digits=3) + 1im*round(eps_imag_interp(hw), digits=3)
    end


   # FCC lattice constant [nm] - 7.710 bohr
    a0 = 0.4080

   # KPM shift and scale
    min_bandval = -10.350917./Ry2eV./2  # min, no fermi [Ha]
    max_bandval =  18.464928./Ry2eV./2  # max, no fermi [Ha]
    fermi = 0.5380./2  # Fermi energy [Ha]

    A = (max_bandval-min_bandval)./2 + min_bandval + fermi
    B = (max_bandval-min_bandval)./2 * 1.05
    println("A = ", round(A, digits=3))
    println("B = ", round(B, digits=3))


    return [onsite, first_neighbor, second_neighbor, A, B, fermi, a0, dielectric]
end
