using Interpolations


function platinum()
   # Pt Slater-Koster 2-center parameters (approximation)
   # -> taken from pp.293 of Papaconstantopoulos' Handbook of the Band Structure of Elemental Solids, 2nd ed. (2015)
   # All quantities are in Rydberg, before being converted to Hartree
    Ry2eV = 13.605704
   # Local energies at each orbital [Ry]
    Es  = 0.77393
    Ep  = 1.48490
    Ed1 = 0.45468
    Ed2 = 0.43701
    onsite = [Es,Ed1,Ed1,Ed1,Ed2,Ed2,Ep,Ep,Ep]

   # First neighbor [Ry]
    ss_sig = -0.07835
    pp_sig =  0.18677
    pp_pi  = -0.02465
    dd_sig = -0.06856
    dd_pi  =  0.03528
    dd_del = -0.00588
    sp_sig =  0.11192
    sd_sig = -0.06197
    pd_sig = -0.08564
    pd_pi  =  0.02446
    first_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Second neighbor [Ry]
    ss_sig =  0.00157
    pp_sig =  0.03487
    pp_pi  = -0.00958
    dd_sig = -0.00300
    dd_pi  =  0.00226
    dd_del = -0.00029
    sp_sig =  0.00074
    sd_sig = -0.00770
    pd_sig = -0.01321
    pd_pi  =  0.00522
    second_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Convert the SK parameters from Ry to Ha
    onsite = onsite./2
    first_neighbor = first_neighbor./2
    second_neighbor = second_neighbor./2



   # Pt optical properties (experimental)
   # -> taken from pp.2212-2213 (12-137--12-138) of Hayne's CRC Handbook of Chemistry and Physics, 97th ed. (2016)
   # Dielectric function calculated from the index of refraction n and the extinction coefficient k

   # Frequency Pt [eV]
    frequencies = [  0.10,     0.15,     0.20,     0.25,     0.30,     0.35,     0.40,    0.45,    0.50,    0.55,
                     0.60,     0.65,     0.70,     0.75,     0.80,     0.85,     0.90,    0.95,    1.00,    1.10,
                     1.20,     1.30,     1.40,     1.50,     1.60,     1.70,     1.80,    1.90,    2.00,    2.10,
                     2.20,     2.30,     2.40,     2.50,     2.60,     2.70,     2.80,    2.90,    3.00,    3.20,
                     3.40,     3.60,     3.80,     4.00,     4.20,     4.40,     4.60,    4.80,    5.00,    5.20,
                     5.40,     5.60,     5.80,     6.00,     6.20,     6.40,     6.60,    6.80,    7.00,    7.20,
                     7.40,     7.60,     7.80,     8.00,     8.20,     8.40,     8.60,    8.80,    9.00,    9.20,
                     9.40,     9.60,     9.80,    10.00,    10.20,    10.40,    10.60,   10.80,   11.00,   11.20,
                    11.40,    11.60,    11.80,    12.00,    12.40,    12.80,    13.20,   13.60,   14.00,   14.40,
                    14.80,    15.20,    15.60,    16.00,    16.50,    17.00,    17.50,   18.00,   18.50,   19.00,
                    19.50,    20.00,    20.50,    21.00,    21.50,    22.00,    22.50,   23.00,   23.50,   24.00,
                    24.50,    25.00,    25.50,    26.00,    26.50,    27.00,    27.50,   28.00,   28.50,   29.00,
                    29.50,    30.00]

   # Real part of dielectric constant of bulk Pt:  eps_real = n^2 - k^2
    eps_real = [-1825.374, -904.033, -538.793, -354.270, -245.779, -175.837, -121.608, -77.495, -44.156, -30.003,
                  -19.246,  -13.885,  -14.045,  -18.256,  -21.366,  -23.218,  -24.995, -25.583, -25.762, -24.038,
                  -22.444,  -20.648,  -18.692,  -17.179,  -15.808,  -14.613,  -13.325, -12.483, -11.275, -10.394,
                   -9.504,   -9.059,   -8.411,   -7.855,   -7.242,   -6.743,   -6.261,  -5.856,  -5.464,  -4.795,
                   -4.208,   -3.654,   -3.276,   -2.842,   -2.477,   -2.117,   -1.870,  -1.518,  -1.248,  -0.939,
                   -0.743,   -0.522,   -0.311,   -0.056,    0.110,    0.352,    0.515,   0.653,   0.712,   0.688,
                    0.707,    0.702,    0.750,    0.769,    0.792,    0.838,    0.861,   0.884,   0.936,   0.988,
                    0.966,    0.943,    0.868,    0.809,    0.699,    0.638,    0.577,   0.568,   0.559,   0.550,
                    0.566,    0.583,    0.644,    0.664,    0.723,    0.781,    0.851,   0.851,   0.851,   0.825,
                    0.748,    0.748,    0.716,    0.746,    0.781,    0.823,    0.890,   0.942,   0.806,   0.658,
                    0.452,    0.269,    0.044,   -0.149,   -0.283,   -0.325,   -0.304,  -0.254,  -0.194,  -0.125,
                   -0.060,    0.000,    0.043,    0.099,    0.125,    0.151,    0.163,   0.203,   0.214,   0.226,
                    0.211,    0.197]

   # Imaginary part of dielectric constant of bulk Pt:  eps_imag = 2nk
    eps_imag = [ 1181.502,  509.778,  282.610,  182.360,  126.694,   89.610,   63.956,  56.419,  60.292,  65.402,
                   69.255,   73.526,   77.999,   78.203,   74.765,   70.498,   65.921,  60.930,  56.270,  48.173,
                   42.032,   36.914,   32.984,   29.609,   26.717,   24.406,   22.239,  20.278,  18.722,  17.483,
                   16.362,   15.414,   14.372,   13.406,   12.606,   11.968,   11.346,  10.776,  10.220,   9.274,
                    8.541,    7.837,    7.252,    6.705,    6.206,    5.834,    5.421,   5.106,   4.787,   4.542,
                    4.379,    4.189,    3.998,    3.864,    3.753,    3.664,    3.654,   3.670,   3.720,   3.750,
                    3.665,    3.611,    3.552,    3.469,    3.440,    3.381,    3.352,   3.322,   3.315,   3.308,
                    3.338,    3.367,    3.404,    3.358,    3.318,    3.220,    3.124,   3.024,   2.926,   2.830,
                    2.756,    2.683,    2.606,    2.580,    2.503,    2.425,    2.437,   2.437,   2.437,   2.418,
                    2.362,    2.362,    2.300,    2.207,    2.158,    2.150,    2.159,   2.306,   2.444,   2.534,
                    2.534,    2.502,    2.420,    2.266,    2.030,    1.810,    1.588,   1.417,   1.305,   1.214,
                    1.124,    1.066,    1.022,    0.992,    0.962,    0.932,    0.918,   0.900,   0.885,   0.870,
                    0.858,    0.847]

   # Interpolation of the optical frequencies
    eps_real_interp = LinearInterpolation(frequencies, eps_real)
    eps_imag_interp = LinearInterpolation(frequencies, eps_imag)
   #
    function dielectric(hw)
        return round(eps_real_interp(hw), digits=3) + 1im*round(eps_imag_interp(hw), digits=3)
    end


   # FCC lattice constant [nm] - 7.408 bohr
    a0 = 0.3920

   # KPM shift and scale
    min_bandval = -10.814462./Ry2eV./2  # min, no fermi [Ha]
    max_bandval =  19.431610./Ry2eV./2  # max, no fermi [Ha]
    fermi = 0.6380./2  # Fermi energy [Ha]

    A = (max_bandval-min_bandval)./2 + min_bandval + fermi
    B = (max_bandval-min_bandval)./2 * 1.05
    println("A = ", round(A, digits=3))
    println("B = ", round(B, digits=3))


    return [onsite, first_neighbor, second_neighbor, A, B, fermi, a0, dielectric]
end
