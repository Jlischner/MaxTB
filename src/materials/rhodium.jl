using Interpolations


function rhodium()
   # Rh Slater-Koster 2-center parameters (approximation)
   # -> taken from pp.232 of Papaconstantopoulos' Handbook of the Band Structure of Elemental Solids, 2nd ed. (2015)
   # All quantities are in Rydberg, before being converted to Hartree
    Ry2eV = 13.605704
   # Local energies at each orbital [Ry]
    Es  = 1.12862
    Ep  = 1.50833
    Ed1 = 0.49824
    Ed2 = 0.48265
    onsite = [Es,Ed1,Ed1,Ed1,Ed2,Ed2,Ep,Ep,Ep]

   # First neighbor [Ry]
    ss_sig = -0.08700
    pp_sig =  0.17058
    pp_pi  = -0.00928
    dd_sig = -0.06303
    dd_pi  =  0.03409
    dd_del = -0.00610
    sp_sig =  0.11907
    sd_sig = -0.05607
    pd_sig = -0.07627
    pd_pi  =  0.02397
    first_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Second neighbor [Ry]
    ss_sig = -0.00094
    pp_sig =  0.04716
    pp_pi  = -0.00153
    dd_sig = -0.00405
    dd_pi  =  0.00206
    dd_del = -0.00013
    sp_sig =  0.00790
    sd_sig = -0.01022
    pd_sig = -0.01194
    pd_pi  =  0.00529
    second_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Convert the SK parameters from Ry to Ha
    onsite = onsite./2
    first_neighbor = first_neighbor./2
    second_neighbor = second_neighbor./2



   # Rh optical properties (experimental)
   # -> taken from pp.2215-2216 (12-140--12-141) of Hayne's CRC Handbook of Chemistry and Physics, 97th ed. (2016)
   # Dielectric function calculated from the index of refraction n and the extinction coefficient k

   # Frequency Rh [eV]
    frequencies = [  0.10,      0.20,     0.30,     0.40,     0.50,     0.60,     0.70,    0.80,    0.90,    1.00,
                     1.10,      1.20,     1.30,     1.40,     1.50,     1.60,     1.70,    1.80,    1.90,    2.00,
                     2.10,      2.20,     2.30,     2.40,     2.50,     2.60,     2.70,    2.90,    3.00,    3.10,
                     3.20,      3.30,     3.40,     3.50,     3.60,     3.70,     3.80,    3.90,    4.00,    4.20,
                     4.40,      4.60,     4.80,     5.00,     5.20,     5.40,     5.60,    5.80,    6.00,    6.20,
                     6.40,      6.60,     6.80,     7.00,     7.20,     7.40,     7.60,    7.80,    8.00,    8.20,
                     8.40,      8.60,     8.80,     9.00,     9.20,     9.40,     9.60,    9.80,   10.00,   10.60,
                    11.00,     11.60,    12.00,    12.60,    13.00,    13.60,    14.00,   14.60,   15.00,   15.60,
                    16.00,     16.50,    17.00,    17.50,    18.00,    18.50,    19.00,   19.50,   20.00,   20.50,
                    21.00,     21.50,    22.00,    22.50,    23.00,    23.50,    24.00,   24.50,   25.00,   25.50,
                    26.00,     26.50,    27.00,    27.50,    28.00,    29.00,    30.00,   31.00,   32.00,   33.00,
                    34.00,     35.00,    36.00,    37.00,    38.00,    39.00]

   # Real part of dielectric constant of bulk Rh:  eps_real = n^2 - k^2
    eps_real = [-4479.015, -1328.256, -638.661, -369.572, -240.605, -167.543, -123.890, -93.739, -74.505, -61.405,
                  -54.759,   -50.724,  -47.589,  -44.376,  -40.853,  -37.330,  -34.213, -30.950, -28.338, -25.866,
                  -23.888,   -22.112,  -20.640,  -19.238,  -18.088,  -17.280,  -16.920, -16.353, -16.063, -15.652,
                  -15.038,   -14.321,  -13.514,  -12.683,  -11.836,  -11.000,  -10.328,  -9.659,  -8.995,  -7.955,
                   -6.978,    -6.152,   -5.428,   -4.852,   -4.349,   -3.940,   -3.604,  -3.376,  -3.147,  -2.890,
                   -2.643,    -2.394,   -2.111,   -1.875,   -1.609,   -1.387,   -1.164,  -0.978,  -0.778,  -0.578,
                   -0.393,    -0.184,    0.000,    0.182,    0.370,    0.516,    0.669,   0.778,   0.893,   1.055,
                    1.087,     1.102,    1.070,    1.070,    1.054,    1.020,    1.003,   0.898,   0.828,   0.753,
                    0.746,     0.739,    0.714,    0.732,    0.739,    0.716,    0.577,   0.290,   0.022,  -0.188,
                   -0.274,    -0.260,   -0.214,   -0.190,   -0.186,   -0.194,   -0.173,  -0.166,  -0.117,  -0.099,
                   -0.054,    -0.013,    0.013,    0.050,    0.074,    0.131,    0.176,   0.170,   0.179,   0.223,
                    0.333,     0.398,    0.460,    0.469,    0.475,    0.500]

   # Imaginary part of dielectric constant of bulk Rh:  eps_imag = 2nk
    eps_imag = [ 2566.133,   648.807,  303.498,  187.704,  134.988,  104.567,   86.025,  75.068,  67.766,  64.331,
                   60.628,    55.739,   49.748,   44.006,   38.753,   34.528,   30.637,  27.692,  25.344,  23.362,
                   21.730,    20.440,   19.167,   18.164,   17.484,   16.835,   16.164,  14.214,  13.127,  11.844,
                   10.634,     9.528,    8.525,    7.717,    7.088,    6.555,    6.079,   5.685,   5.366,   4.880,
                    4.416,     4.056,    3.887,    3.697,    3.523,    3.424,    3.296,   3.160,   2.934,   2.701,
                    2.478,     2.298,    2.144,    2.006,    1.888,    1.782,    1.702,   1.632,   1.546,   1.477,
                    1.436,     1.388,    1.378,    1.355,    1.387,    1.434,    1.477,   1.546,   1.615,   1.840,
                    1.961,     2.112,    2.165,    2.165,    2.191,    2.244,    2.270,   2.314,   2.304,   2.250,
                    2.207,     2.165,    2.147,    2.123,    2.165,    2.300,    2.430,   2.478,   2.398,   2.180,
                    1.911,     1.720,    1.577,    1.490,    1.422,    1.305,    1.226,   1.134,   1.063,   0.992,
                    0.924,     0.858,    0.832,    0.793,    0.767,    0.702,    0.673,   0.627,   0.537,   0.444,
                    0.390,     0.386,    0.394,    0.414,    0.400,    0.375]

   # Interpolation of the optical frequencies
    eps_real_interp = LinearInterpolation(frequencies, eps_real)
    eps_imag_interp = LinearInterpolation(frequencies, eps_imag)
   #
    function dielectric(hw)
        return round(eps_real_interp(hw), digits=3) + 1im*round(eps_imag_interp(hw), digits=3)
    end


   # FCC lattice constant [nm] - 7.181 bohr
    a0 = 0.3800

   # KPM shift and scale
    min_bandval = -7.544613./Ry2eV./2  # min, no fermi [Ha]
    max_bandval = 21.375996./Ry2eV./2  # max, no fermi [Ha]
    fermi = 0.6335./2  # Fermi energy [Ha]

    A = (max_bandval-min_bandval)./2 + min_bandval + fermi
    B = (max_bandval-min_bandval)./2 * 1.05
    println("A = ", round(A, digits=3))
    println("B = ", round(B, digits=3))


    return [onsite, first_neighbor, second_neighbor, A, B, fermi, a0, dielectric]
end
