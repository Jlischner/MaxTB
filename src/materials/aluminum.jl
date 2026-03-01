using Interpolations


function aluminum()
   # Al Slater-Koster 2-center parameters (approximation)
   # -> taken from pp.312 of Papaconstantopoulos' Handbook of the Band Structure of Elemental Solids, 2nd ed. (2015)
   # All quantities are in Rydberg, before being converted to Hartree
    Ry2eV = 13.605704
   # Local energies at each orbital [Ry]
    Es  = 0.48220
    Ep  = 1.04964
    Ed1 = 1.73912
    Ed2 = 1.62982
    onsite = [Es,Ed1,Ed1,Ed1,Ed2,Ed2,Ep,Ep,Ep]

   # First neighbor [Ry]
    ss_sig = -0.05831
    pp_sig =  0.17149
    pp_pi  = -0.01020
    dd_sig = -0.17004
    dd_pi  =  0.07394
    dd_del = -0.00490
    sp_sig = -0.09315
    sd_sig =  0.08296
    pd_sig = -0.16259
    pd_pi  =  0.03088
    first_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Second neighbor [Ry]
    ss_sig = -0.00163
    pp_sig =  0.00026
    pp_pi  = -0.00316
    dd_sig = -0.01163
    dd_pi  =  0.00315
    dd_del = -0.00283
    sp_sig =  0.00519
    sd_sig = -0.03032
    pd_sig =  0.02447
    pd_pi  = -0.00266
    second_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Convert the SK parameters from Ry to Ha
    onsite = onsite./2
    first_neighbor = first_neighbor./2
    second_neighbor = second_neighbor./2



   # Al optical properties (experimental)
   # -> taken from pp.2200-2201 (12-125--12-126) of Hayne's CRC Handbook of Chemistry and Physics, 97th ed. (2016)
   # Dielectric function calculated from the index of refraction n and the extinction coefficient k

   # Frequency Al [eV]
    frequencies = [   0.04,       0.05,       0.06,       0.07,       0.08,       0.09,      0.10,      0.13,      0.15,      0.18,
                      0.20,       0.25,       0.30,       0.35,       0.40,       0.50,      0.60,      0.70,      0.80,      0.90,
                      1.00,       1.10,       1.20,       1.30,       1.40,       1.50,      1.60,      1.70,      1.80,      1.90,
                      2.00,       2.20,       2.40,       2.60,       2.80,       3.00,      3.20,      3.40,      3.60,      3.80,
                      4.00,       4.20,       4.40,       4.60,       4.80,       5.00,      6.00,      6.50,      7.00,      7.50,
                      8.00,       8.50,       9.00,       9.50,      10.00,      10.50,     11.00,     11.50,     12.00,     12.50,
                     13.00,      13.50,      14.00,      14.20,      14.40,      14.60,     14.80,     15.00,     15.20,     15.40,
                     15.60,      15.80,      16.00,      16.20,      16.40,      16.75,     17.00,     17.25,     17.50,     17.75,
                     18.00,      18.50,      19.00,      19.50,      20.00,      20.50,     21.00,     21.50,     22.00,     22.50,
                     23.00,      23.50,      24.00,      24.50,      25.00,      25.50,     26.00,     27.00,     28.00,     29.00,
                     30.00,      35.00,      40.00,      45.00,      50.00,      55.00,     60.00,     65.00,     70.00,     72.50,
                     75.00,      77.50,      80.00,      85.00,      90.00,      95.00,    100.00,    110.00,    120.00,    130.00,
                    140.00,     150.00,     160.00,     170.00,     180.00,     190.00,    200.00,    220.00,    240.00,    260.00,
                    280.00,     300.00]

   # Real part of dielectric constant of bulk Al:  eps_real = n^2 - k^2
    eps_real = [-31773.123, -24027.946, -18789.964, -15466.886, -13213.928, -11447.065, -9963.593, -7342.311, -5577.922, -4275.878,
                 -3387.134,  -2252.896,  -1632.038,  -1237.717,   -971.467,   -644.950,  -452.922,  -332.783,  -252.477,  -194.991,
                  -153.882,   -123.572,    -98.613,    -77.930,    -62.433,    -61.504,   -67.018,   -68.904,   -64.291,   -58.954,
                   -54.235,    -45.831,    -38.794,    -33.157,    -28.641,    -24.967,   -21.954,   -19.424,   -17.291,   -15.465,
                   -13.901,    -12.545,    -11.365,    -10.332,     -9.420,     -8.619,    -5.700,    -4.710,    -3.923,    -3.284,
                    -2.760,     -2.328,     -1.962,     -1.651,     -1.386,     -1.156,    -0.957,    -0.779,    -0.625,    -0.489,
                    -0.369,     -0.266,     -0.172,     -0.136,     -0.104,     -0.070,    -0.037,    -0.008,     0.020,     0.021,
                     0.073,      0.097,      0.120,      0.141,      0.163,      0.199,     0.223,     0.246,     0.269,     0.290,
                     0.310,      0.348,      0.384,      0.417,      0.445,      0.474,     0.499,     0.524,     0.546,     0.567,
                     0.586,      0.605,      0.622,      0.638,      0.654,      0.667,     0.682,     0.705,     0.729,     0.748,
                     0.767,      0.837,      0.884,      0.916,      0.939,      0.958,     0.974,     0.990,     1.012,     1.051,
                     1.022,      1.015,      1.013,      1.013,      1.009,      0.997,     0.981,     0.987,     0.982,     0.974,
                     0.978,      0.980,      0.978,      0.978,      0.980,      0.980,     0.982,     0.984,     0.986,     0.986,
                     0.988,      0.990]

   # Imaginary part of dielectric constant of bulk Al:  eps_imag = 2nk
    eps_imag = [ 40167.800,  25828.817,  18956.037,  14577.090,  11330.075,   9048.517,  7278.797,  4456.253,  2858.602,  1910.718,
                  1393.176,    828.291,    553.697,    387.175,    280.468,    157.170,    97.298,    64.881,    46.078,    35.445,
                    30.213,     26.857,     25.225,     26.274,     36.740,     45.616,    45.134,    36.744,    28.570,    23.275,
                    19.505,     13.938,     10.380,      8.062,      6.440,      5.255,     4.331,     3.603,     3.030,     2.573,
                     2.199,      1.897,      1.649,      1.437,      1.261,      1.118,     0.622,     0.478,     0.377,     0.297,
                     0.239,      0.192,      0.157,      0.126,      0.104,      0.086,     0.070,     0.058,     0.052,     0.048,
                     0.046,      0.042,      0.040,      0.040,      0.038,      0.037,     0.036,     0.038,     0.038,     0.086,
                     0.041,      0.041,      0.042,      0.042,      0.041,      0.040,     0.040,     0.040,     0.040,     0.039,
                     0.039,      0.038,      0.037,      0.036,      0.036,      0.034,     0.034,     0.033,     0.033,     0.032,
                     0.032,      0.031,      0.030,      0.029,      0.029,      0.028,     0.026,     0.025,     0.024,     0.024,
                     0.023,      0.018,      0.015,      0.013,      0.012,      0.010,     0.008,     0.008,     0.008,     0.008,
                     0.049,      0.050,      0.048,      0.056,      0.062,      0.072,     0.059,     0.050,     0.048,     0.041,
                     0.032,      0.030,      0.028,      0.022,      0.020,      0.018,     0.014,     0.012,     0.010,     0.008,
                     0.006,      0.004]

   # Interpolation of the optical frequencies
    eps_real_interp = LinearInterpolation(frequencies, eps_real)
    eps_imag_interp = LinearInterpolation(frequencies, eps_imag)
   #
    function dielectric(hw)
        return round(eps_real_interp(hw), digits=3) + 1im*round(eps_imag_interp(hw), digits=3)
    end


   # FCC lattice constant [nm] - 7.653 bohr
    a0 = 0.4050

   # KPM shift and scale
    min_bandval = -11.521277./Ry2eV./2  # min, no fermi [Ha]
    max_bandval =  22.488803./Ry2eV./2  # max, no fermi [Ha]
    fermi = 0.6195./2  # Fermi energy [Ha]

    A = (max_bandval-min_bandval)./2 + min_bandval + fermi
    B = (max_bandval-min_bandval)./2 * 1.05
    println("A = ", round(A, digits=3))
    println("B = ", round(B, digits=3))


    return [onsite, first_neighbor, second_neighbor, A, B, fermi, a0, dielectric]
end
