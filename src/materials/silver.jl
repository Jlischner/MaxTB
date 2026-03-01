using Interpolations


function silver()
   # Ag Slater-Koster 2-center parameters (approximation)
   # -> taken from pp.242 of Papaconstantopoulos' Handbook of the Band Structure of Elemental Solids, 2nd ed. (2015)
   # All quantities are in Rydberg, before being converted to Hartree
    Ry2eV = 13.605704
   # Local energies at each orbital [Ry]
    Es  = 0.68297
    Ep  = 1.13432
    Ed1 = 0.12249
    Ed2 = 0.12006
    onsite = [Es,Ed1,Ed1,Ed1,Ed2,Ed2,Ep,Ep,Ep]

   # First neighbor [Ry]
    ss_sig = -0.06581
    pp_sig =  0.15752
    pp_pi  =  0.00649
    dd_sig = -0.03151
    dd_pi  =  0.01757
    dd_del = -0.00336
    sp_sig =  0.09781
    sd_sig = -0.03110
    pd_sig = -0.03905
    pd_pi  =  0.01519
    first_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Second neighbor [Ry]
    ss_sig =  0.00143
    pp_sig =  0.03971
    pp_pi  =  0.00434
    dd_sig = -0.00282
    dd_pi  =  0.00171
    dd_del = -0.00038
    sp_sig =  0.00545
    sd_sig = -0.00462
    pd_sig = -0.00065
    pd_pi  =  0.00172
    second_neighbor = [sp_sig, ss_sig, pp_sig, pp_pi, sd_sig, pd_sig, pd_pi, dd_sig, dd_pi, dd_del]

   # Convert the SK parameters from Ry to Ha
    onsite = onsite./2
    first_neighbor = first_neighbor./2
    second_neighbor = second_neighbor./2



   # Ag optical properties (experimental)
   # -> taken from pp.2218-2219 (12-143--12-144) of Hayne's CRC Handbook of Chemistry and Physics, 97th ed. (2016)
   # Dielectric function calculated from the index of refraction n and the extinction coefficient k

   # Frequency Ag [eV]
    frequencies = [  0.10,      0.20,     0.30,     0.40,     0.50,    1.00,    1.50,    2.00,   2.50,   3.00,
                     3.25,      3.50,     3.60,     3.70,     3.77,    3.80,    3.90,    4.00,   4.10,   4.20,
                     4.30,      4.50,     4.75,     5.00,     5.50,    6.00,    6.50,    7.00,   7.50,   8.00,
                     9.00,     10.00,    11.00,    12.00,    13.00,   14.00,   14.50,   15.00,  16.00,  17.00,
                    18.00,     19.00,    20.00,    21.00,    21.50,   22.00,   22.50,   23.00,  23.50,  24.00,
                    24.50,     25.00,    25.50,    26.00,    26.50,   27.00,   27.50,   28.00,  28.50,  29.00,
                    30.00,     31.00,    32.00,    33.00,    34.00,   35.00,   36.00,   38.00,  40.00,  42.00,
                    44.00,     46.00,    48.00,    50.00,    52.00,   54.00,   56.00,   58.00,  60.00,  62.00,
                    64.00,     66.00,    68.00,    70.00,    72.00,   74.00,   76.00,   78.00,  80.00,  85.00,
                    90.00,     95.00,   100.00]

   # Real part of dielectric constant of bulk Ag:  eps_real = n^2 - k^2
    eps_real = [-8050.465, -2080.424, -928.872, -523.124, -335.174, -81.463, -33.451, -17.400, -9.491, -5.100,
                   -3.407,    -1.972,   -1.224,   -0.503,    0.121,   0.443,   1.560,   2.232,  2.270,  1.939,
                    1.716,     1.218,    0.797,    0.553,    0.307,   0.157,   0.170,   0.269,  0.472,  0.783,
                    1.455,     1.818,    1.997,    2.244,    2.346,   2.350,   1.915,   1.587,  1.188,  1.029,
                    0.998,     1.050,    1.160,    1.260,    1.237,   1.039,   0.723,   0.485,  0.345,  0.272,
                    0.223,     0.214,    0.220,    0.262,    0.316,   0.370,   0.408,   0.462,  0.503,  0.533,
                    0.573,     0.584,    0.566,    0.550,    0.534,   0.537,   0.599,   0.640,  0.673,  0.688,
                    0.701,     0.708,    0.696,    0.690,    0.714,   0.746,   0.689,   0.699,  0.709,  0.730,
                    0.730,     0.730,    0.713,    0.649,    0.690,   0.694,   0.697,   0.700,  0.703,  0.710,
                    0.716,     0.736,    0.755]

   # Imaginary part of dielectric constant of bulk Ag:  eps_imag = 2nk
    eps_imag = [ 1789.151,   259.576,   86.038,   41.660,   24.549,   5.057,   3.127,   2.257,  1.483,  1.044,
                    0.856,     0.596,    0.520,    0.462,    0.424,   0.438,   0.936,   1.932,  2.941,  3.710,
                    3.910,     4.326,    4.315,    4.216,    3.886,   3.430,   2.950,   2.502,  2.075,  1.740,
                    1.490,     1.635,    1.702,    1.900,    2.125,   2.683,   2.886,   2.870,  2.584,  2.288,
                    2.048,     1.905,    1.832,    2.025,    2.192,   2.332,   2.344,   2.200,  2.046,  1.872,
                    1.723,     1.577,    1.420,    1.332,    1.228,   1.157,   1.104,   1.062,  1.037,  1.030,
                    1.004,     0.986,    0.975,    0.918,    0.862,   0.774,   0.783,   0.694,  0.666,  0.630,
                    0.594,     0.576,    0.552,    0.510,    0.498,   0.299,   0.452,   0.418,  0.383,  0.370,
                    0.370,     0.370,    0.365,    0.332,    0.306,   0.289,   0.272,   0.255,  0.238,  0.187,
                    0.136,     0.103,    0.070]

   # Interpolation of the optical frequencies
    eps_real_interp = LinearInterpolation(frequencies, eps_real)
    eps_imag_interp = LinearInterpolation(frequencies, eps_imag)
   #
    function dielectric(hw)
        return round(eps_real_interp(hw), digits=3) + 1im*round(eps_imag_interp(hw), digits=3)
    end


   # FCC lattice constant [nm] - 7.722 bohr
    a0 = 0.4086

   # KPM shift and scale
    min_bandval = -7.641894./Ry2eV./2  # min, no fermi [Ha]
    max_bandval = 19.722771./Ry2eV./2  # max, no fermi [Ha]
    fermi = 0.4635./2  # Fermi energy [Ha]

    A = (max_bandval-min_bandval)./2 + min_bandval + fermi
    B = (max_bandval-min_bandval)./2 * 1.05
    println("A = ", round(A, digits=3))
    println("B = ", round(B, digits=3))


    return [onsite, first_neighbor, second_neighbor, A, B, fermi, a0, dielectric]
end
