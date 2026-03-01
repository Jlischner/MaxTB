### Gnuplot material hot carrier generation --- [MAXTB]


### INPUTS ###
 material = "Au"
 Fermi_level = 7.3198
 optical_freq = 2.4
 overwrite_x = 1  # (FALSE_0;TRUE_1)
 xov_shift = 6.5
 Ef_to_zero = 0  # (FALSE_0;TRUE_1)
 title_output = "Plot_hcg_MTB"
 filename_mtb = 'hcg.dat'



### ACTUAL CODE ###

## Display choices ###

 # Set plot title
  set title sprintf("{/*3.5 Hot carrier generation rate (%s)  [MaxTB]}", material) font "Arial,5"

 # Set x-y-axis label
  set xlabel "Carrier energy (eV)" enhanced font "Arial,16" enhanced offset 0, -0.5
  set ylabel "Intensity  (arb. units)" enhanced font "Arial,16" enhanced offset -1.9, 0
  set bmargin screen 0.12 ; set lmargin screen 0.12

 # Set x-y-axis range
  set xrange [-8.0:28.0] ; set xtics font "Arial,13"
  set yrange [-0.002:0.012] ; set ytics font "Arial,13"

 # Enable/disable grid
  set grid

 # Enable/disable key (legend) and set position
  set key top right ; set key font "Arial,16"


## Plotting choices ###

 # Set line styles and colors
  set style line 1 lc rgb "orange" lw 3.0
  set style line 2 lc rgb "dark-orange" lw 3.0
  #set style line 1 lc rgb "cyan" lw 3.0
  #set style line 2 lc rgb "dark-blue" lw 3.0
  set style line 3 lc rgb "red" lw 2.0
  set style line 4 lc rgb "blue" lw 2.0
  set style line 5 lc rgb "green" lw 2.0
  set style line 6 lc rgb "dark-green" lw 2.0
  set style line 7 lc rgb "purple" lw 2.0
  set style line 8 lc rgb "cyan" lw 2.0
  set style line 9 lc rgb "slategrey" lw 2.0
  set style line 10 lc rgb "navy" lw 2.0
  set style line 11 lc rgb "dark-yellow" lw 2.0

 # Set dashed line style (Fermi energy)
  set style line 15 lc rgb "black" dashtype 2 lw 2.0

  shift_x = 0.0
  if (Ef_to_zero) {shift_x = Fermi_level}
  set arrow from (Fermi_level-shift_x), graph 0 to (Fermi_level-shift_x), graph 1 nohead ls 15 front

 # Initialization boudaries min/max
  xmin = 1e10; xmax = -1e10  ;  ymin = 1e10; ymax = -1e10
  xmin = real(system("awk 'NR==1{print $1}' ".filename_mtb))
  xmax = real(system("awk 'END{print $1}' ".filename_mtb))
  ymin = real(system("awk 'NR==1{min=$2;max=$2} {for(i=2;i<=NF;i++){if($i<min) min=$i; if($i>max) max=$i}} END{print min}' " . filename_mtb))
  ymax = real(system("awk 'NR==1{min=$2;max=$2} {for(i=2;i<=NF;i++){if($i<min) min=$i; if($i>max) max=$i}} END{print max}' " . filename_mtb))
  rangeX = xmax - xmin  ;  rangeY = ymax - ymin
  percent_padX = 0.02  ; percent_padY = 0.05
  padX = (rangeX == 0) ? 1 : rangeX * percent_padX  ;  padY = (rangeY == 0) ? 1 : rangeY * percent_padY
  set xrange [xmin - padX - shift_x : xmax + padX - shift_x]  ;  set yrange [ymin - padY : ymax + padY]
  if (overwrite_x) {set xrange [Fermi_level - xov_shift - shift_x : Fermi_level + xov_shift - shift_x]}

 # Filling zone added
  set object 1 rectangle \
      from (Fermi_level - optical_freq - shift_x),(ymin - padY) to (Fermi_level + optical_freq - shift_x),(ymax + padY) \
      fc rgb "purple" fillstyle solid 0.3 \
      noborder \
      behind


## Actual data plotting ###
  set output
  set format y "%.3f"
  plot filename_mtb using ($1-shift_x):($2) with lines ls 1 title "electrons"
  replot filename_mtb using ($1-shift_x):($3) with lines ls 2 title "holes"
#  replot filename_mtb using ($1-shift_x):($4) with lines ls 3 title "2*E_{points} opt. mat."
#  replot filename_mtb using ($1-shift_x):($5) with lines ls 4 title "2*E_{points} integrat."
#  replot filename_mtb using ($1-shift_x):($6) with lines ls 5 title "70% polys 1^{st} index"
#  replot filename_mtb using ($1-shift_x):($7) with lines ls 6 title "70% polys 2^{nd} index"
#  replot filename_mtb using ($1-shift_x):($8) with lines ls 7 title "70% polys both ind."
#  replot filename_mtb using ($1-shift_x):($9) with lines ls 8 title "2*gph"
#  replot filename_mtb using ($1-shift_x):($10) with lines ls 9 title "2*Temp"
#  replot filename_mtb using ($1-shift_x):($11) with lines ls 10 title "E_{F}+1"
#  replot filename_mtb using ($1-shift_x):($12) with lines ls 11 title "hv+1"

  set terminal pdfcairo size 16cm,10cm enhanced font "Arial,10"
  fileout = sprintf("%s__%s.pdf", title_output, material)
  set output fileout
  replot
  set terminal pngcairo size 800,600 enhanced font "Arial,10"
  fileout = sprintf("%s__%s.png", title_output, material)
  set output fileout
  replot
  set output


 # Uncomment the following line if you want to pause the script
 # pause -1 "Hit any key to continue"
