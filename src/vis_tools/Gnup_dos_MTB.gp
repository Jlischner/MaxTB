### Gnuplot material density of states --- [MAXTB]


### INPUTS ###
 material = "Au"
 Fermi_level = 7.3198
 optical_freq = 2.4
 overwrite_x = 1  # (FALSE_0;TRUE_1)
 xov_shift = 6.5
 Ef_to_zero = 1  # (FALSE_0;TRUE_1)
 title_output = "Plot_dos_MTB"
 filename_mtb = 'dos.dat'



### ACTUAL CODE ###

## Display choices ###

 # Set plot title
  set title sprintf("{/*3.5 Density of states (%s)  [MaxTB]}", material) font "Arial,5"

 # Set x-y-axis label
  set xlabel "State energy (eV)" enhanced font "Arial,16" enhanced offset 0, -0.5
  set ylabel "Intensity  (arb. units)" enhanced font "Arial,16" enhanced offset -1.9, 0
  set bmargin screen 0.12 ; set lmargin screen 0.12

 # Set x-y-axis range
  set xrange [-8.0:28.0] ; set xtics font "Arial,13"
  set yrange [-0.02:0.18] ;  set ytics font "Arial,13"

 # Enable/disable grid
  set grid

 # Enable/disable key (legend) and set position
  set key horizontal
  set key top center ; set key font "Arial,16"
  set key at screen 0.20,0.06


## Plotting choices ###

 # Set line styles and colors
  set style line 1 lc rgb "dark-green" lw 3.0

 # Set dashed line style (Fermi energy)
  set style line 15 lc rgb "black" dashtype 2 lw 2.0

  shift_x = 0.0
  if (Ef_to_zero) {shift_x = Fermi_level}
  set arrow from (Fermi_level-shift_x), graph 0 to (Fermi_level-shift_x), graph 1 nohead ls 15 front

 # Initialization boudaries min/max
  xmin = 1e10; xmax = -1e10  ;  ymin = 1e10; ymax = -1e10
  xmin = real(system("awk 'NR==1{print $1}' ".filename_mtb))
  xmax = real(system("awk 'END{print $1}' ".filename_mtb))
  ymin = real(system("awk 'NR>0{if(min==\"\") min=$2; if($2<min) min=$2} END{print min}' ".filename_mtb))
  ymax = real(system("awk 'NR>0{if(max==\"\") max=$2; if($2>max) max=$2} END{print max}' ".filename_mtb)) 
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
  set format y "%.2f"
  plot filename_mtb using ($1-shift_x):($2) with lines ls 1 title "MaxTB"

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
