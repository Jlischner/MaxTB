### Gnuplot material band structure --- [MAXTB]
#   (run "python extract_bands_MTB.py" first)

### INPUTS ###
 material = "Au"
 Fermi_level = 7.3198
 Nb = 9
 Nk = 9
 shift_k = 30  +(-1)
 array lb[Nk] = [ "G" , "X" , "W" , "K" , "G" , "L" , "U" , "W" , "L" ]
 optical_freq = 2.4
 overwrite_y = 1  # (FALSE_0;TRUE_1)
 yov_shift = 6.5
 title_output = "Plot_bands_MTB"
 filename_mtb = 'bands_parsed_mtb.dat'



### ACTUAL CODE ###

## Display choices ###

 # Set plot title
  set title sprintf("{/*3.5 Band structure (%s)  [MaxTB]}", material) font "Arial,5"

 # Set x-y-axis label
  set xlabel " " enhanced font "Arial,16" enhanced offset 0, -1.5
  set ylabel "E_{v} (eV)" enhanced font "Arial,16" enhanced offset -1.5, 0
  set bmargin screen 0.15 ; set lmargin screen 0.15

 # Set x-y-axis range
  set xrange [0.0:1.0] ; set noxtics
  set yrange [-10.0:20.0] ; set ytics font "Arial,13"

 # Enable/disable grid
  set grid

 # Enable/disable key (legend) and set position
  set key horizontal
  set key top center ; set key font "Arial,16"
  set key at screen 0.54,0.07


## Plotting choices ###

 # Set line styles and colors
  set style line 1 lc rgb "dark-green" lw 3.0

 # Set dashed line style (high-symmetry kpts)
  set style line 15 lc rgb "black" dashtype 2 lw 0.5


 # High-symmetry kpt positions in data
  array kpts[Nk]
  kpts[1] = 0
  do for [i=2:Nk] {\
    kpts[i] = kpts[i-1] + shift_k}

 # Read high-symmetry kpt values in data
  array xvals[Nk]
  do for [i=1:Nk] {\
      xvals[i] = real(system(sprintf("awk 'NR==%d{print $1}' %s", kpts[i]+1, filename_mtb)))}
  maxvalx_mtb = xvals[Nk]
  do for [i=1:Nk] {\
      xvals[i] = xvals[i] / xvals[Nk]}

 # Set kpt labels and placements
  do for [i=1:Nk] {\
      set label i lb[i] at xvals[i], graph -0.03 center font "Arial,16"
      if (i > 1 && i < Nk) {set arrow from xvals[i], graph 0 to xvals[i], graph 1 nohead ls 15}}

 # Initialization boudaries min/max
  xmin = 1e10; xmax = -1e10  ;  ymin = 1e10; ymax = -1e10
  ymin = real(system("awk -v N=".(Nb+1)." '{for(i=2;i<=N;i++){if(NR==1 && i==2){min=$i; max=$i}; \
      if($i<min) min=$i; if($i>max) max=$i}} END{print min}' ".filename_mtb))
  ymax = real(system("awk -v N=".(Nb+1)." '{for(i=2;i<=N;i++){if(NR==1 && i==2){min=$i; max=$i}; \
      if($i<min) min=$i; if($i>max) max=$i}} END{print max}' ".filename_mtb))
  rangeY = ymax - ymin
  percent_padY = 0.05
  padY = (rangeY == 0) ? 1 : rangeY * percent_padY
  set xrange [xvals[1] : xvals[Nk]]  ;  set yrange [ymin - padY - Fermi_level : ymax + padY - Fermi_level]
  if (overwrite_y) {set yrange [-yov_shift : yov_shift]}

 # Filling zone added
  set object 1 rectangle \
      from (xvals[1]),(-optical_freq) to (xvals[Nk]),(optical_freq) \
      fc rgb "purple" fillstyle solid 0.3 \
      noborder \
      behind


## Actual data plotting ###
  set output
  set format y "%.1f"
  plot filename_mtb using ($1/maxvalx_mtb):($2-Fermi_level) with lines ls 1 title "MaxTB"
  replot for [i=3:Nb+1] filename_mtb using ($1/maxvalx_mtb):(column(i)-Fermi_level) with lines ls 1 notitle

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
