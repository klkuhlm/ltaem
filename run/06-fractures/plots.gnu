set dgrid3d 40,40,4
set hidden3d
set contour
set size 0.8
unset key
set view 60,60

set terminal postscript eps enhanced color
set output 'frac01.eps'
splot 'fractures_contour.dat' index 0 using 1:2:3 with lines

set output 'frac02.eps'
splot 'fractures_contour.dat' index 1 using 1:2:3 with lines

set output 'frac03.eps'
splot 'fractures_contour.dat' index 2 using 1:2:3 with lines

set output 'frac04.eps'
splot 'fractures_contour.dat' index 3 using 1:2:3 with lines

set logscale x
set logscale y
set xlabel 'time'
set ylabel 'flowrate'
set output 'flowrate.eps'
plot 'fractures_timeseries.Q' using 1:2 with lines

unset output
unset terminal
