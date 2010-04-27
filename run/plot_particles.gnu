set xrange [-2:3]
set yrange [-2:2]
set autoscale x
set autoscale y
unset key
set terminal postscript eps enhanced color linewidth 1.5
set output 'particle_grid_results.eps'
set size ratio -1
set xlabel 'X'
set ylabel 'Y'
set grid 
plot 'circles.dat' with lines, \
     'well.dat' using 2:3 with points, \
     'path_grid' using 2:3:(sprintf("%6.4g", $1)) with labels, \
     'path_grid' using 2:3 with linespoints 
