set xrange [-1.7:2.6]
set yrange [-2.1:2.1]
unset key
unset xtics
unset ytics
set cbrange [0:85]
color(t)=t/84.0
set size ratio -1
set terminal postscript eps enhanced color linewidth 3.0 
set output 'grid_particles_agu08.eps'
plot 'circles.dat' using 1:2 with lines, \
     'path_grid' using 2:3:1 with points pt 7 ps 1.0 lc palette z, \
     'well.dat' using 2:3 with points pt 7 ps 2, \
     'path_grid' using 2:3 with lines lt 0 lw 0.5
