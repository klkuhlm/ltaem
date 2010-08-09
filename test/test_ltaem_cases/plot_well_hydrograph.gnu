set logscale x
set logscale y
set xlabel 'time'
set ylabel 'abs(head)'
set yrange [1.0E-10:0.1]
set terminal postscript eps color enhanced
set output 'well_hydrograph_head.eps'
unset key
plot 'well_hydrograph.dat' u 1:(abs($2)) w l
unset output
#
set xlabel 'time'
set yrange [1.0E-10:1]
set ylabel 'abs(velocity)'
set key top left
set output 'well_hydrograph_vectors.eps'
plot 'well_hydrograph.dat' u 1:(abs($3)) t 'v_x'w l, \
     'well_hydrograph.dat' u 1:(abs($4)) t 'v_y'w l
unset output
