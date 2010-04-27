set dgrid3d 30,30,6
set hidden3d
set contour base
set cntrparam level incr 0.5,0.1,1.2
unset key
set xlabel 'y'
set ylabel 'x'
#set xrange [-1.5:1.0]
#set yrange [-1.0:1.0]
set size ratio -1
set terminal postscript eps enhanced color linewidth 3
#set title 'Mathieu line'
unset title
set view ,75
set output 'line_surface_jem_01.eps'
splot "bench_line_compare_jem_250_16_32.out" using 1:2:6 with lines
#
##set title 'Bessel point -> line'
##set output 'line_surface_jem_02.eps'
##splot 'line_bench.out' using 1:2:3 with lines
#
set title 'difference between solutions'
set output 'line_surface_jem_03.eps'
splot "bench_line_compare_jem_250_16_32.out" using 1:2:(($3-$6)) with lines
