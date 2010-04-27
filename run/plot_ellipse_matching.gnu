set dgrid3d 30,30,4
set hidden3d
unset key
set contour base
set cntrparam levels incr 0,0.1,2
set xtics -1.5,0.5,1.5
set ytics -1.5,0.5,1.5
set ztics 0.0,0.1,1.2
set size ratio -1
set xlabel 'X'
set ylabel 'Y'
set terminal postscript eps enhanced color linewidth 2.0
set output 'ellipse_well_highk.eps'
splot 'line_check03.out' using 1:2:3 with lines
set cntrparam levels incr 0,0.05,1
set output 'ellipse_well_lowk.eps'
splot 'line_check04.out' using 1:2:3 with lines
