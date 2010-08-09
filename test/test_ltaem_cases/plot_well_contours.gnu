set term postscript eps enhanced color
set xlabel 'x'
set ylabel 'y'
set xlabel 'head'
set xrange [0:1.1]
set yrange [0:1.1]
unset key
set dgrid3d 100,100,4
set contour base
set output 'well_contours_head.eps'
splot 'well_contours.xyz' using 1:2:3 w l
unset output
#
set output 'well_contours_vectors.eps'
set size ratio -1
plot 'well_contours.xyz' using 1:2:(2*$4):(2*$5) w vec
unset output
