set dgrid3d 21,21,1
set contour base
unset surface
set view 0,0,1,1
set xrange [-1.0:1.0]
set yrange [-1.0:1.0]
unset key
set cntrparam levels incr -1.0,0.05,1.0
set table 'contours.tab'
splot 'test_ellipses' u 1:2:3 with line
#
#
set term postscript eps enhanced color 'Times Roman' 20
set output 'check_ellipses_contours.eps'
plot 'test_ellipses' u 1:2:($4/5.0):($5/5.0) with vector, \
     'contours.tab' with line, 'ellipse.dat' with line
#
#
set output 'check_ellipses_surface.eps'
set term postscript eps enhanced color 'Times Roman' 20
set surface
set hidden3d
set cntrparam levels incr -0.8,0.1,0.8
set view 60,30,1,1
splot 'test_ellipses' u 1:2:3 with line
unset output
set term X11
replot
