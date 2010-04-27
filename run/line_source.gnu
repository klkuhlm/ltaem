set xrange [-1:1]
set xtics -1,0.5,1
set xlabel "X"
set yrange [-1:1]
set ytics -1,0.5,1
set ylabel "Y"
set zrange [0.0:0.5]
set ztics 0.0,0.1,0.5
unset zlabel
##set zlabel '{/Symbol F}'
unset key
set dgrid3d 20,20,2
set hidden3d
set contour
set cntrparam levels incremental 0,0.1,0.5
set output 'line_source_surface.eps'
set term postscript eps enhanced color linewidth 2.5 'Times-Roman' 22
splot 'line_mathieu.out' u 1:2:3 w l
set output 'line_source_vector_field.eps'
plot 'line_mathieu.out' u 1:2:(-$4/8):(-$5/8) w vec
