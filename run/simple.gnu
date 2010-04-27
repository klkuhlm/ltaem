set dgrid3d 20,20,2
set hidden3d
#unset surface
set contour base
set output 'check.eps'
set terminal postscript eps enhanced color
set xlabel 'X'
set ylabel 'Y'
#set view 0,0,1,1
splot 'arc_source.out' using 1:2:3 with lines
set output 'arc.eps'
plot  'arc.bdry' i 0 w lp, 'arc.bdry' i 1 w lp
set output 'check_vec.eps'
plot 'arc_source.out' using 1:2:($4/5):($5/5) with vector
