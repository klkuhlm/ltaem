set terminal gif animate delay 100
set output 'nested-ic.gif'
set dgrid3d 80,80,2
set zrange [-1:2]
do for [i=0:15] {
  splot 'initial_Condition_contour.dat' index (i) using 1:2:3 with lines
}
unset output
