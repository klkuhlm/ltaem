set logscale x
set logscale y
set yrange [1.0E-5:100]
set xlabel 't'
set ylabel 's'
set key bottom right
set terminal postscript eps enhanced color linewidth 2.0
set output 'storage_comparison.eps'
plot 'test_theis' using 1:2 with lines title 'Theis' lw 2, \
     'test_finite' using 1:2 with lines title 'r_w' lw 2, \
     'test_stor' using 1:2 with lines title 'storage' lw 2, \
     'test_skin' using 1:2 with lines title 'skin + storage inside' lw 2, \
     'test_skin2' using 1:2 with lines title 'skin + storage outside' lw 2, \
     'test_skin3' using 1:2 with lines title 'skin w/o storage' lw 2
