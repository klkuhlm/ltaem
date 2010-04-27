set xrange [0:pi]
set xlabel "{/Symbol y}"
set ylabel "q_{BC}"
set xtics ("0" 0, "{/Symbol p}/4" pi/4,"{/Symbol p}/2" pi/2,"3{/Symbol p}/4" 3*pi/4,"{/Symbol p}" pi)
set size 0.7,0.7
set xrange [0:pi/2]
set yrange [0:1.1]
unset logscale y
set key inside bottom right
unset label 1
unset label 2
unset label 3
unset arrow
set label 1 "N=4" at pi/4-0.05,0.5
set label 2 "N=8" at 0.02,0.35
set arrow from 0.06,0.32 to 0.02,0.09 lw 0.5
set label 3 "N=16" at 0.14,0.43
set arrow from 0.2,0.4 to 0.04,0.06 lw 0.5
set output 'line_compare_flux.eps'
set terminal postscript eps enhanced color linewidth 2.5
plot 'line_compare_flux_4_36.out.plot' using 1:3 title "|sin({/Symbol y})|" with lines lw 2,\
     'line_compare_flux_4_36.out.plot' using 1:2 notitle 'N=4' with lines lt 2 lw 0.75,\
     'line_compare_flux_8_36.out.plot' using 1:2 notitle 'N=8' with lines lt 2 lw 0.75,\
     'line_compare_flux_12_36.out.plot' using 1:2 notitle 'N=12' with lines lt 2 lw 0.75,\
     'line_compare_flux_16_36.out.plot' using 1:2 notitle 'N=16' with lines lt 2 lw 0.75,\
     'line_compare_flux_20_36.out.plot' using 1:2 notitle 'N=20' with lines lt 2 lw 0.75, \
     'line_compare_flux_24_36.out.plot' using 1:2 notitle 'N=24' with lines lt 2 lw 0.75
#set key inside top center horizontal 
unset key
unset label 1
unset label 2
unset label 3
unset arrow
set label 1 "N=4" at pi/32,0.15
set label 2 "N=8" at pi/16,0.075
set arrow from pi/16,0.075 to 0.05,0.05  lw 0.5
set label 3 "N=16" at pi/8 + 0.01,0.05
set arrow from pi/8+0.05,0.04 to pi/8+0.05,0.005 lw 0.5
set yrange [-0.1:0.2]
set ylabel "q_{BC} - |sin({/Symbol y})|"     
set output 'line_compare_flux_error.eps'
plot 'line_compare_flux_4_36.out.plot' using 1:($2-$3) title 'N=4' with lines lt 1 lw 1.25,\
     'line_compare_flux_8_36.out.plot' using 1:($2-$3) title 'N=8' with lines lt 1 lw 1.25 ,\
     'line_compare_flux_12_36.out.plot' using 1:($2-$3) title 'N=12' with lines lt 1 lw 1.25 ,\
     'line_compare_flux_16_36.out.plot' using 1:($2-$3) title 'N=16' with lines lt 1 lw 1.25 ,\
     'line_compare_flux_20_36.out.plot' using 1:($2-$3) title 'N=20' with lines lt 1 lw 1.25 , \
     'line_compare_flux_24_36.out.plot' using 1:($2-$3) title '  N=24' with lines lt 1 lw 1.25  
