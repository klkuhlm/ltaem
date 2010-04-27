set style line 1 linetype 1 linecolor rgb "red" linewidth 4
set style line 2 linetype 2 linecolor rgb "red" linewidth 0.5
set style line 3 linetype 1 linecolor rgb "orange" linewidth 4
set style line 4 linetype 2 linecolor rgb "orange" linewidth 0.5
set style line 5 linetype 1 linecolor rgb "red" linewidth 6
set style line 6 linetype 2 linecolor rgb "black" linewidth 0.5
set style line 7 linetype 1 linecolor rgb "green" linewidth 4 
set style line 8 linetype 2 linecolor rgb "green" linewidth 0.5
set style line 9 linetype 1 linecolor rgb "cyan" linewidth 4
set style line 10 linetype 2 linecolor rgb "cyan" linewidth 0.5
set style line 11 linetype 1 linecolor rgb "blue" linewidth 4
set style line 12 linetype 2 linecolor rgb "blue" linewidth 0.5
set style line 13 linetype 1 linecolor rgb "violet" linewidth 4
set style line 14 linetype 2 linecolor rgb "violet" linewidth 0.5
set logscale x
set logscale y
set format y "10^{%T}"
set format x "10^{%T}"
set terminal postscript eps enhanced color linewidth 2.0
set size 0.7,0.7
set xlabel 't'
set ylabel 's' 2,0
set yrange [0.0001:0.04]
set label 'y=4' at 3.0E-5,2.0E-4
set label 'y=2' at 2.0E-5,1.25E-3 rotate by 55
set label 'y=1' at 1.2E-5,4.0E-3 rotate by 40
set label 'y=0.25' at 1.0E-4,3.0E-2
set arrow from 2.0E-4,2.7E-2 to 1.0E-4,1.1E-2 front
set label 'y=0' at 1.25E-5,1.2E-2 rotate by 15 front
set label 'away from line source' at 0.004,0.03
set arrow from 0.02,0.025 to 0.03,0.006 front
set key bottom right
unset key
set output 'line_unconf_dist.eps'
plot 'line_unconf_0.00_0.25_0.80.out' title 'y=0.0' with lines ls 5, \
     'line_unconf_0.25_0.25_0.80.out' title 'y=0.25' with lines ls 7, \
     'line_unconf_1.00_0.25_0.80.out' title 'y=1.0' with lines ls 9, \
     'line_unconf_2.00_0.25_0.80.out' title 'y=2.0' with lines ls 11, \
     'line_unconf_4.00_0.25_0.80.out' title 'y=4.0' with lines ls 13
set output 'line_unconf_early_late.eps'
set yrange [0.001:0.05]
set xrange [1.0E-5:100]
unset label
unset arrow
unset key
set label "{/Symbol b}=0.1" at 1.0E-3,5.0E-3
set label "{/Symbol b}=0.01" at 4.0E-3,9.0E-3
set label "{/Symbol b}=0.001" at 1.0E-1,1.8E-2
set label 'S_S' at 0.001,0.014 rotate by 22
set label 'S_y' at 0.004,0.0016 rotate by 59
plot 'line_confSS_1.00_0.25_0.8.out' notitle with lines ls 5, \
     'line_confSy_1.00_0.25_0.8.out' notitle with lines ls 5, \
     'line_unconf_1.00_0.25_8.00.out' notitle  with lines ls 7, \
      'line_unconf_1.00_0.25_0.80.out' notitle  with lines ls 9, \
     'line_unconf_1.00_0.25_0.08.out' notitle  with lines ls 11
