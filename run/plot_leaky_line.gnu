set style line 1 linetype 1 linecolor rgb "red" linewidth 4
set style line 2 linetype 2 linecolor rgb "red" linewidth 0.5
set style line 3 linetype 1 linecolor rgb "orange" linewidth 4
set style line 4 linetype 2 linecolor rgb "orange" linewidth 0.5
set style line 5 linetype 1 linecolor rgb "black" linewidth 6
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
set label 'away from source' at 0.003,0.03
set arrow from 0.02,0.025 to 0.03,0.006
set key bottom right
set output 'line_leaky1_dist.eps'
plot 'line_leaky_1_0.00_0.08_0.05_5.0.out' title 'y=0.0' with lines, \
     'line_leaky_1_0.25_0.08_0.05_5.0.out' title 'y=0.25' with lines, \
     'line_leaky_1_1.00_0.08_0.05_5.0.out' title 'y=1.0' with lines, \
     'line_leaky_1_2.00_0.08_0.05_5.0.out' title 'y=2.0' with lines, \
     'line_leaky_1_4.00_0.08_0.05_5.0.out' title 'y=4.0' with lines
set output 'line_leaky2_dist.eps'
plot 'line_leaky_2_0.00_0.08_0.05_5.0.out' title 'y=0.0' with lines, \
     'line_leaky_2_0.25_0.08_0.05_5.0.out' title 'y=0.25' with lines, \
     'line_leaky_2_1.00_0.08_0.05_5.0.out' title 'y=1.0' with lines, \
     'line_leaky_2_2.00_0.08_0.05_5.0.out' title 'y=2.0' with lines, \
     'line_leaky_2_4.00_0.08_0.05_5.0.out' title 'y=4.0' with lines
set output 'line_leaky3_dist.eps'
plot 'line_leaky_3_0.00_0.08_0.05_5.0.out' title 'y=0.0' with lines, \
     'line_leaky_3_0.25_0.08_0.05_5.0.out' title 'y=0.25' with lines, \
     'line_leaky_3_1.00_0.08_0.05_5.0.out' title 'y=1.0' with lines, \
     'line_leaky_3_2.00_0.08_0.05_5.0.out' title 'y=2.0' with lines, \
     'line_leaky_3_4.00_0.08_0.05_5.0.out' title 'y=4.0' with lines     
#
set output 'line_leaky_normalized.eps'
plot 'line_leaky_1_0.25_0.08_0.05_5.0.out' using ($1/(0.25*0.25)):2 title 'y=0.25' with lines, \
     'line_leaky_1_1.00_0.08_0.05_5.0.out' using ($1/(1*1)):2 title 'y=1.0' with lines, \
     'line_leaky_1_2.00_0.08_0.05_5.0.out' using ($1/(2*2)):2 title 'y=2.0' with lines, \
     'line_leaky_1_4.00_0.08_0.05_5.0.out' using ($1/(4*4)):2 title 'y=4.0' with lines
set output 'line_leaky_all_types.eps'
set yrange [0.0003:0.05]
set xrange [1.0E-5:100]
unset label
unset arrow
set label 'at source (y=0)' at 0.0001,0.017 rotate by 17
set label 'y=4' at 3.0E-5,0.001 rotate by 65
plot 'line_leaky_0_0.00_0.08_0.05_1.0.out' title 'non-leaky' with lines ls 5, \
     'line_leaky_0_4.00_0.08_0.05_1.0.out' notitle with lines ls 5, \
     'line_leaky_1_0.00_0.08_0.05_1.0.out' title 'type I' with lines linetype 2 linewidth 3, \
     'line_leaky_2_0.00_0.08_0.05_1.0.out' title 'type II' with lines linetype 3 linewidth 3, \
     'line_leaky_3_0.00_0.08_0.05_1.0.out' title 'b_2 {/Symbol \256 \245}' with lines linetype 4 linewidth 3, \
     'line_leaky_1_4.00_0.08_0.05_1.0.out' notitle with lines linetype 2 linewidth 3, \
     'line_leaky_2_4.00_0.08_0.05_1.0.out' notitle with lines linetype 3 linewidth 3, \
     'line_leaky_3_4.00_0.08_0.05_1.0.out' notitle with lines linetype 4 linewidth 3
