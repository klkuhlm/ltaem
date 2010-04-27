set style line 1 linetype 1 linecolor rgb "red" linewidth 4
set style line 2 linetype 2 linecolor rgb "red" linewidth 4
set style line 3 linetype 1 linecolor rgb "orange" linewidth 4
set style line 4 linetype 2 linecolor rgb "orange" linewidth 4
set style line 5 linetype 1 linecolor rgb "black" linewidth 4
set style line 6 linetype 2 linecolor rgb "black" linewidth 4
set style line 7 linetype 1 linecolor rgb "green" linewidth 4 
set style line 8 linetype 2 linecolor rgb "green" linewidth 4
set style line 9 linetype 1 linecolor rgb "cyan" linewidth 4
set style line 10 linetype 2 linecolor rgb "cyan" linewidth 4
set style line 11 linetype 1 linecolor rgb "blue" linewidth 4
set style line 12 linetype 2 linecolor rgb "blue" linewidth 4
set style line 13 linetype 1 linecolor rgb "violet" linewidth 4
set style line 14 linetype 2 linecolor rgb "violet" linewidth 4
set format y "%0.1f"
set format x "10^{%T}"
set logscale x
set logscale y
set xlabel 't_D' 0.0,0.25
set ylabel 'h_D' 2.5,0.0 
set xrange [1.0E-1:1.0E+12]
set yrange [1.0E-1:50]
unset key
#unset arrow
set label 1 "E_{1}(t_{D}/4)" at 0.3,8
set label 2 "leaky case I" at 7.0E+6,1.0
set label 3 "leaky case II" at 9.0E+3,30
set arrow from 1.0E+5,22  to 5.0E+6,12 
set label 4 "leaky b_2 {/Symbol \256 \245}" at 5.0E+7,6 rotate by 10
set output 'leaky_method_comparison2.eps'
set terminal postscript eps enhanced color 'Times-Roman' 22
set size 0.75,0.75
##
## change bounding box command in eps file to 60 60 315 229
##
plot 'beta_10-3_b1.out' using 1:2  with lines ls 5, \
     'beta_10-3_b1.out' using 1:3  with lines ls 1, \
     'beta_10-3_b1.out' using 1:4  with lines ls 7, \
     'beta_10-3_b1.out' using 1:7  with lines ls 9, \
     'beta_10-3_b0.1.out' using 1:3  with lines ls 2, \
     'beta_10-3_b0.1.out' using 1:4  with lines ls 8

