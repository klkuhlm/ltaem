set style line 1 linetype 1 linecolor rgb "red" linewidth 3
set style line 2 linetype 2 linecolor rgb "red" linewidth 0.5
set style line 3 linetype 1 linecolor rgb "orange" linewidth 3
set style line 4 linetype 2 linecolor rgb "orange" linewidth 0.5
set style line 5 linetype 1 linecolor rgb "black" linewidth 3
set style line 6 linetype 2 linecolor rgb "black" linewidth 0.5
set style line 7 linetype 1 linecolor rgb "green" linewidth 3 
set style line 8 linetype 2 linecolor rgb "green" linewidth 0.5
set style line 9 linetype 1 linecolor rgb "cyan" linewidth 3
set style line 10 linetype 2 linecolor rgb "cyan" linewidth 0.5
set style line 11 linetype 1 linecolor rgb "blue" linewidth 3
set style line 12 linetype 2 linecolor rgb "blue" linewidth 0.5
set style line 13 linetype 1 linecolor rgb "violet" linewidth 3
set style line 14 linetype 2 linecolor rgb "violet" linewidth 0.5
set format y "%0.1f"
set format x "10^{%T}"
set logscale x
set logscale y
set xlabel 't_D
set ylabel 'h_D
set xrange [1.0E-2:1.0E+11]
set yrange [1.0E-1:50]
#unset label 1
#unset label 2
#unset label 3
#unset label 4
#unset label 5
#unset label 6
set label 1 'E(S_s)' at 500,11
set label 2 'E(S_y)' at 5.0E5,0.2
##set label 3 '{/Symbol b}=1' at 10,0.25
##set label 4 '{/Symbol b}=0.5' at 100,0.8
##set label 5 '{/Symbol b}=0.1' at 1000,2.0
##set label 6 '{/Symbol b}=0.01' at 10000,4.5
##set label 7 '{/Symbol b}=0.001' at 1.0E+5,7.5
set output 'neuman_test.eps'
unset key 
set terminal postscript eps enhanced color 'Times Roman' 20
plot 'neuman_test.out' using 1:2 title 'E(S_S)' with lines ls 5, \
     'neuman_test.out' using 1:5 title 'N72?' with lines ls 1, \
     'neuman_test.out' using 1:10 title 'E(S_Y)' with lines ls 5   
