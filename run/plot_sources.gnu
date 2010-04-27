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
unset label 1
unset label 2
unset label 3
unset label 4
unset label 5
unset label 6
set label 1 'E(S_s)' at 900,11
set label 2 'E(S_y)' at 1.0E6,0.2
set label 3 '{/Symbol b}=1' at 100,0.25
set label 4 '{/Symbol b}=0.5' at 100,0.8
set label 5 '{/Symbol b}=0.1' at 1000,1.7
set label 6 '{/Symbol b}=0.01' at 10000,5
set output 'unconfined_source.eps'
unset key 
set terminal postscript eps enhanced color 'Times Roman' 20
plot 'unconfined_ratio_1.out' using 1:2 title 'E(S_S)' with lines ls 1, \
     'unconfined_ratio_1.out' using 1:5 title '{/Symbol b}=1' with lines lw 2, \
     'unconfined_ratio_2.out' using 1:5 title '{/Symbol b}=0.5' with lines lw 2, \
     'unconfined_ratio_10.out' using 1:5 title '{/Symbol b}=0.1' with lines lw 2, \
     'unconfined_ratio_100.out' using 1:5 title '{/Symbol b}=0.01' with lines lw 2, \
     'unconfined_ratio_1.out' using 1:9 title 'E(S_Y)' with lines        
