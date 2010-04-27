##set style line 1 linetype 1 linecolor rgb "red" linewidth 4
##set style line 2 linetype 2 linecolor rgb "red" linewidth 4
##set style line 3 linetype 1 linecolor rgb "orange" linewidth 4
##set style line 4 linetype 2 linecolor rgb "orange" linewidth 4
##set style line 5 linetype 1 linecolor rgb "black" linewidth 4
##set style line 6 linetype 2 linecolor rgb "black" linewidth 4
##set style line 7 linetype 1 linecolor rgb "green" linewidth 4 
##set style line 8 linetype 2 linecolor rgb "green" linewidth 4
##set style line 9 linetype 1 linecolor rgb "cyan" linewidth 4
##set style line 10 linetype 2 linecolor rgb "cyan" linewidth 4
##set style line 11 linetype 1 linecolor rgb "blue" linewidth 4
##set style line 12 linetype 2 linecolor rgb "blue" linewidth 4
##set style line 13 linetype 1 linecolor rgb "violet" linewidth 4
##set style line 14 linetype 2 linecolor rgb "violet" linewidth 4
set style line 1 linetype 1 linewidth 1
set style line 2 linetype 2 linewidth 4
set style line 3 linetype 1 linewidth 4
set style line 4 linetype 2 linewidth 4
set style line 5 linetype 1 linewidth 4
set style line 6 linetype 2 linewidth 4
set style line 7 linetype 1 linewidth 1 
set style line 8 linetype 2 linewidth 4
set style line 9 linetype 1 linewidth 4
set style line 10 linetype 2 linewidth 4
set style line 11 linetype 1 linewidth 4
set style line 12 linetype 2 linewidth 4
set style line 13 linetype 1 linewidth 4
set style line 14 linetype 2 linewidth 4
set format y "10^{%T}"
set format x "10^{%T}"
set logscale x
set logscale y
set yrange [0.005:10]
set xrange [1.0E-4:400]
unset key
set xlabel 't/r^2' 0.0,0.25
set ylabel 's' 2.0,0.0
set label 1 'E_{1}(t/4)' at 0.0003,1
set label 2 'circles B'  at 0.002,0.09 rotate by 33
set label 3 'circles A'   at 0.003,0.015 rotate by 60
set label 4 'uniform B'  at 0.105,0.025 rotate by 48
set label 5 'uniform A'   at 0.01,0.015 rotate by 53
set terminal postscript eps enhanced mono 'Times-Roman' 20
set size 0.85,0.85
set output 'leaky_circles_time.eps'
plot 'h1'       index 0 using ($1/0.5):(-$2) with lines ls 2, \
     'h1_leak' index 0 using  ($1/0.5):(-$2) with lines ls 1, \
     'h1'      index 1 using  ($1/1.75):(-$2) with lines ls 8, \
     'h1_leak' index 1 using  ($1/1.75):(-$2) with lines ls 7, \
     'h1_bench' index 0 using ($1/0.5):(-$2) with lines ls 5

