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
set style line 1 linetype 2 linewidth 4
set style line 2 linetype 2 linewidth 1
set style line 3 linetype 3 linewidth 2
set style line 4 linetype 2 linewidth 4
set style line 5 linetype 1 linewidth 2
set style line 6 linetype 2 linewidth 4
set style line 7 linetype 4 linewidth 4 
set style line 8 linetype 2 linewidth 1
set style line 9 linetype 5 linewidth 4
set style line 10 linetype 2 linewidth 4
set style line 11 linetype 6 linewidth 2
set style line 12 linetype 2 linewidth 4
set style line 13 linetype 1 linewidth 4
set style line 14 linetype 2 linewidth 4
set format y "10^{%T}"
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
set label 2 "leaky case I" at 9.0E+6,1.6
set arrow from 2.0E+10,1.6 to 2.0E+11,2.6
set arrow from 2.0E+10,1.6 to 1.0E+11,0.85
#set label 5 'b_2=1' at 5.0E+4,2.05 font "Times-Roman, 16"
#set label 6 'b_2=1/10' at 2.0E+2,0.65 font "Times-Roman, 16"
set label 3 "leaky case II" at 9.0E+3,30
set arrow from 3.0E+7,30 to 7.0E+8,18
set arrow from 3.0E+7,30 to 7.0E+8,15
#set label 7 'b_2=1' at 3.0E+2,18.0 font "Times-Roman, 16"
set arrow from 1.0E+3,15.0 to 3.0E+4,5.0
#set label 8 'b_2=1/10' at 70,2 rotate by 45 font "Times-Roman, 16"
set label 4 "leaky b_2 {/Symbol \256 \245}" at 5.0E+7,6 rotate by 10
set output 'leaky_method_comparison.eps'
set terminal postscript eps enhanced mono 'Times-Roman' 22
set size 0.75,0.75
##
## change bounding box command in eps file to 60 60 315 229
##
plot 'leaky_unit.out' using 1:2  with lines ls 5, \
     'leaky_unit.out' using 1:3  with lines ls 1, \
     'leaky_unit.out' using 1:4  with lines ls 1, \
     'leaky_unit.out' using 1:7  with lines ls 9  #, \
 #    'leaky_thin.out' using 1:3  with lines ls 2, \
 #    'leaky_thin.out' using 1:4  with lines ls 2

