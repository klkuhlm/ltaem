set style line 1 linetype 1 linecolor rgb "red" linewidth 2
set style line 2 linetype 2 linecolor rgb "red" linewidth 0.5
set style line 15 linetype 3 linecolor rgb "red" linewidth 2
set style line 3 linetype 1 linecolor rgb "orange" linewidth 2
set style line 4 linetype 2 linecolor rgb "orange" linewidth 0.5
set style line 16 linetype 3 linecolor rgb "orange" linewidth 2
set style line 5 linetype 1 linecolor rgb "black" linewidth 2
set style line 6 linetype 2 linecolor rgb "black" linewidth 0.5
set style line 7 linetype 1 linecolor rgb "green" linewidth 2 
set style line 8 linetype 2 linecolor rgb "green" linewidth 0.5
set style line 17  linetype 3 linecolor rgb "green" linewidth 2
set style line 9 linetype 1 linecolor rgb "cyan" linewidth 2
set style line 10 linetype 2 linecolor rgb "cyan" linewidth 0.5
set style line 18 linetype 3 linecolor rgb "cyan" linewidth 2
set style line 11 linetype 1 linecolor rgb "blue" linewidth 2
set style line 12 linetype 2 linecolor rgb "blue" linewidth 0.5
set style line 19 linetype 3 linecolor rgb "blue" linewidth 2
set style line 13 linetype 1 linecolor rgb "violet" linewidth 2
set style line 14 linetype 2 linecolor rgb "violet" linewidth 0.5
set style line 20 linetype 3 linecolor rgb "violet" linewidth 2
set format y "%0.1f"
set format x "10^{%T}"
set logscale x
set logscale y
set xlabel 't_D' 0.0,0.0
set ylabel 'h_D' 2.0,0.0
set xrange [1.0E-1:1.0E+7]
set yrange [1.0E-1:10]
unset key
set label 1 'E(S_s)' at 1.0,4
set label 2 'E(S_y)' at 6.0E5,0.2
set label 3 '{/Symbol b}=1' at 10,0.25
set label 4 '{/Symbol b}=0.5' at 100,0.8
set label 5 '{/Symbol b}=0.1' at 1000,2.0
set label 6 '{/Symbol b}=0.01' at 10000,4.5
set label 7 '{/Symbol b}=0.001' at 1.0E+5,7.5
##set title 'comparison of unconfined methods r=1'
set output 'unconfined_r1_comparison.eps'
set terminal postscript eps enhanced color 'Times-Roman' 24
plot 'r1unconf_beta2.out' using 1:9  title 'E(S_s)' with lines ls 5, \
     'r1unconf_beta1.out' using 1:10  title 'E(S_y)' with lines ls 5, \
     'r1unconf_beta1.out' using 1:8 title '{/Symbol b}=1' with lines ls 1, \
     'r1unconf_beta1.out' using 1:5 notitle  with lines ls 2, \
     'dimensionless_neuman_beta1.dat' using 1:2 notitle with lines ls 15, \
     'r1unconf_beta2.out' using 1:8 title '{/Symbol b}=2' with lines ls 3, \
     'r1unconf_beta2.out' using 1:5 notitle  with lines ls 4, \
     'dimensionless_neuman_beta2.dat' using 1:2 notitle with lines ls 16, \
     'r1unconf_beta10.out' using 1:8 title '{/Symbol b}=10' with lines ls 7, \
     'r1unconf_beta10.out' using 1:5 notitle  with lines ls 8, \
     'dimensionless_neuman_beta10.dat' using 1:2 notitle with lines ls 17, \
     'r1unconf_beta100.out' using 1:8 title '{/Symbol b}=100' with lines ls 9, \
     'r1unconf_beta100.out' using 1:5 notitle  with lines ls 10, \
     'dimensionless_neuman_beta100.dat' using 1:2 notitle with lines ls 18, \
     'r1unconf_beta1000.out' using 1:8 title '{/Symbol b}=1000' with lines ls 11, \
     'r1unconf_beta1000.out' using 1:5 notitle  with lines ls 12, \
     'dimensionless_neuman_beta1000.dat' using 1:2 notitle with lines ls 19
set key bottom right
unset label 1
unset label 2
unset label 3
unset label 4
unset label 5
unset label 6
unset label 7
set title 'comparison of unconfined methods r=0.1'
set output 'unconfined_r0.1_comparison.eps'
set terminal postscript eps enhanced color 'Times-Roman' 24
plot 'r0.1unconf_beta1.out' using 1:8 title "{/Symbol b}=1" with lines ls 1, \
     'r0.1unconf_beta1.out' using 1:5 notitle  with lines ls 2, \
     'r0.1unconf_beta2.out' using 1:8 title '{/Symbol b}=2' with lines ls 3, \
     'r0.1unconf_beta2.out' using 1:5 notitle  with lines ls 4, \
     'r0.1unconf_beta10.out' using 1:8 title '{/Symbol b}=10' with lines ls 7, \
     'r0.1unconf_beta10.out' using 1:5 notitle  with lines ls 8, \
     'r0.1unconf_beta100.out' using 1:8 title '{/Symbol b}=100' with lines ls 9, \
     'r0.1unconf_beta100.out' using 1:5 notitle  with lines ls 10, \
     'r0.1unconf_beta1000.out' using 1:8 title '{/Symbol b}=1000' with lines ls 11, \
     'r0.1unconf_beta1000.out' using 1:5 notitle  with lines ls 12, \
     'r0.1unconf_beta2.out' using 1:9  title 'E(S_s)' with lines ls 5, \
     'r0.1unconf_beta2.out' using 1:10  title 'E(S_y)' with lines ls 5
set format y "10^{%T}"
set yrange [1.0E-2:50]
set title 'comparison of unconfined methods r=2'
set output 'unconfined_r2_comparison.eps'          
plot 'r2unconf_beta1.out' using 1:8 title "{/Symbol b}=1" with lines ls 1, \
     'r2unconf_beta1.out' using 1:5 notitle  with lines ls 2, \
     'r2unconf_beta2.out' using 1:8 title '{/Symbol b}=2' with lines ls 3, \
     'r2unconf_beta2.out' using 1:5 notitle  with lines ls 4, \
     'r2unconf_beta10.out' using 1:8 title '{/Symbol b}=10' with lines ls 7, \
     'r2unconf_beta10.out' using 1:5 notitle  with lines ls 8, \
     'r2unconf_beta100.out' using 1:8 title '{/Symbol b}=100' with lines ls 9, \
     'r2unconf_beta100.out' using 1:5 notitle  with lines ls 10, \
     'r2unconf_beta1000.out' using 1:8 title '{/Symbol b}=1000' with lines ls 11, \
     'r2unconf_beta1000.out' using 1:5 notitle  with lines ls 12, \
     'r2unconf_beta2.out' using 1:9  title 'E(S_s)' with lines ls 5, \
     'r2unconf_beta2.out' using 1:10  title 'E(S_y)' with lines ls 5
set format y "10^{%T}"
set yrange [1.0E-7:50]
set title 'comparison of unconfined methods r=10'
set output 'unconfined_r10_comparison.eps'          
plot 'r10unconf_beta1.out' using 1:8 title "{/Symbol b}=1" with lines ls 1, \
     'r10unconf_beta1.out' using 1:5 notitle  with lines ls 2, \
     'r10unconf_beta2.out' using 1:8 title '{/Symbol b}=2' with lines ls 3, \
     'r10unconf_beta2.out' using 1:5 notitle  with lines ls 4, \
     'r10unconf_beta10.out' using 1:8 title '{/Symbol b}=10' with lines ls 7, \
     'r10unconf_beta10.out' using 1:5 notitle  with lines ls 8, \
     'r10unconf_beta100.out' using 1:8 title '{/Symbol b}=100' with lines ls 9, \
     'r10unconf_beta100.out' using 1:5 notitle  with lines ls 10, \
     'r10unconf_beta1000.out' using 1:8 title '{/Symbol b}=1000' with lines ls 11, \
     'r10unconf_beta1000.out' using 1:5 notitle  with lines ls 12, \
     'r10unconf_beta2.out' using 1:9  title 'E(S_s)' with lines ls 5, \
     'r10unconf_beta2.out' using 1:10  title 'E(S_y)' with lines ls 5
