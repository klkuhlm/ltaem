set logscale x
set xlabel 'time since pumping began (min)'
set ylabel 'drawdown (ft)'
set key outside below
set terminal postscript eps enhanced color
set size 0.75,0.75
set xrange [0.05:400]
set title 'Pumping well C6 and injection well C5'
set output 'results_pumping.eps'
plot 'boise06' index 11 using 1:2 notitle with l, \
     'test1_both.dat' using 1:14 title 'C5' with p, \
     'boise06' index 12 using 1:2 notitle with l, \
     'test1_both.dat' using 1:16 title 'C6' with p
#     
set output 'obs01_results.eps'
set title 'observation wells A1-B3'
plot 'boise06' index 0 using 1:2 notitle with l, \
     'test1_both.dat' using 1:10 title 'A1' with p, \
     'boise06' index 1 using 1:2 notitle with l, \
     'test1_both.dat' using 1:9 title 'B1' with p, \
     'boise06' index 2 using 1:2 notitle with l, \
     'test1_both.dat' using 1:15 title 'B2' with p, \
     'boise06' index 3 using 1:2 notitle with l, \
     'test1_both.dat' using 1:11 title 'B3' with p
##     
set output 'obs02_results.eps'
set title 'observation wells B4-C1'
plot 'boise06' index 4 using 1:2 notitle with l, \
     'test1_both.dat' using 1:7 title 'B4' with p, \
     'boise06' index 5 using 1:2 notitle with l, \
     'test1_both.dat' using 1:12 title 'B5' with p, \
     'boise06' index 6 using 1:2 notitle with l, \
     'test1_both.dat' using 1:5 title 'B6' with p, \
     'boise06' index 7 using 1:2 notitle with l, \
     'test1_both.dat' using 1:6 title 'C1' with p
##     
set output 'obs03_results.eps'
set title 'observations wells C2-C4, X2 & X4'
plot 'boise06' index 8 using 1:2 notitle with l, \
     'test1_both.dat' using 1:2 title 'C2' with p, \
     'boise06' index 9  using 1:2 notitle with l, \
     'test1_both.dat' using 1:13 title 'C3' with p, \
     'boise06' index 10 using 1:2 notitle with l, \
     'test1_both.dat' using 1:8 title 'C4' with p, \
     'boise06' index 14 using 1:2 notitle with l, \
     'test1_both.dat' using 1:4 title 'X2' with p, \
     'boise06' index 16 using 1:2 notitle with l, \
     'test1_both.dat' using 1:3 title 'X4' with p     
r = 248.00
set xrange [0.3:100]
set output 'pumping_recov_results.eps'
plot 'boise06' index 11 using ($1-r):2 notitle with l, \
     'test1_both.dat' using ($1-r):14 title 'C5' with p, \
     'boise06' index 12 using ($1-r):2 notitle with l, \
     'test1_both.dat' using ($1-r):16 title 'C6' with p
##     
set output 'obs01_recov_results.eps'
set title 'observation wells A1-B3'
plot 'boise06' index 0 using ($1-r):2 notitle with l, \
     'test1_both.dat' using ($1-r):10 title 'A1' with p, \
     'boise06' index 1 using ($1-r):2 notitle with l, \
     'test1_both.dat' using ($1-r):9 title 'B1' with p, \
     'boise06' index 2 using ($1-r):2 notitle with l, \
     'test1_both.dat' using ($1-r):15 title 'B2' with p, \
     'boise06' index 3 using ($1-r):2 notitle with l, \
     'test1_both.dat' using ($1-r):11 title 'B3' with p
##     
set output 'obs02_recov_results.eps'
set title 'observation wells B4-C1'
plot 'boise06' index 4 using ($1-r):2 notitle with l, \
     'test1_both.dat' using ($1-r):7 title 'B4' with p, \
     'boise06' index 5 using ($1-r):2 notitle with l, \
     'test1_both.dat' using ($1-r):12 title 'B5' with p, \
     'boise06' index 6 using ($1-r):2 notitle with l, \
     'test1_both.dat' using ($1-r):5 title 'B6' with p, \
     'boise06' index 7 using ($1-r):2 notitle with l, \
     'test1_both.dat' using ($1-r):6 title 'C1' with p
##     
set output 'obs03_recov_results.eps'
set title 'observations wells C3,C4,C6,X2 & X4'
plot 'boise06' index 8 using ($1-r):2 notitle with l, \
     'test1_both.dat' using ($1-r):2 title 'C2' with p, \
     'boise06' index 9 using ($1-r):2 notitle with l, \
     'test1_both.dat' using ($1-r):13 title 'C3' with p, \
     'boise06' index 10 using ($1-r):2 notitle with l, \
     'test1_both.dat' using ($1-r):8 title 'C4' with p, \
     'boise06' index 14 using ($1-r):2 notitle with l, \
     'test1_both.dat' using ($1-r):4 title 'X2' with p, \
     'boise06' index 16 using ($1-r):2 notitle with l, \
     'test1_both.dat' using ($1-r):3 title 'X4' with p          
