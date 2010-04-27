##set style line 1 linetype 1 linecolor rgb "red" linewidth 4
##set style line 2 linetype 2 linecolor rgb "red" linewidth 0.5
##set style line 3 linetype 1 linecolor rgb "orange" linewidth 4
##set style line 4 linetype 2 linecolor rgb "orange" linewidth 0.5
##set style line 5 linetype 1 linecolor rgb "black" linewidth 4
##set style line 6 linetype 2 linecolor rgb "black" linewidth 0.5
##set style line 7 linetype 1 linecolor rgb "green" linewidth 4 
##set style line 8 linetype 2 linecolor rgb "green" linewidth 0.5
##set style line 9 linetype 1 linecolor rgb "cyan" linewidth 4
##set style line 10 linetype 2 linecolor rgb "cyan" linewidth 0.5
##set style line 11 linetype 1 linecolor rgb "blue" linewidth 4
##set style line 12 linetype 2 linecolor rgb "blue" linewidth 0.5
##set style line 13 linetype 1 linecolor rgb "violet" linewidth 4
##set style line 14 linetype 2 linecolor rgb "violet" linewidth 0.5
set style line 1 linetype 2 linewidth 2
set style line 2 linetype 2 linewidth 4
set style line 3 linetype 3 linewidth 2
set style line 4 linetype 2 linewidth 4
set style line 5 linetype 1 linewidth 2
set style line 6 linetype 2 linewidth 4
set style line 7 linetype 4 linewidth 2 
set style line 8 linetype 2 linewidth 4
set style line 9 linetype 5 linewidth 2
set style line 10 linetype 2 linewidth 4
set style line 11 linetype 6 linewidth 2
set style line 12 linetype 2 linewidth 4
set style line 13 linetype 1 linewidth 4
set style line 14 linetype 2 linewidth 4
set format y "10^{%T}"
set format x "10^{%T}"
set logscale x
set logscale y
unset title
set xlabel 't_D' 0.0,0.5
set ylabel 's_D' 2.5,0.0
#set autoscale x
#set autoscale y
set xrange [1.0E-2:1.0E+6]
set yrange [1.0E-1:50]
set label 1 '{/Symbol t}=0.1' at 20,0.4
set label 2 '{/Symbol t}=0.01' at 3,6
set label 3 '{/Symbol t}=10^{-3}' at 0.3,3
set label 4 '{/Symbol t}=10^{-4}' at 0.06,1
set label 5 '{/Symbol t}=10^{-5}' at 0.02,0.3
set label 6 'E_{1}(t_{D}/4)' at 1E+4,6
set output 'transient_sources.eps'
unset key
set size 0.75,0.75
##
## change bounding box command in eps file to 60 60 315 229
##
set terminal postscript eps enhanced mono 'Times-Roman' 22
plot 'trans_tau-5_dimless.out' using 1:6  with lines ls 11, \
     'trans_tau-4_dimless.out' using 1:6  with lines ls 9, \
     'trans_tau-3_dimless.out' using 1:6  with lines ls 7, \
     'trans_tau-2_dimless.out' using 1:6  with lines ls 3, \
     'trans_tau-1_dimless.out' using 1:6  with lines ls 1, \
     'trans_tau-1_dimless.out' using 1:2  with lines ls 5
set logscale x
set logscale y
set xrange [0.01:2]
set yrange [0.001:10]
set ylabel 's_D' 1.0, 0.0
set xlabel 'r'
set format x "10^{%T}"
set format y "10^{%T}"
unset label 1
unset label 2
unset label 3
unset label 4
unset label 5
unset label 6
set label 1 '{/Symbol t}=0.1' at 0.05,1.0
set label 2 '{/Symbol t}=0.01' at 0.18,0.05
set label 3 '{/Symbol t}=10^{-3}' at 0.8, 0.02
set label 4 'E_{1}(t_{D}/4)' at 0.1,4.0
set size 0.75,0.75
set arrow from 1.0,0.015 to 0.8,0.002 front
set output 'transient_sources_surface_earlyt.eps'
plot 'wave_surface_plot_59_1.0E-4.dat' using 1:3 with lines ls 9, \
     'wave_surface_plot_59_1.0E-3.dat' using 1:3 with lines ls 7, \
     'wave_surface_plot_59_1.0E-2.dat' using 1:3 with lines ls 3, \
     'wave_surface_plot_59_1.0E-1.dat' using 1:3 with lines ls 1, \
     'wave_surface_plot_59_1.0E-4.dat' using 1:2 with lines ls 5
#set autoscale
#set yrange [0.001:20]
#set output 'transient_sources_surface_midt.eps'
#plot 'wave_surface_plot_81_1.0E-4.dat' using 1:2 with lines ls 5, \
#     'wave_surface_plot_81_1.0E-4.dat' using 1:3 with lines ls 9, \
#     'wave_surface_plot_81_1.0E-3.dat' using 1:3 with lines ls 7, \
#     'wave_surface_plot_81_1.0E-2.dat' using 1:3 with lines ls 3, \
#     'wave_surface_plot_81_1.0E-1.dat' using 1:3 with lines ls 1
#set autoscale
#set yrange [0.001:30]
#set output 'transient_sources_surface_latet.eps'
#plot 'wave_surface_plot_103_1.0E-4.dat' using 1:2 with lines ls 5, \
#     'wave_surface_plot_103_1.0E-4.dat' using 1:3 with lines ls 9, \
#     'wave_surface_plot_103_1.0E-3.dat' using 1:3 with lines ls 7, \
#     'wave_surface_plot_103_1.0E-2.dat' using 1:3 with lines ls 3, \
#     'wave_surface_plot_103_1.0E-1.dat' using 1:3 with lines ls 1
    
