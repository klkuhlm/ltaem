set xrange [0:2*pi]
set logscale y
set output 'line_compare_flux_22.eps'
set terminal postscript eps enhanced color linewidth 2.0
plot 'line_compare_flux_4_22.out.plot' using 1:(abs($2-$3)) title 'N=4' with lines ,\
     'line_compare_flux_8_22.out.plot' using 1:(abs($2-$3)) title 'N=8' with lines ,\
     'line_compare_flux_12_22.out.plot' using 1:(abs($2-$3)) title 'N=12' with lines ,\
     'line_compare_flux_16_22.out.plot' using 1:(abs($2-$3)) title 'N=16' with lines
##
unset logscale y
set output 'line_compare_flux_22a.eps'     
plot 'line_compare_flux_24_42.out.plot' using 1:2 title 'N=4' with lines ,\
     'line_compare_flux_24_42.out.plot' using 1:3 notitle 'N=4' with lines
##
set logscale y
set output 'line_compare_flux_32.eps'
plot 'line_compare_flux_4_32.out.plot' using 1:(abs($2-$3)) title 'N=4' with lines ,\
     'line_compare_flux_8_32.out.plot' using 1:(abs($2-$3)) title 'N=8' with lines ,\
     'line_compare_flux_12_32.out.plot' using 1:(abs($2-$3)) title 'N=12' with lines ,\
     'line_compare_flux_16_32.out.plot' using 1:(abs($2-$3)) title 'N=16' with lines, \
     'line_compare_flux_20_32.out.plot' using 1:(abs($2-$3)) title 'N=20' with lines, \
     'line_compare_flux_24_32.out.plot' using 1:(abs($2-$3)) title 'N=24' with lines
##
set output 'line_compare_flux_42.eps'
plot 'line_compare_flux_4_42.out.plot' using 1:(abs($2-$3)) title 'N=4' with lines ,\
     'line_compare_flux_8_42.out.plot' using 1:(abs($2-$3)) title 'N=8' with lines ,\
     'line_compare_flux_12_42.out.plot' using 1:(abs($2-$3)) title 'N=12' with lines ,\
     'line_compare_flux_16_42.out.plot' using 1:(abs($2-$3)) title 'N=16' with lines, \
     'line_compare_flux_20_42.out.plot' using 1:(abs($2-$3)) title 'N=20' with lines, \
     'line_compare_flux_24_42.out.plot' using 1:(abs($2-$3)) title 'N=24' with lines          

