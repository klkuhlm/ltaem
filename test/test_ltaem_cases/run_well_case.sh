#!/bin/bash
set -e  # terminate on failed step
export OMP_NUM_THREADS=6

exe=../../run/ltaem

# one well at x=0.55,y=0.55

# run contours test case
${exe} well_contours.in #2>&1 | tee well_contours.screen
gnuplot plot_well_contours.gnu

# run hydrograph test case
${exe} well_hydrograph.in #2>&1 | tee well_hydrograph.screen
gnuplot plot_well_hydrograph.gnu

# run particle test case
