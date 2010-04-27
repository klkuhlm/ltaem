set xlabel 'index of singular value, j'
set ylabel 's_j/s_{1}'
set key outside below
set xrange [1:74]
set title 'relative singular values from SVD least-squares solution of elliptical problem'
set output 'mf_singular_values_allq.eps'
set terminal postscript eps enhanced color
plot 'ev.dat' using 1 title 'q_1: 2->1' with lines, \
     'ev.dat' using 3 title 'q_2: 2->1' with lines, \
     'ev.dat' using 5 title 'q_3: 2->1' with lines, \
     'ev.dat' using 7 title 'q_4: 2->1' with lines, \
     'ev.dat' using 9 title 'q_5: 2->1' with lines, \
     'ev.dat' using 11 title 'q_6: 2->1' with lines, \
     'ev.dat' using 13 title 'q_7: 2->1' with lines, \
     'ev.dat' using 15 title 'q_8: 2->1' with lines, \
     'ev.dat' using 17 title 'q_9: 2->1' with lines, \
     'ev.dat' using 19 title 'q_{10}: 2->1' with lines, \
     'ev.dat' using 21 title 'q_{11}: 2->1' with lines, \
     'ev.dat' using 23 title 'q_{12}: 2->1' with lines, \
     'ev.dat' using 25 title 'q_{13}: 2->1' with lines, \
     'ev.dat' using 27 title 'q_{14}: 2->1' with lines, \
     'ev.dat' using 29 title 'q_{15}: 2->1' with lines, \
     'ev.dat' using 31 title 'q_{16}: 2->1' with lines, \
     'ev.dat' using 33 title 'q_{17}: 2->1' with lines
set output 'mf_singular_values_allq2.eps'
plot 'ev.dat' using 2 title 'q_1: 1->2' with lines, \
     'ev.dat' using 4 title 'q_2: 1->2' with lines, \
     'ev.dat' using 6 title 'q_3: 1->2' with lines, \
     'ev.dat' using 8 title 'q_4: 1->2' with lines, \
     'ev.dat' using 10 title 'q_5: 1->2' with lines, \
     'ev.dat' using 12 title 'q_6: 1->2' with lines, \
     'ev.dat' using 14 title 'q_7: 1->2' with lines, \
     'ev.dat' using 16 title 'q_8: 1->2' with lines, \
     'ev.dat' using 18 title 'q_9: 1->2' with lines, \
     'ev.dat' using 20 title 'q_{10}: 1->2' with lines, \
     'ev.dat' using 22 title 'q_{11}: 1->2' with lines, \
     'ev.dat' using 24 title 'q_{12}: 1->2' with lines, \
     'ev.dat' using 26 title 'q_{13}: 1->2' with lines, \
     'ev.dat' using 28 title 'q_{14}: 1->2' with lines, \
     'ev.dat' using 30 title 'q_{15}: 1->2' with lines, \
     'ev.dat' using 32 title 'q_{16}: 1->2' with lines, \
     'ev.dat' using 34 title 'q_{17}: 1->2' with lines
