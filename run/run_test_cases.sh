#!/bin/bash

EXE="./ltaem"
#EXE="valgrind --leak-check=full --show-leak-kinds=all ./ltaem"

cd 01-pumping-well
echo "now in $(pwd)"
ln -sf ../ltaem .
./run-cases.sh
cd ..

#export LD_LIBRARY_PATH=/usr/local/stow/gcc-trunk/lib64  

# 06-

for d in {01a-,02-,03-,04-,06-,07-,08-,09-,10-,11-,12-}; do
    cd ${d}*
    echo "now in $(pwd)"
    ln -sf ../ltaem .
    for i in *.in; do
        echo "running ${i}"
        ${EXE} ${i} 1> ${i}_screen.out 2> ${i}_screen.err
        if [ -s ${i}_screen.err ] ; then
            echo "********************************** WARNING: ${t} produced an error **********************************"
            tail ${i}_screen.err
        fi            
    done
    cd ..
done

cd 05-particle-sensitivity
echo "now in $(pwd)"
ln -sf ../ltaem .
python particles-input.py > python_particle_input_screen.out
python plot-results.py > python_plot_results_screen.out
cd ..

