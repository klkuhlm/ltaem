#!/bin/bash

cd 01-pumping-well
ln -sf ../ltaem .
./run-cases.sh
cd ..

for d in {01a-,02-,03-,04-,05-,06-,07-,08-,09-,10-,11-}; do
    cd ${d}*
    echo "now in $(pwd)"
    ln -sf ../ltaem .
    for i in *.in; do
        echo "running ${i}"
        ./ltaem ${i} > ${i}_screen.out
    done
    cd ..
done

            
