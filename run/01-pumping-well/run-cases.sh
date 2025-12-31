# general types
# "particle_track"  <<<< this was blowing up and producing gigabytes of output

EXE="./ltaem"
#export LD_LIBRARY_PATH=/home/klkuhlm/src/gcc/local/lib64
#EXE="valgrind --leak-check=full --show-leak-kinds=all ./ltaem"

for t in  "contour" "timeseries" ; do
    echo "pumping_well $t output"
    ${EXE} ./pumping_well_${t}.in \
        1> pumping_well_${t}.screen \
        2> pumping_well_${t}.err 
    if [ -s pumping_well_${t}.err ] ; then
        echo "********************************** WARNING: ${t} produced an error **********************************"
        tail pumping_well_${t}.err
    fi
done

ln -sf ../plot_contour.py .
python plot_contour.py pumping_well_contour.dat

ln -sf ../plot_timeseries.py .
python plot_timeseries.py pumping_well_timeseries.dat

# specific PDEs associated with hydrograph types
for t in  "dualporosity" "leaky1" "leaky2" "leaky3" "kazemi" "warrenroot" "stepwiseQ" "linearQ"; do 
    echo "pumping_well ${t} timeseries"
    ${EXE} ./pumping_well_${t}_timeseries.in \
        1> pumping_well_${t}_timeseries.screen \
        2> pumping_well_${t}_timeseries.err
    if [ -s pumping_well_${t}_timeseries.err ] ; then
        echo "********************************** WARNING: ${t} produced an error **********************************"
        tail pumping_well_${t}_timeseries.err
    else
        python plot_timeseries.py pumping_well_${t}_timeseries.dat
    fi
done
