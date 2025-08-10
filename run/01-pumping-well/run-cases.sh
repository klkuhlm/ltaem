# general types
# "particle_track"  <<<< this was blowing up and producing gigabytes of output

EXE="./ltaem"
# EXE="valgrind --leak-check=full --show-leak-kinds=all ./ltaem"

for t in  "contour" "hydrographs" ; do
    echo "pumping_well $t output"
    ${EXE} ./pumping_well_${t}.in \
        1> pumping_well_${t}.screen \
        2> pumping_well_${t}.err 
    if [ -s pumping_well_${t}.err ] ; then
        echo "**********************************"
        echo "WARNING: ${t} produced an error"
        echo "**********************************"
        cat pumping_well_${t}.err
    fi
done

# specific PDEs associated with hydrograph types
for t in  "dualporosity" "leaky1" "leaky2" "leaky3" "stepwiseQ" "linearQ"; do 
    echo "pumping_well ${t} hydrographs"
    ${EXE} ./pumping_well_${t}_hydrographs.in \
        1> pumping_well_${t}_hydrographs.screen \
        2> pumping_well_${t}_hydrographs.err
    if [ -s pumping_well_${t}_hydrographs.err ] ; then
        echo "**********************************"
        echo "WARNING: ${t} produced an error"
        echo "**********************************"
        cat pumping_well_${t}_hydrographs.err
    fi
done
