# general types
for t in  "contour" "hydrographs" "particle_track" ; do
    echo "pumping_well $t output"
    ./ltaem ./pumping_well_${t}.in \
        1> pumping_well_${t}.screen \
        2> pumping_well_${t}.err 
    if [ -s pumping_well_${t}.err ] ; then
        echo "**********************************"
        echo "WARNING: ${t} produced an error"
        echo "**********************************"
        head pumping_well_${t}.err
    fi
done

# specific PDEs associated with hydrograph types
for t in  "dualporosity" "leaky1" "leaky2" "leaky3" ; do 
    echo "pumping_well ${t} hydrographs"
    ./ltaem ./pumping_well_${t}_hydrographs.in \
        1> pumping_well_${t}_hydrographs.screen \
        2> pumping_well_${t}_hydrographs.err
    if [ -s pumping_well_${t}_hydrographs.err ] ; then
        echo "**********************************"
        echo "WARNING: ${t} produced an error"
        echo "**********************************"
        head pumping_well_${t}_hydrographs.err
    fi
done