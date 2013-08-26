# general types
for t in  "contour" "hydrograph" "particle" ; do
    echo "pumping_well $t output"
    ./ltaem ./pumping_well_${t}.in \
        1> pumping_well_${t}.screen \
        2> pumping_well_${t}.err
done

# specific PDEs associated with hydrograph types
for t in  "dualporosity" "leaky1" "leaky2" "leaky3" ; do 
    echo "pumping_well ${t} hydrographs"
    ./ltaem ./pumping_well_${t}_hydrographs.in \
        1> pumping_well_${t}_hydrographs.screen \
        2> pumping_well_${t}_hydrographs.err
done
