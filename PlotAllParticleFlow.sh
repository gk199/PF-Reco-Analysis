#!/bin/bash

for label in standardPF cellTimingPF seedTimingPF depth1SeedTimingPF; do
    python3 Plotting/plot_pfObjects.py \
        --input  pfObjectsNtuple_${label}.root \
        --output pfObjectsHistos_${label}.root
done

python3 Plotting/compare_hcal_clusters.py --inputdir . --output hcal_comparison.root --pdf hcal_comparison.pdf