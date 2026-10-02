#!/bin/bash
# Label both PF versions and make the comparison PDF for each sample and timing setting.
set -e

BASE=/eos/user/c/chtong/Public/Rereco
SAMPLES="LLP_MH350_MS160_CTau10000 LLP_MH125_MS50_CTau3000"
TIMINGS="2 4"

for sample in $SAMPLES; do
    for ns in $TIMINGS; do
        dir=$BASE/$sample/${ns}ns_N2000
        # echo "===== $sample, ${ns} ns ====="

        # for pf in standardPF seedTimingPF; do
        #     python3 label_llp_clusters.py \
        #         --input  $dir/pfObjectsNtuple_${pf}.root \
        #         --output $dir/labels_${pf}.root
        # done

        python3 compare_timing_vs_standard.py \
            --standard $dir/labels_standardPF.root \
            --timing   $dir/labels_seedTimingPF.root \
            --timingLabel "seed timing PF (${ns} ns)" \
            --output   $dir/timing_vs_standard_LLP_${ns}ns_${sample}.pdf
    done
done