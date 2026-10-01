# 1. Set up the environment (ROOT and python3)
cd /afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src
cmsenv
cd PF-Reco-Analysis/PF_Timing_Evaluation

# 2. Label the clusters in each PF version
D=/eos/user/c/chtong/Public/Rereco/LLP/20260930T081201Z_19373
D=/eos/user/c/chtong/Public/Rereco/LLP_MH350_MS160_CTau10000/20260930T230839Z_427
python3 label_llp_clusters.py --input $D/pfObjectsNtuple_standardPF.root   --output labels_standardPF.root
python3 label_llp_clusters.py --input $D/pfObjectsNtuple_seedTimingPF.root --output labels_seedTimingPF.root

# 3. Make the comparison plots
python3 compare_timing_vs_standard.py \
    --standard labels_standardPF.root \
    --timing   labels_seedTimingPF.root \
    --output   timing_vs_standard_LLP.pdf