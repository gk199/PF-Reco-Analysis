#!/bin/bash

echo " "
echo "Running all ParticleFlow tests, on 100 events MC sample"
echo " "
echo "-----------------------------------"
echo "Standard Particle Flow:"
echo "-----------------------------------"
echo " "
cd ..
cmsenv
cp PF-Reco-Analysis/PFTestingAlgos/Basic2DGenericTopoClusterizer_original.cc.edit RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericTopoClusterizer.cc 
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHBHE_original_cfi.py RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHBHE_cfi.py 
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_original.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_original.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc 
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_original_cfi.py RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py 
scram b -j 8
cd PF-Reco-Analysis
cmsRun MyPFStudy_ReReco_MC_Sim_DIGI_RAW2DIGI_L1Reco_RECO.py
mv pf_only_reReco_MC_Sim.root pf_only_reReco_MC_Sim_standardPF.root

echo " "
echo "-----------------------------------"
echo "Modified Particle Flow: cell level timing cut vs global highest energy seed"
echo "-----------------------------------"
echo " "
cd ..
cp PF-Reco-Analysis/PFTestingAlgos/Basic2DGenericTopoClusterizer_timing.cc.edit RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericTopoClusterizer.cc 
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHBHE_timing_cfi.py RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHBHE_cfi.py 
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_timing.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_timing.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc 
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_timing_cfi.py RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py 
scram b -j 8
cd PF-Reco-Analysis
cmsRun MyPFStudy_ReReco_MC_Sim_DIGI_RAW2DIGI_L1Reco_RECO.py
mv pf_only_reReco_MC_Sim.root pf_only_reReco_MC_Sim_cellTimingPF.root

echo " "
echo "-----------------------------------"
echo "Modified Particle Flow: seed level timing cut vs global highest energy seed"
echo "-----------------------------------"
echo " "
cd ..
cp PF-Reco-Analysis/PFTestingAlgos/Basic2DGenericTopoClusterizer_original.cc.edit RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericTopoClusterizer.cc 
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHBHE_original_cfi.py RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHBHE_cfi.py 
echo "Note: only timing cut on seeds, not on all cells, so only the clusterizer and cluster producer are modified"
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_timing.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_seedTiming.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc 
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_seedTiming_cfi.py RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py 
scram b -j 8
cd PF-Reco-Analysis
cmsRun MyPFStudy_ReReco_MC_Sim_DIGI_RAW2DIGI_L1Reco_RECO.py
mv pf_only_reReco_MC_Sim.root pf_only_reReco_MC_Sim_seedTimingPF.root

echo " "
echo "-----------------------------------"
echo "Modified Particle Flow: seed level timing cut vs first cluster seed"
echo "-----------------------------------"
echo " "
cd ..
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_timing.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_depth1Timing.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc 
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_depth1Timing_cfi.py RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py 
scram b -j 8
cd PF-Reco-Analysis
cmsRun MyPFStudy_ReReco_MC_Sim_DIGI_RAW2DIGI_L1Reco_RECO.py
mv pf_only_reReco_MC_Sim.root pf_only_reReco_MC_Sim_depth1SeedTimingPF.root

echo " "
echo "All tests completed, output files are in PF-Reco-Analysis directory"
echo " "
echo "-----------------------------------"
echo "Reverting to standard Particle Flow configuration"
echo "-----------------------------------"
echo " "
cd ..
cp PF-Reco-Analysis/PFTestingAlgos/Basic2DGenericTopoClusterizer_original.cc.edit RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericTopoClusterizer.cc 
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHBHE_original_cfi.py RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHBHE_cfi.py 
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_original.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_original.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc 
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_original_cfi.py RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py 
scram b -j 8
echo "Reverted to standard Particle Flow configuration"
echo "Done!"