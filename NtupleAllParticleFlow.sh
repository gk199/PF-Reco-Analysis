#!/bin/bash

cmsenv

cmsRun PFObjectsNtupler/python/runPFObjectsNtupler_cfg.py \
    inputFiles=file:pf_only_reReco_MC_Sim_standardPF.root \
    outputFile=pfObjectsNtuple_standardPF.root

cmsRun PFObjectsNtupler/python/runPFObjectsNtupler_cfg.py \
    inputFiles=file:pf_only_reReco_MC_Sim_cellTimingPF.root \
    outputFile=pfObjectsNtuple_cellTimingPF.root

cmsRun PFObjectsNtupler/python/runPFObjectsNtupler_cfg.py \
    inputFiles=file:pf_only_reReco_MC_Sim_seedTimingPF.root \
    outputFile=pfObjectsNtuple_seedTimingPF.root

cmsRun PFObjectsNtupler/python/runPFObjectsNtupler_cfg.py \
    inputFiles=file:pf_only_reReco_MC_Sim_depth1SeedTimingPF.root \
    outputFile=pfObjectsNtuple_depth1SeedTimingPF.root
