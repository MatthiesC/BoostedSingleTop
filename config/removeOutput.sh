#!/bin/bash

read -p "Do you really want do delete ANALYSIS workdirs and output data? (Type e/m for electron/muon channel) " answer

if [ "$answer" == "m" ]
    then
    rm -r /nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Analysis/Muon/workdir_analysis_muo/
    rm /nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Analysis/Muon/uhh2.AnalysisModuleRunner.*
    rm -r /nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/config/workdir_analysis_muo
elif [ "$answer" == "e" ]
    then
    rm -r /nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Analysis/Electron/workdir_analysis_ele/
    rm /nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Analysis/Electron/uhh2.AnalysisModuleRunner.*
    rm -r /nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/config/workdir_analysis_ele
else
    echo "No valid answer. Did nothing."
    exit 1
fi

echo "Done."
