#!/bin/bash

read -p "Do you really want do delete workdirs and output data? (y/n) " answer

if [ "$answer" == "n" ]
    then
    exit 1
elif [ "$answer" == "y" ]
    then
    #rm -r /nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Analysis/Electron/workdir_analysis_ele/
    rm -r /nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/PreSelection/Muon/workdir_PreSel_muo/
    #rm /nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Analysis/Electron/uhh2.AnalysisModuleRunner.*
    #rm /nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Analysis/Muon/uhh2.AnalysisModuleRunner.*
    #rm -r /nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/config/workdir_analysis*
else
    echo "No valid answer. Did nothing."
    exit 1
fi

echo "Done."