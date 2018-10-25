#!/bin/bash

#bash removeOutput.sh
#sframe_batch.py BoostedSingleTopAnalysisEle.xml
#sframe_batch.py BoostedSingleTopAnalysisMuo.xml
#sframe_batch.py -s BoostedSingleTopAnalysisEle.xml
#sframe_batch.py -s BoostedSingleTopAnalysisMuo.xml
#sframe_batch.py -a BoostedSingleTopAnalysisEle.xml
#sframe_batch.py -a BoostedSingleTopAnalysisMuo.xml

#sframe_batch.py BoostedSingleTopPreSelectionEle-SYS.xml
sframe_batch.py -a BoostedSingleTopPreSelectionEle-SYS.xml
#sframe_batch.py BoostedSingleTopPreSelectionEle-SYSdown.xml
sframe_batch.py -a BoostedSingleTopPreSelectionEle-SYSdown.xml
#sframe_batch.py BoostedSingleTopPreSelectionMuo-SYS.xml
sframe_batch.py -a BoostedSingleTopPreSelectionMuo-SYS.xml
#sframe_batch.py BoostedSingleTopPreSelectionMuo-SYSdown.xml
sframe_batch.py -a BoostedSingleTopPreSelectionMuo-SYSdown.xml

echo "Done."

#htop -u matthies
