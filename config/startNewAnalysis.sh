#!/bin/bash

bash removeOutput.sh
#sframe_batch.py BoostedSingleTopAnalysisEle.xml
sframe_batch.py BoostedSingleTopAnalysisMuo.xml
#sframe_batch.py -s BoostedSingleTopAnalysisEle.xml
sframe_batch.py -s BoostedSingleTopAnalysisMuo.xml
#sframe_batch.py -a BoostedSingleTopAnalysisEle.xml
sframe_batch.py -a BoostedSingleTopAnalysisMuo.xml

echo "Done."

#htop -u matthies
