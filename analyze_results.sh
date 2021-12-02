#############################################
### WARPT pipeline - results analysis
#############################################


#!/bin/bash

###################################### Help ################################################
Help()
{
   # Display Help
   echo
   echo "Workflow for Association of Receptor Pairs from TREK-seq"
   echo
   echo "Help:"
   echo "Syntax: analyze_results.sh [-t -s -a -p -r -d]"
   echo "options:"
   echo "   -h     Print this Help."
   echo "   -t     table (.tsv) of alignments created with warpt"
   echo "   -s     sample name"
   echo "   -a     RNA seq annotations full path and name"
   echo "   -p     cluster proportion threshold"
   echo "   -r     cluster ratio threshold"
   echo "   -d     working directory for WARPT analysis"
   echo
}

################################ Default variables ########################################

MyTable="false"
MySample="mysample"
scRNAannotation="false"
MinProportion=0.5
MinRatio=2
BaseFolder=$(pwd)

################################ Options ###########################################

while getopts ":ht:s:a:p:r:d:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      t) 
         MyTable=$OPTARG;;
      s) 
         MySample=$OPTARG;;
      a) 
         scRNAannotation=$OPTARG;;
      p) 
         MinProportion=$OPTARG;;
      r) 
         MinRatio=$OPTARG;;
      d) 
         BaseFolder=$OPTARG;;
      :)    # If expected argument omitted:
         echo "Error: -${OPTARG} requires an argument."
         exit;;
      \?) # Invalid option
         echo "Error: Invalid option $OPTARG. Try analize_results.sh -h"
         exit;;
   esac
done


# create folders
if [ ! -e ${BaseFolder}/results/ ]; then
	mkdir ${BaseFolder}/results/
	mkdir ${BaseFolder}/results/plots/
fi


## -----------------------------------------------------------------------
### Results analysis

cp ${MyTable} ${BaseFolder}/results/
cp ${BaseFolder}/warpt/stats.log ${BaseFolder}/results/
cp ${BaseFolder}/warpt/QC/BC_UMI_cluster_metrics.txt ${BaseFolder}/results/

# with RNAseq annotations
if [ -f $scRNAannotation ]; then
	cp ${scRNAannotation} ${BaseFolder}/results/
   sum_up_results.r ${BaseFolder} ${MyTable} ${MySample} ${MinProportion} ${MinRatio} ${scRNAannotation}
   plots.r ${BaseFolder} ${MySample} ${scRNAannotation}

fi

# with NO RNAseq annotations
if [ ! -f $scRNAannotation ]; then
   sum_up_results.r ${BaseFolder} ${MyTable} ${MySample} ${MinProportion} ${MinRatio}
fi





