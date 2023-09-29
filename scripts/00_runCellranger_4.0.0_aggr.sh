#!/bin/sh
#$ -N aggr
#$ -S /bin/sh
#$ -cwd
#$ -q ngsc,ngsclargememory
#$ -l mem_free=55G,h_vmem=60G
#$ -pe ngsc-pe 8

. /etc/profile.d/modules.sh

module load sharedapps 

export PATH=/data/ngsc2data/NGSC_scripts/cellranger/cellranger-4.0.0:$PATH

cellranger aggr --id=aggr --csv=NZ01JHU501_libraries.csv --normalize=mapped --jobmode local --localcores=$NSLOTS --localmem=55
