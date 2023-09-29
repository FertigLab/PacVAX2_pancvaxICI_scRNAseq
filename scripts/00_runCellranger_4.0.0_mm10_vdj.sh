#!/bin/sh
#$ -N vdj
#$ -S /bin/sh
#$ -cwd
#$ -q ngsc,ngsclargememory
#$ -l mem_free=110G,h_vmem=120G
#$ -pe ngsc-pe 8,16

. /etc/profile.d/modules.sh

module load sharedapps 

export PATH=/data/ngsc2data/NGSC_scripts/cellranger/cellranger-4.0.0:$PATH

samplename=`echo $in | sed -e s/^..._//`

echo $TMPDIR
cd $TMPDIR

cp -r /data/ngsc2data/${projectid}/${projectid}_${in} .

ls -l

cellranger vdj --id=${samplename} --sample=${samplename} --fastqs=${projectid}_${in} --reference=/data/ngsc2data/NGSC_scripts/cellranger/refdata-cellranger-vdj-GRCm38-alts-ensembl-4.0.0 --localcores=$NSLOTS --localmem=110

cp -r ${samplename}/outs /data/ngsc2data/${projectid}/${projectid}_000_analysis/cellranger/${samplename}
