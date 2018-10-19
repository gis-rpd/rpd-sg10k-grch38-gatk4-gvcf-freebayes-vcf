#!/bin/bash

# check if nextflow is available
which nextflow >/dev/null || exit 1
# make sure we are in the test directory
ls ../main.nf ../tests  >& /dev/null || exit 1

PARAMS=params.yaml

if [ -z "$WORKDIR" ]; then
    echo "ERROR: export WORKDIR firsti (e.g /seq/astar/gis/rpd/testing/output/rpd-sg10k-grch38-gatk4-gvcf-freebayes)" 1>&2
    exit 1
fi 

suffix=run-$(date +%Y%m%d-%H%M)
tmpdir=$(mktemp --tmpdir -d ${suffix}-XXXXXX) || exit  1
echo "INFO: Cloning setup with rsync to $tmpdir"
rsync -q -av --exclude tests --exclude .git  ../ $tmpdir/
cp $PARAMS $tmpdir
echo "INFO: This copy makes messing around easier; just don't forget to incorporate changes meant to become permanent"
cd $tmpdir

workdir=$WORKDIR/$suffix/work
publishdir=$WORKDIR/$suffix/results
echo "INFO: Running in newly created $tmpdir. workdir is $workdir. publishdir is $publishdir"
nohup nextflow main.nf -c nextflow.config -params-file $(basename $PARAMS) -w $workdir -profile nscc --publishdir $publishdir --keep_workdir -resume &

