#!/bin/bash

# check if nextflow is available
which nextflow >/dev/null || exit 1
# make sure we are in the test directory
ls ../main.nf ../tests  >& /dev/null || exit 1

PARAMS=params.yaml

basedir="$1"
if [ -z "$basedir" ]; then
    echo "ERROR: Missing basedir argument. This is where results and workdir will be created." 1>&2
    echo "You can for example use: basedir=\$(mktemp -d /data/users/astar/gis/rpd/testing/output/rpd-sg10k-grch38-gatk4-gvcf-freebayes/\$(date -I)-XXXXXX)" 1>&2
    exit 1
fi

suffix=run-$(date +%Y%m%d-%H%M)
tmpdir=$(mktemp --tmpdir -d ${suffix}-XXXXXX) || exit  1
echo "INFO: Cloning setup with rsync to $tmpdir"
rsync -q -av --exclude tests --exclude .git  ../ $tmpdir/
cp $PARAMS $tmpdir
echo "INFO: This copy makes messing around easier; just don't forget to incorporate changes meant to become permanent"
cd $tmpdir

workdir=$basedir/$suffix/work
publishdir=$basedir/$suffix/results
echo "INFO: Running in newly created $tmpdir. workdir is $workdir. publishdir is $publishdir"
nohup nextflow main.nf -c nextflow.config -params-file $(basename $PARAMS) -w $workdir -profile nscc --publishdir $publishdir --keep_workdir -resume &

