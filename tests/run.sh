#!/bin/bash

# check if nextflow is available
which nextflow >/dev/null || exit 1
# make sure we are in the test directory
ls ../main.nf ../tests  >& /dev/null || exit 1

tmpdir=$(mktemp -d run-$(date +%Y%m%d-%H%M)-XXXXXX)
echo "INFO: Cloning setup wth rsync to $tmpdir"
rsync -q -av --exclude tests --exclude .git  ../ $tmpdir/
echo "INFO: This copy makes messing around easier; just don't forget to incorporate changes meant to become permanent"
echo "INFO: Running in newly created $tmpdir"
cd $tmpdir

nohup nextflow main.nf  -c nextflow.config -params-file params.yaml -profile nscc --publishdir results -resume &

