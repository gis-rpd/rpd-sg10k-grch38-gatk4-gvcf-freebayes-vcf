#!/bin/bash
which nextflow >/dev/null || exit 1
tmpdir=$(mktemp -d run-$(date +%Y%m%d-%H%M)-XXXXXX)
echo "Running in newly created $tmpdir"
cd $tmpdir
nohup nextflow ../../main.nf -c ../../nextflow.config -params-file ../../tests/params.yaml -profile nscc --publishdir results -resume &

