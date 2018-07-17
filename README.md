# SG10K: GRCh38 GATK4-gVCF Freebayes-VCF

This workflow processes SG10K samples from FastQ to CRAM and
creates a GATK4 gVCF as well as a filtered Freebayes VCF.

## Software dependencies

- Nextflow (min version 0.29)
- Singularity (min. version 2.5.1)

## How to run

A typical invocation would be 

```
nohup nextflow main.nf -c nextflow.config -params-file your-params.yaml -profile nscc --publishdir results -resume &
```

Only the `nscc` profile is maintained. Take special note, of the
singularity module dependency there and the singularity container
path.

## Parameters

Parameters listed under `params` in `nextflow.config` can be changed
on the command-line (see Nextflow manual)


Most if not all parameters should go into a file called `params.yaml`
and passed to Nextflow with `-params-file params.yaml`.  To create a
this file, concatenate `references.yaml` with the usual (GIS)
`samples.yaml`. The latter has to be of the new format, where
`readunits` are listed under `samples`. The GIS RPD snakemake framework has a
[converter script](https://github.com/gis-rpd/pipelines/blob/devel/tools/sample_conf_converter.py).

Examples can be found in the
`tests` directory.

## Output

**FIXME**
