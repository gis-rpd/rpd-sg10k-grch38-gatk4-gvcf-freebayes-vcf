# SG10K Health: GRCh38 GATK4-gVCF Freebayes-VCF

This is the main genome analytics workflow powering the production analysis of whole genome samples
for the Singapore National Precision Medicine (NPM) Program Phase 1A, sometimes also referred to SG10K
Health. It processes samples from FastQ to lossless CRAM, computes multiple QC metrics as well as Freebayes
variant calls and GATK4 gvcfs.

To ensure reproducibility, scalability and mobility the workflow is implemented as [Nextflow](https://www.nextflow.io/) recipe and uses containers
([Singularity](https://www.sylabs.io/docs/) on [NSCC's Aspire 1](https://www.nscc.sg/about-nscc/our-facilityaspire-1/) and [Docker](https://www.docker.com) on
[AWS Batch](https://aws.amazon.com/batch/)). Container building is simplified by the use of
[Bioconda](https://bioconda.github.io/).




## Output

All results can be found in the `results` folder of a pipeline
execution. Results there are grouped per sample, with the exception of
Goleft indexcov, which summarises over the sample set.

### Main results

- [GATK4](https://software.broadinstitute.org/gatk/gatk4) gVCF (indexed): `{sample}/{sample}.g.vcf.gz`
- [Freebayes](https://github.com/ekg/freebayes) VCF (Q>=20; indexed): `{sample}/{sample}.fb.vcf.gz`
- CRAM (lossless, with OQ, indexed): `{sample}/{sample}.bqsr.cram`

### QC etc.

- [Goleft](https://github.com/brentp/goleft) indexcov: `indexcov/all/` (main file `indexcov/all/all.html`)
- [Samtools](http://www.htslib.org/doc/samtools.html) stats: `{sample}/stats/` (main files: `{sample}/stats/{sample}.stats` and `{sample}/stats/{sample}.html`)
- [Verifybamid](https://genome.sph.umich.edu/wiki/VerifyBamID) for the three ethnicities: `{sample}/verifybamid/` (main files: `{sample}/verifybamid/{sample}.SGVP_MAF0.01.{ethnicity}.selfSM`)
- Coverage as per SOP: `{sample}/{sample}.cov-062017.txt`


## Notes

- GATK commandline parameters are based on [the official WDL implementation](https://github.com/broadinstitute/wdl/tree/develop/scripts/broad_pipelines/germline-short-variant-discovery/gvcf-generation-per-sample/1.0.0)
- This repo was migrated from Bitbucket, from where some issues/requests are not yet moved over
- We share this code for transparency. This is not meant to a generic whole genome workflow for wider use, but rather specific to the program's needs.
 For the same reason this documentation is rudimentary.

## Authors

The workflow was implemented in the [Genome Institute of Singapore
(GIS)](https://www.a-star.edu.sg/gis) by:

- Lavanya VEERAVALLI <veeravallil@gis.a-star.edu.sg>
- Andreas WILM <wilma@gis.a-star.edu.sg>




