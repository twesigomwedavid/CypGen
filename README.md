# CypGen
Calling human cytochrome P450 star alleles by leveraging genome graph-based variant detection.

Model gene: CYP2D6

## Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow management system that facilitates parallelisation, scalability, reproducibility and portability of pipelines via [`Docker`](https://docs.docker.com) and [`Singularity`](https://sylabs.io/) technologies.

## Quick Start

The following are required to run the CypGen pipeline;

1. Software
    - [`Nextflow`](https://nf-co.re/usage/installation)
    - [`Singularity`](https://sylabs.io/) (especially for HPC environments running Linux OS) or [`Docker`](https://docs.docker.com) (recommended for MacOS users)

2. Whole genome sequence (WGS) data
    - Indexed BAM/CRAM files
    
3. Reference genome
    - hg19, b37, or hg38
    
Note: For a full description of the differences among reference genomes, please check out this excellent [`Documentation`] (https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies). For the purpose of using this pipeline, if the GRCh37 reference genome you are using has contigs that start with 'chr' (i.e. chr1, chr2, ..., chrX, chrM, ...), use the hg19 option. You should use the b37 option if the contigs in the reference genome do not have 'chr' (i.e. 1, 2, ..., X, MT). For GRCh38, the hg38 option is sufficient.

### Downloading the pipeline
```bash
git clone https://github.com/twesigomwedavid/CypGen.git
```

