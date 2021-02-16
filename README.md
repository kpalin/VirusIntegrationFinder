# Snakemake workflow: VirusIntegrationFinder

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/virusintfind.svg?branch=master)](https://travis-ci.org/snakemake-workflows/virusintfind)

This is the Snakemake workflow for VirusIntegrationFinder. The main configuration in the `config.yaml` file.

## Authors

* Kimmo Palin (@kpalin)

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL
of this (original) repository and, if available, its DOI (see above).

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a
   template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to
   your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing `config/config.yml` file.  The important
parts are `INSERTED_SEQUENCE_FASTA` containing the inserted viral sequence, `INPUT_FASTQ` containing
the long whole genome reads and `GENOME_REFERENCE` for getting genomic context and coordinates of
the insertions.

### Step 3: Install conda execution environment

Install Snakemake using
[conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):
 
    conda create  -n vfind --file workflow/envs/myenv.yaml

For installation details, see the [instructions in the Snakemake
documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate vfind

Test your configuration by performing a dry-run via

    snakemake -s workflow/Snakefile -n

Execute the workflow locally via

    snakemake -s workflow/Snakefile--cores $N

using `$N` cores or run it in a cluster environment via

    snakemake -s workflow/Snakefile --cluster qsub --jobs 100

or

    snakemake -s workflow/Snakefile  --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for
further details.

### Step 5: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results
via:

    conda install jinja2 pygraphviz pygments
    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators.  *THIS IS NOT VERY USEFUL. JUST TECHNICAL
RUNTIME ETC*
### Step 6: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the
repository:

    git commit -a
    git push
