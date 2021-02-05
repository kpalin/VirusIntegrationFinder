from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")
INSERT_NAMES=[x[1:].strip().split()[0] for x in open(config["INSERTED_SEQUENCE_FASTA"]) if x.startswith(">")]
assert len(INSERT_NAMES)>0,"Couldn't find any insert names from "+ config["INSERTED_SEQUENCE_FASTA"]