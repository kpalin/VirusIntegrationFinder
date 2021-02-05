from Bio import SeqIO,bgzf
import gzip
import pandas as pd

fastq_fname = str(snakemake.input.fastq)
inserts_fname = str(snakemake.input.insert_paf)
insert_name = str(snakemake.wildcards["insert_name"])

outfile_name = str(snakemake.output[0])
slop_distance = snakemake.params.SLOP


inserts = pd.read_table(inserts_fname, sep="\t", header=None)
inserts = inserts[[0, 1, 2, 3, 4, 5, 6, 7, 8,9]]
inserts.columns = [
    "read_id",
    "read_len",
    "start",
    "end",
    "strand",
    "insert_name",
    "insert_len",
    "insert_start",
    "insert_end",
    "matches"
]
inserts = inserts.set_index(["read_id","insert_name"]).sort_index()


from snakemake.logging import logger


with bgzf.BgzfWriter(outfile_name) as output_handle:
    for record in SeqIO.parse(gzip.open(fastq_fname, "rt"), "fastq"):
        insert_pos = inserts.loc[(record.id,insert_name)]
        if isinstance(insert_pos,pd.DataFrame):
            insert_pos = insert_pos.iloc[insert_pos.matches.argmax()]
        start = max(0,insert_pos.start - slop_distance)
        end = min(insert_pos.read_len, insert_pos.end + slop_distance)

        #record.seq = record.seq[start:end]
        record = record[start:end]
        #record.annotations["chunk"] = f"{start}:{end}"
        record.description += f" chunk={start}:{end}"
        #logger.info(record.description)
        #logger.info(str(insert_pos))
        SeqIO.write(record, handle=output_handle, format="fasta")