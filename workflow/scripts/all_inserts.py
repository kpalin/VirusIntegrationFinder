from Bio import SeqIO,bgzf
import gzip
import re


infile_names = snakemake.input
outfile_name = snakemake.output[0]

file_pat = re.compile(r"results/inserts/([^ /]+)/consensus_([ ]*)[.]fasta")
file_pat = re.compile(r"results/inserts/(.*)/consensus_(.*).fasta")
with bgzf.BgzfWriter(outfile_name) as output_handle:
    for fasta_file in infile_names:
        m = file_pat.match(fasta_file)
        assert m,fasta_file
        for record in SeqIO.parse(gzip.open(fasta_file, "rt"), "fasta"):
            if record.id == "Consensus":
                record.id = f"{m.group(1)}_{m.group(2)}_consensus"
                record.seq = record.seq.ungap()
                SeqIO.write(record, handle=output_handle, format="fasta")