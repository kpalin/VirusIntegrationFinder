"""Read fasta multiple alignment, filter consensus sites with less than 51% coverage and output a
fasta sequence of the consensus without the gap characters.
"""

from Bio import SeqIO, bgzf
import gzip
import re
import numpy as np
import pdb
from Bio.Seq import Seq

from snakemake.logging import logger

infile_names = snakemake.input
outfile_name = snakemake.output[0]

file_pat = re.compile(r"results/inserts/([^ /]+)/consensus_([ ]*)[.]fasta")
file_pat = re.compile(r"results/inserts/(.*)/consensus_(.*).fasta")

with bgzf.BgzfWriter(outfile_name) as output_handle:
    for fasta_file in infile_names:
        m = file_pat.match(fasta_file)
        assert m, fasta_file
        n_reads_counted = 0
        for record in SeqIO.parse(gzip.open(fasta_file, "rt"), "fasta"):
            if record.id == "Consensus":
                record.id = f"{m.group(1)}_{m.group(2)}_reads{n_reads_counted}_consensus"
                # pdb.set_trace()
                mseq = record.seq.tomutable()
                for i in np.nonzero(coverages < min(3.5,(0.51 * n_reads_counted)))[0]:
                    mseq[int(i)] = "-"
                record.seq = mseq.toseq()
                record.seq = record.seq.ungap()
                if len(record.seq) > 10:
                    SeqIO.write(record, handle=output_handle, format="fasta")
                else:
                    logger.info(str(record))
            else:
                if n_reads_counted == 0:
                    coverages = np.zeros(len(record.seq))
                seq_a = np.array(str(record.seq), dtype="|c")
                coverages[seq_a != b"-"] += 1
                n_reads_counted += 1
