from builtins import KeyError
from Bio import SearchIO,SeqIO,bgzf
import gzip

input_psl = snakemake.input.inserts_psl
input_fasta = snakemake.input.fasta[0]

output_bed = snakemake.output.insert_bed
output_fasta = snakemake.output.masked_fasta

qresult = SearchIO.read(input_psl, "blat-psl")



with bgzf.BgzfWriter(output_fasta) as output_handle, open(output_bed,"wt") as output_bed_handle:
    for record in SeqIO.parse(gzip.open(input_fasta, "rt"), "fasta"):
        try:
            qr =qresult[record.id]
            best_hit = sorted(qr.hsps,key=lambda x:-x.match_num)[0]
            seq_id,start,end = qr.id,best_hit.hit_start,best_hit.hit_end
            print(seq_id,start,end,len(record),sep="\t",file=output_bed_handle)


            s = record.seq.tomutable()
            s[start:end]="N"*(end-start)
            record.seq = s
            SeqIO.write(record, handle=output_handle, format="fasta")
        except KeyError:
            pass
#%%
# from Bio import SearchIO 
# psl="output.psl"

# qresult = SearchIO.read(psl, "blat-psl")
# qresult