import pandas as pd

ava = pd.read_table(
    snakemake.input.ava_paf, sep="\t", usecols=list(range(10)), header=None
)
ava.columns=["read_id","length","start","end","strand","insert","insert_length","insert_start","insert_end","matches"]

ava = ava.set_index("read_id")
out_files = {x.split("/")[-1].split(".")[-2]:open(x,"wt") for x in snakemake.output.insert_bed}

for read_id,read_data in ava.iterrows():
    if (read_data.insert_end - read_data.insert_start)>snakemake.params.insert_frac*read_data.insert_length:
        m_start,m_end = read_data.start,read_data.end

        if read_data.strand == "+":
            m_start -=read_data.insert_start
            m_end += (read_data.insert_length-read_data.insert_end)
        else:
            m_start -=(read_data.insert_length-read_data.insert_end)
            m_end += read_data.insert_start

        m_start = max(m_start,0)
        m_end = min(m_end,read_data.length)
        out_bed = out_files[read_data.insert]
        print(read_id,m_start,m_end,sep="\t",file=out_bed)

