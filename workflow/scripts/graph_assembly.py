#%%
import pandas as pd
import networkx as nx
import subprocess as sp
from pathlib import Path

FNAME = snakemake.input.ava_paf
INSERT_PAF = snakemake.input.insert_paf
INSERT_NAME = snakemake.wildcards.insert_name
OUTPUT_DIR = Path(snakemake.output.cluster_dir)
SLACK = snakemake.params.SLACK
SLOP = snakemake.params.SLOP

from snakemake.logging import logger

print = logger.info
if not OUTPUT_DIR.exists():
    OUTPUT_DIR.mkdir()


# Reading ava
p = sp.Popen(f"zcat {FNAME} |cut -f-10", stdout=sp.PIPE, shell=True)
ava = pd.read_table(p.stdout, header=None, sep="\t")

ava.columns = ["r1", "l1", "s1", "e1", "strand", "r2", "l2", "s2", "e2", "matches"]


inserts = pd.read_table(
    INSERT_PAF,
    sep="\t",
    header=None,
)

x = [
    "read_id",
    "length",
    "start",
    "end",
    "strand",
    "insert",
    "insert_length",
    "insert_start",
    "insert_end",
    "matches",
]
inserts.columns = x + list(inserts.columns[len(x) :])

all_reads = set(inserts.read_id)
read_status = pd.Series("OK",index=list(all_reads))
reads_remaning = all_reads.copy()

print("Total reads:", len(all_reads))


# Filter reads

idx_ins = (
    inserts.query("matches>(0.1*insert_length)")
    .sort_values("matches")
    .drop_duplicates("read_id", keep="last")
    .set_index("read_id")
)
print("Reads with insert", inserts.read_id.nunique())
print("Reads with insert of sensible (10%) size", idx_ins.index.nunique())


insensible_reads = all_reads - set(idx_ins.index)
reads_remaning -= insensible_reads
read_status[list(insensible_reads)] = "ShortInsert"
print(f"Filtered {len(insensible_reads)} reads")


# Filtering edges:

# Reads with good inserts:
ava_f = ava.sort_values("matches").drop_duplicates(["r1", "r2"], keep="last").copy()
ava_f = ava_f.loc[ava_f.r1.isin(idx_ins.index) & ava_f.r2.isin(idx_ins.index)]


ava_f = ava_f.set_index("r1")
ava_f["i1_start"] = idx_ins.loc[ava_f.index, "start"]
ava_f["i1_end"] = idx_ins.loc[ava_f.index, "end"]

ava_f = ava_f.reset_index().set_index("r2")
ava_f["i2_start"] = idx_ins.loc[ava_f.index, "start"]
ava_f["i2_end"] = idx_ins.loc[ava_f.index, "end"]
ava_f = ava_f.reset_index()

import numpy as np

ava_f["extra1"] = np.maximum((ava_f.i1_start - ava_f.s1), 0) + np.maximum(
    ava_f.e1 - ava_f.i1_end, 0
)
ava_f["extra2"] = np.maximum(ava_f.i2_start - ava_f.s2, 0) + np.maximum(
    ava_f.e2 - ava_f.i2_end, 0
)


ava_dat_w_contained = ava_f.copy()

print("Pairwise alignments between reads with sensible inserts", len(ava_f))


ava_contained_a = ava_f.loc[((ava_f.s1 < SLACK) & ((ava_f.l1 - ava_f.e1) < SLACK)) & ( (ava_f.l2-ava_f.l1)>3*SLACK ) ]
ava_contained_b = ava_f.loc[((ava_f.s2 < SLACK) & ((ava_f.l2 - ava_f.e2) < SLACK)) & ( (ava_f.l1-ava_f.l2)>3*SLACK )]


# Graph for later assignment to proper clusters of reads
contained_graph = nx.DiGraph()
contained_graph.add_edges_from((x.r1, x.r2) for _, x in ava_contained_a.iterrows())
contained_graph.add_edges_from((x.r2, x.r1) for _, x in ava_contained_b.iterrows())

contained = set(ava_contained_a.r1)
contained.update(set(ava_contained_b.r2))

ava_contained = pd.concat([ava_contained_a, ava_contained_b])

reads_remaning -= contained
read_status[list(contained)] = "PreliminaryContained"
print("contained", len(contained))
ava_f = ava_f.loc[~(ava_f.r1.isin(contained) | ava_f.r2.isin(contained))]


print("Reads contained in some other", len(contained))
print(
    "Pairwise alignments between reads with sensible inserts and not fully contained",
    len(ava_f),
)

read_end_overlap_idx = ((ava_f.s1 < SLACK) | ((ava_f.l1 - ava_f.e1) < SLACK)) & (
    (ava_f.s2 < SLACK) | ((ava_f.l2 - ava_f.e2) < SLACK)
)
ava_f = ava_f.loc[read_end_overlap_idx]
print(
    "Pairwise alignments between read ends with sensible inserts and not fully contained",
    len(ava_f),
)

reads_with_edge_overlap = set(ava_f.r1) | set(ava_f.r2)
reads_without_edge_overlap = reads_remaning - reads_with_edge_overlap
print("Reads without edge overlap", len(reads_without_edge_overlap))
reads_remaning -= reads_without_edge_overlap

read_status[list(reads_without_edge_overlap)] = "NoEdgeOverlap"
print("Remaining reads", len(reads_remaning))

ava_dat = ava_f.copy()

## Main filter for alignments:  Need at least SLOP bases aligned on one side of the insert
ava_f = ava_f.loc[
    ((ava_f.s1 < ava_f.i1_start - SLOP) | (ava_f.e1 > ava_f.i1_end + SLOP))
    & ((ava_f.s2 < ava_f.i2_start - SLOP) | (ava_f.e2 > ava_f.i2_end + SLOP))
]


reads_with_covered_insert = set(ava_f.r1) | set(ava_f.r2)
reads_with_insufficient_insert_cover = reads_remaning - reads_with_covered_insert

read_status[list(reads_with_insufficient_insert_cover)] = "NoAlignmentMoreThanInsert"

print(
    "Reads with insufficent post insert extension",
    len(reads_with_insufficient_insert_cover),
)
reads_remaning -= reads_with_insufficient_insert_cover

print(
    "Pairwise alignments overlapping insert by",
    SLOP,
    ",between reads with sensible inserts",
    len(ava_f),
)


# %%

G = nx.Graph()

G.add_edges_from(
    (
        x.r1,
        x.r2,
        {
            "color": "red",
            "weight": (x.matches / max(x.e1 - x.s1, x.e2 - x.s2) - 0.5) * 3,
        },
    )
    # (x.r1,x.r2)
    for _, x in ava_f.drop_duplicates(["r1", "r2"]).iterrows()
    if x.r1 != x.r2
)

components = list(nx.connected_components(G))

print("Reads remaining", G.number_of_nodes())
print("Inserts observed", len(components))
#%%


# Assigning cluster identities:
cluster_names = {}
for i, C in enumerate(components):
    cluster_names.update({x: i for x in C})
    for x in C:
        G.nodes[x]["cluster"] = i
cluster_names = pd.Series(cluster_names)


#%%

for n in contained_graph.nodes:
    try:
        contained_graph.nodes[n]["cluster"] = G.nodes[n]["cluster"]
    except KeyError:
        pass

for source_node in contained_graph.nodes:
    if "cluster" in contained_graph.nodes[source_node]:
        continue
    this_cluster = None
    for t_node in nx.dfs_postorder_nodes(contained_graph, source=source_node):
        t_cluster = contained_graph.nodes[t_node].get("cluster", -1)
        if this_cluster is None and t_cluster != -1:
            # If we found first node with known cluster status
            this_cluster = t_cluster
        elif this_cluster is not None:  # If we know which cluster we are in
            if t_cluster == -1:
                contained_graph.nodes[t_node]["cluster"] = this_cluster
            elif (
                t_cluster != this_cluster
            ):  # If we found conflicting cluster. Cant proceed
                break


found_cluster = set()
ambiguous_cluster = set()

contained_G = G.copy()

for r in contained:
    if contained_graph.nodes[r].get("cluster", -1) != -1:
        found_cluster.add(r)
        read_status[r] = "RecoveredContained"
        contained_G.add_edges_from(
            (r, n, {"weight": 1, "color": "green"})
            for n in contained_graph.neighbors(r)
        )
    else:
        ambiguous_cluster.add(r)
        read_status[r] = "AmbiguousContained"

print(
    f"Found clusters for {len(found_cluster)} edges. Ambiguous {len(ambiguous_cluster)}."
)


missing_nodes = set(inserts.read_id) - set(contained_G.nodes)

OUTPUT_DIR.joinpath("excluded.lst").open("wt").write("\n".join(missing_nodes))
for i, C in enumerate(nx.connected_components(contained_G)):
    OUTPUT_DIR.joinpath(f"insert_{i}.lst").open("wt").write("\n".join(C))
