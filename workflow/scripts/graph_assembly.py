#%%
from builtins import KeyError
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

inserts_columns = [
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
inserts.columns = inserts_columns + list(inserts.columns[len(inserts_columns) :])
Q = "15741ece-00aa-4814-bff1-6f34b5a3198b"
Q_node = ("15741ece-00aa-4814-bff1-6f34b5a3198b",
        "Lenti-Cas9-2A-Blast",
        25367)
assert Q in set(inserts.read_id), "About line 48"

#%%
UNSPECIFIED = ""
read_annot = pd.Series(UNSPECIFIED, index=inserts.read_id.drop_duplicates())
read_annot.name = "faith"
# Filter reads that are almost completely the inserts

_insert_reads = set(inserts.read_id)
filtered_inserts = inserts.query(f"(end-start)< (length - 4*{SLACK})")
assert (
    read_annot[list(_insert_reads - (set(filtered_inserts.read_id)))] == UNSPECIFIED
).all()
read_annot[
    list(_insert_reads - (set(filtered_inserts.read_id)))
] = "Almost completely insert. "
_insert_reads = set(filtered_inserts.read_id)

filtered_inserts = inserts.query("length>insert_length")
assert (
    read_annot[list(_insert_reads - (set(filtered_inserts.read_id)))] == UNSPECIFIED
).all()
read_annot[
    list(_insert_reads - (set(filtered_inserts.read_id)))
] = "Shorter than insert. "
_insert_reads = set(filtered_inserts.read_id)
assert Q in set(filtered_inserts.read_id)
# Filter matches that match less than 10% of the insert
filtered_inserts = filtered_inserts.query("matches>(0.1*insert_length)")[
    inserts_columns
]
assert (
    read_annot[list(_insert_reads - (set(filtered_inserts.read_id)))] == UNSPECIFIED
).all()
read_annot[
    list(_insert_reads - (set(filtered_inserts.read_id)))
] = "Matches less than 10% of insert. "
_insert_reads = set(filtered_inserts.read_id)

# filter bad alignments

filtered_inserts["match_prop"] = filtered_inserts.matches / (
    filtered_inserts.end - filtered_inserts.start
)
filtered_inserts = filtered_inserts.query("match_prop>0.85")

assert (
    read_annot[list(_insert_reads - (set(filtered_inserts.read_id)))] == UNSPECIFIED
).all()
read_annot[
    list(_insert_reads - (set(filtered_inserts.read_id)))
] = "Match identity < 85%. "
_insert_reads = set(filtered_inserts.read_id)


filtered_inserts = filtered_inserts.sort_values("matches").drop_duplicates(
    "read_id", keep="last"
)


cr1 = ava.query(f"(r1!=r2) &(l1<(e1-s1+ 2*{SLACK}) )")
cr2 = ava.query(f"(r1!=r2) &(l2<(e2-s2+ 2*{SLACK}) )")
contained_reads = set(cr1.r1)
contained_reads |= set(cr2.r2)


# r1 is contained in r2  :   edge r2->r1
contains_graph = nx.DiGraph()
contains_graph.add_edges_from((x.r2, x.r1) for _, x in cr1.iterrows())
contains_graph.add_edges_from((x.r1, x.r2) for _, x in cr2.iterrows())

filtered_inserts_all = filtered_inserts.copy()
filtered_inserts = filtered_inserts.loc[~filtered_inserts.read_id.isin(contained_reads)]
#%%
# Merging with all-vs-all alignments. Require insert _overlap_ alignment:
assert Q in set(filtered_inserts.read_id)
ava_m = pd.merge(
    ava.query("(r1!=r2)"),
    filtered_inserts,
    left_on="r1",
    right_on="read_id",
    suffixes=("_ava", ""),
).query("(s1<end) & (e1>start)")

assert Q in (set(ava_m.r1) | set(ava_m.r2))
ava_m = pd.merge(
    ava_m, filtered_inserts, left_on="r2", right_on="read_id", suffixes=("1", "2")
).query("(s2<end2) & (e2>start2)")

# Require that insert is in the same strand as the a-v-a alignment

ava_m = ava_m.loc[
    ava_m.strand_ava
    == (ava_m.strand1 == ava_m.strand2).replace({True: "+", False: "-"})
]


# Lift over read2 insert locations to read1
ava_m["r1_start2"] = pd.Series((ava_m.s1 - ava_m.s2) + ava_m.start2).where(
    ava_m.strand_ava == "+", (ava_m.e1 - (ava_m.end2 - ava_m.s2))
)
ava_m["r1_end2"] = pd.Series((ava_m.s1 - ava_m.s2) + ava_m.end2).where(
    ava_m.strand_ava == "+", (ava_m.e1 - (ava_m.start2 - ava_m.s2))
)

# Require overlap in lifted over coordinates also.
ava_m = ava_m.loc[(ava_m.start1 < ava_m.r1_end2) & (ava_m.end1 > ava_m.r1_start2)]
ava_m["insert_overlap_bp"] = ava_m[["end1", "r1_end2"]].min(axis=1) - ava_m[
    ["start1", "r1_start2"]
].max(axis=1)


# Extra space in alignment before and after the insert
extra_pre = ava_m.start1 - ava_m.s1
extra_post = ava_m.e1 - ava_m.end1
ava_m["extra_5p"] = extra_pre.where(ava_m.strand1 == "+", extra_post)
ava_m["extra_3p"] = extra_post.where(ava_m.strand1 == "+", extra_pre)


# Length of unaligned tips
_tips = ava_m[["strand_ava", "s1", "s2"]], (
    ava_m[["l1", "l2"]] - ava_m[["e1", "e2"]].values
)
_tips = pd.concat(_tips, axis=1)
ava_m["tip3"] = (
    _tips[["s1", "s2"]]
    .min(axis=1)
    .where(ava_m.strand_ava == "+", _tips[["s1", "l2"]].min(axis=1))
)
ava_m["tip5"] = (
    _tips[["l1", "l2"]]
    .min(axis=1)
    .where(ava_m.strand_ava == "+", _tips[["s2", "l1"]].min(axis=1))
)
assert Q in (set(ava_m.r1) | set(ava_m.r2))
#%%
ava_back = ava_m.copy()
#%%

# ava_m = ava_back.copy()

# Limit length of unaligned tips
ava_m = ava_m.loc[ava_m[["tip3", "tip5"]].max(axis=1) < SLACK]


# Require significant overlap of the insert
ava_m = ava_m.loc[ava_m.insert_overlap_bp > SLOP]

# Require sigificant overlap of the non-insert
ava_m = ava_m.loc[ava_m[["extra_3p", "extra_5p"]].max(axis=1) > SLOP]

# Maybe something..
_n_pre = len(ava_m)
ava_m = ava_m.loc[ava_m.matches_ava >= SLOP + 1.1 * ava_m.insert_overlap_bp]

logger.info("Filtered: {}".format(_n_pre - len(ava_m)))
#%%


G = nx.Graph()
G.add_nodes_from(
    ((x.read_id, x.insert, x.start), x.to_dict())
    for _, x in filtered_inserts.iterrows()
)
G.add_edges_from(
    (
        (x.r1, x.insert1, x.start1),
        (x.r2, x.insert2, x.start2),
        dict(color="red", weight=1, **x),
    )
    for _, x in ava_m.iterrows()
    if x.r1 != x.r2
)

components = list(sorted(nx.connected_components(G), key=len, reverse=True))
non_unit_components = {}
for i, C in enumerate(components):
    non_unit_components[i] = C
    for n in C:
        G.nodes[n]["cluster_id"] = i


        
#%%
# Propagate cluster id through the contained reads:
for cluster_id, C in non_unit_components.items():
    for n in C:
        read_id = n[0]
        if read_id in contains_graph.nodes:
            # logger.info(f"Read {read_id} in contains_graph.")
            contains_graph.nodes[read_id].setdefault("cluster_id", set()).add(
                cluster_id
            )
            for contained_read in nx.dfs_preorder_nodes(contains_graph, read_id):
                contains_graph.nodes[contained_read].setdefault(
                    "cluster_id", set()
                ).add(cluster_id)
                # logger.info(f"Added cluster {cluster_id} for {contained_read}")


cluster_includes = {}
for n in contains_graph.nodes:
    cluster_ids = contains_graph.nodes[n].get("cluster_id", set())
    # logger.info(f"contains_graph {n}  {cluster_ids}")

    if len(cluster_ids) == 1:
        cluster_id = next(iter(cluster_ids))
        cluster_includes.setdefault(cluster_id, set()).add(n)
        if contains_graph.in_degree(n) >0:
            read_annot[n] += "Contained in a read in a cluster. "
    elif len(cluster_ids) == 0:
        if contains_graph.in_degree(n) >0:
            read_annot[n] += "Contained in a dangling/unclustered read. "
    else:
        if contains_graph.in_degree(n) >0:
            read_annot[n] += "Contained in many clusters. "
    # if len(cluster_ids)==0:
    #    print(n,"MISSING")
    # elif len(cluster_ids)>1:
    #    print(n,",".join(map(str,cluster_ids)))
#%%
all_nodes = set(inserts.set_index(["read_id", "insert", "start"]).index)

import functools

missing_nodes = all_nodes - functools.reduce(
    lambda x, y: x | y, non_unit_components.values(), set()
)
OUTPUT_DIR.joinpath("excluded.lst").open("wt").write(
    "\n".join("\t".join(map(str, x)) for x in missing_nodes)
)

filtered_inserts.index = filtered_inserts.set_index(
    ["read_id", "insert", "start"]
).index

for i, C in non_unit_components.items():
    # OUTPUT_DIR.joinpath(f"reads_cluster_{i}.lst").open("wt").write("\n".join(C))

    # C_inserts = inserts.loc[insert_idx.get_loc(list(C))]

    # C_inserts = filtered_inserts.loc[list(C)]
    read_annot[[n for n, _, _ in C]] = "Included in a cluster"
    C_inserts = (
        pd.concat(
            [
                filtered_inserts.loc[list(C)],
                filtered_inserts_all.loc[
                    filtered_inserts_all.read_id.isin(cluster_includes.get(i, set()))
                ],
            ]
        )
        .sort_values("matches", ascending=False)
        .drop_duplicates("read_id")
    )
    assert (
        not C_inserts.read_id.duplicated().any()
    ), f"Cluster {i} has duplicated reads: {','.join(C_inserts.read_id)}"
    if len(C_inserts.read_id)==1:
        read_annot[next(iter(C))[0]] += "Single read component. "
    else:
        C_inserts.to_csv(
            OUTPUT_DIR.joinpath(f"reads_insert_{i}.paf"),
            sep="\t",
            index=False,
            header=False,
        )

from networkx.readwrite import json_graph
import json

json.dump(
    json_graph.node_link_data(G), OUTPUT_DIR.joinpath("merge_graph.json").open("wt")
)
logger.info(list(G.neighbors(Q_node)))
# assert Q in contains_graph.nodes
# assert read_annot[Q]!=UNSPECIFIED
read_annot.to_csv(OUTPUT_DIR.joinpath("read_annotation.tsv"), sep="\t")
