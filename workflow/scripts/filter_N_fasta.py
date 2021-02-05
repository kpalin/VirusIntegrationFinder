#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" 

Created at Tuesday 02 February 2021  18:00 by Kimmo Palin <kpalin@helsinki.fi>
"""

__version__ = "0.1"

from os import EX_OK

import sys


def parse_args(args):
    import argparse

    description = "\n".join(__doc__.splitlines()[:-1]).strip()

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-i", "--input", help="Input file [default:%(default)s]", default="/dev/stdin"
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file [default:%(default)s]",
        default="/dev/stdout",
    )

    version = f"%(prog)s {__version__}"
    parser.add_argument("--version", action="version", version=version)

    parser.add_argument(
        "-V",
        "--verbose",
        default=False,
        action="store_true",
        help="Be more verbose with output",
    )

    args = parser.parse_args()

    import logging

    if args.verbose:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s:%(funcName)s:%(levelname)s:%(message)s",
        )

    return args


from Bio import SeqIO
from Bio.Seq import Seq
import gzip
maskN=str.maketrans("acgt","NNNN")

try:
    infile_name = snakemake.input[0]
    outfile_name = snakemake.output[0]
except NameError:
    args = parse_args(sys.argv)
    infile_name = args.input
    outfile_name = args.output


with open(outfile_name, "w") as output_handle:
   
    for record in SeqIO.parse(gzip.open(infile_name, "rt"), "fasta"):
        record.seq = Seq(str(record.seq).translate(maskN))
        if len(record.seq) > 1000:
            SeqIO.write(record, output_handle, "fasta")
