import csv
import glob
import os
import re

from collections import defaultdict

from Bio import SeqIO

# Map sequence to group for output-writing to ease later decomposition.
seq_group = dict()
for group in glob.iglob("ref/pool-groups/*"):
    group_name = os.path.basename(group).removesuffix(".fa")

    for record in SeqIO.parse(group, 'fasta'):
        seq_group[record.id] = group_name


# Use gate_enrichment lib. as a proxy for integration efficiencies.
enrichment_seq_repr = dict()
enrichment_lib_size = 0

with open("gate_enrichment-count.tsv", 'r') as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader)

    for row in reader:
        seq, count = row
        if seq:
            count = int(count)
            if count > 0:
                enrichment_seq_repr[seq] = count
                enrichment_lib_size += count

sorted_seq_repr = defaultdict(lambda: defaultdict(dict))
sorted_lib_size = defaultdict(lambda: defaultdict(lambda: 0))

for path in glob.iglob("data/*"):
    if os.path.isdir(path):
        gate = int(re.search(r'\d', os.path.basename(path)).group(0))

        with open(f"{path}/gate-{gate}_rep-shared_counts.tsv", 'r') as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)

            for row in reader:
                seq, group, *counts = row
                counts = map(int, counts)

                for rep, count in enumerate(counts):
                    sorted_seq_repr[gate][rep+1][seq] = count
                    sorted_lib_size[gate][rep+1] += count

with open("normalized-count.tsv", 'w') as f:
    gates = [f"gate-{i}_rep-{j}" for i in range(1, 6) for j in range(1, 4)]
    f.write(f"seqname\tgroup\t" + "\t".join(gates) + "\n")

    for seq in enrichment_seq_repr:
        group = seq_group[seq]
        row = [seq, group]

        for gate in range(1, 6):
            for rep in range(1, 4):
                # For each sequence's count, `X`, we take
                #     `(X_s / L_s) / (X_e / L_e)`
                # where the `s` subscript indicates a sorted lib., the `e`
                # subscript indicates the gate_enrichment lib., and `L` is the lib.
                # size -- the total number of mapped reads in that lib.
                X_s = sorted_seq_repr[gate][rep].get(seq, 0)
                L_s = sorted_lib_size[gate][rep]
                X_e = enrichment_seq_repr[seq]
                L_e = enrichment_lib_size

                row.append(str(round((X_s / L_s) / (X_e / L_e), 4)))

        f.write("\t".join(row) + "\n")
