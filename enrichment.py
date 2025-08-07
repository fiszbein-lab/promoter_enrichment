import csv

from statistics import mean

import numpy as np

from scipy.stats import chisquare
from statsmodels.stats.multitest import multipletests


for cut_off in range(50, 100, 5):
    print(cut_off)
    THRESH = cut_off

    seq_repr = dict()
    seq_group = dict()

    with open("out/normalized-count.tsv", 'r') as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)

        for row in reader:
            seq, group, *counts = row
            counts = map(float, counts)

            seq_repr[seq] = [mean(d) for d in np.array_split(list(counts), 5)]
            seq_group[seq] = group

    seqs = list()
    freqs = list()
    ps = list()

    for seq in seq_repr:
        seq_freq = seq_repr[seq]

        if sum(seq_freq) > 0:
            stat, p = chisquare(seq_freq)

            seqs.append(seq)
            freqs.append(seq_freq)
            ps.append(p)

    _, corrected_ps, *_ = multipletests(ps, method='bonferroni')

    with open(f"out/enrichment/uniform-exp_enrichment>{THRESH}.tsv",
              'w') as f:
        f.write(f"seqname\t"
                f"group\t"
                f"gate-1_norm-percent\t"
                f"gate-2_norm-percent\t"
                f"gate-3_norm-percent\t"
                f"gate-4_norm-percent\t"
                f"gate-5_norm-percent\t"
                f"gate-1_norm-count\t"
                f"gate-2_norm-count\t"
                f"gate-3_norm-count\t"
                f"gate-4_norm-count\t"
                f"gate-5_norm-count\t"
                f"corrected-p\t"
                f"enriched-gate\n")

        for seq, freq, corrected_p in zip(seqs, freqs, corrected_ps):
            group = seq_group[seq]

            enriched_gate = 0
            if corrected_p < 0.05:
                for gate, count in enumerate(freq):
                    if (count / sum(freq) * 100) > THRESH:
                        enriched_gate = gate + 1

            f.write(f"{seq}\t{group}\t" +
                    f"\t".join([str(f/sum(freq)*100) for f in freq]) + "\t" +
                    f"\t".join([str(f) for f in freq]) +
                    f"\t{corrected_p}"
                    f"\t{enriched_gate}\n")
