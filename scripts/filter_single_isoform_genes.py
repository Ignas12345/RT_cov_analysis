"""
Parse a GTF annotation file and write out all genes that have exactly one
annotated transcript (isoform).

Each output line contains three space-separated fields:
    gene_id  gene_name  transcript_id

The output file can be used downstream to restrict read-distance analysis
to genes whose exon/intron structure is unambiguous.
"""

import re
from collections import defaultdict

# gene_id -> set of transcript_ids seen for that gene
gene_transcripts = defaultdict(set)
# gene_id -> gene_name (last value wins, but it's stable per gene)
gene_names = {}

with open(snakemake.input.gtf) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "transcript":
            continue
        attr = fields[8]

        m = re.search(r'gene_id "([^"]+)"', attr)
        if not m:
            continue
        gene_id = m.group(1)

        m = re.search(r'transcript_id "([^"]+)"', attr)
        if not m:
            continue
        transcript_id = m.group(1)

        m = re.search(r'gene_name "([^"]+)"', attr)
        gene_name = m.group(1) if m else gene_id

        gene_transcripts[gene_id].add(transcript_id)
        gene_names[gene_id] = gene_name

with open(snakemake.output[0], "w") as out:
    for gene_id, transcripts in sorted(gene_transcripts.items()):
        if len(transcripts) == 1:
            transcript_id = next(iter(transcripts))
            gene_name = gene_names[gene_id]
            out.write(f"{gene_id} {gene_name} {transcript_id}\n")
