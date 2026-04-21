configfile: "config.yaml"

SAMPLES = config["samples"]
DATA_DIR = config["data_dir"]
OUT_DIR = config["out_dir"]
TSO = config["tso"]
GTF = config["gtf"]

# Input fastq filename suffixes (relative to DATA_DIR/{sample})
CDNA_SUFFIX = "_cdna_001.fastq.gz"
BC_SUFFIX = "_bc_001.fastq.gz"

# Output fastq filename suffixes for TSO-containing reads (relative to OUT_DIR/{sample}/reads_with_tso)
CDNA_TRIMMED_SUFFIX = "_cdna_001_trimmed.fastq"
BC_OUT_SUFFIX = "_bc_001.fastq"


rule all:
    input:
        expand(
            "{out_dir}/{sample}/reads_with_tso/{sample}" + CDNA_TRIMMED_SUFFIX,
            out_dir=OUT_DIR,
            sample=SAMPLES,
        ),
        expand(
            "{out_dir}/{sample}/reads_with_tso/{sample}" + BC_OUT_SUFFIX,
            out_dir=OUT_DIR,
            sample=SAMPLES,
        ),
        OUT_DIR + "/genome_and_annotations/single_isoform_genes.txt",
        expand(
            "{out_dir}/{sample}/star/{sample}_Aligned.sortedByCoord.out.bam",
            out_dir=OUT_DIR,
            sample=SAMPLES,
        ),


rule detect_tso:
    """
    Use cutadapt to detect the TSO at the 5' end of cDNA reads (mismatches allowed).
    Only read pairs where the cDNA read contained the TSO are written to the output
    fastq files. The TSO is trimmed from the cDNA reads; the paired barcode reads are
    written unchanged. The output files are intended for downstream inspection,
    trimming, and remapping.
    """
    input:
        cdna=DATA_DIR + "/{sample}" + CDNA_SUFFIX,
        bc=DATA_DIR + "/{sample}" + BC_SUFFIX,
    output:
        cdna=OUT_DIR + "/{sample}/reads_with_tso/{sample}" + CDNA_TRIMMED_SUFFIX,
        bc=OUT_DIR + "/{sample}/reads_with_tso/{sample}" + BC_OUT_SUFFIX,
    params:
        tso=TSO,
        error_rate=config["tso_error_rate"],
    threads: 4
    conda:
        "envs/cutadapt.yaml"
    log:
        os.path.join(OUT_DIR, "logs/detect_tso/{sample}.log")
    shell:
        """
        cutadapt \
            -g {params.tso} \
            --discard-untrimmed \
            -m 30 \
            -j {threads} \
            -o {output.cdna} \
            -p {output.bc} \
            {input.cdna} {input.bc}
        """


rule star_mapping:
    """
    Map TSO-containing cDNA reads with STARsolo.
    Accepts the filtered and TSO-trimmed fastq files produced by the detect_tso rule.
    """
    input:
        cdna=rules.detect_tso.output.cdna,
        bc=rules.detect_tso.output.bc,
    output:
        bam=OUT_DIR + "/{sample}/star/{sample}_Aligned.sortedByCoord.out.bam",
    params:
        genome_dir=config["genome_dir"],
        bc1_whitelist=os.path.join(config["whitelist_dir"],"bc1_list.txt"),
        bc2_whitelist=os.path.join(config["whitelist_dir"],"bc2_list.txt"),
        bc3_whitelist=os.path.join(config["whitelist_dir"],"bc3_list.txt"),
        prefix=OUT_DIR + "/{sample}/star/{sample}_",
    threads: 14
    conda:
        "envs/star.yaml"
    log:
        os.path.join(OUT_DIR, "logs/star_mapping/{sample}.log")
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {params.genome_dir} \
            --readFilesIn {input.cdna} {input.bc} \
            --soloCBwhitelist {params.bc1_whitelist} {params.bc2_whitelist} {params.bc3_whitelist} \
            --readFilesCommand cat \
            --outFileNamePrefix {params.prefix} \
            --alignIntronMax 0 \
            --outSAMunmapped Within \
            --outSAMtype BAM SortedByCoordinate \
            --outBAMsortingBinsN 20 \
            --outSAMattributes NH HI nM AS CR UR CB UB sS sQ sM GX GN \
            --soloType CB_UMI_Complex \
            --soloUMIdedup Exact \
            --soloCBposition 0_0_0_7 0_8_0_17 0_18_0_25 \
            --soloUMIposition 0_26_0_35 \
            --soloBarcodeReadLength 1 \
            --soloCBmatchWLtype EditDist_2 \
            --soloMultiMappers Uniform

            samtools index {output.bam}
        """
        
rule filter_single_isoform_genes:
    """
    Parse a GTF annotation file and write out all genes that have exactly one
    annotated transcript (isoform).

    Each output line contains three space-separated fields:
        gene_id  gene_name  transcript_id

    The output file can be used downstream to restrict read-distance analysis
    to genes whose exon/intron structure is unambiguous.
    """
    input:
        gtf=GTF,
    output:
        OUT_DIR + "/genome_and_annotations/single_isoform_genes.txt",
    run:
        import re
        from collections import defaultdict
        import os

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)

        # gene_id -> set of transcript_ids seen for that gene
        gene_transcripts = defaultdict(set)
        # gene_id -> gene_name (last value wins, but it's stable per gene)
        gene_names = {}

        with open(input.gtf) as fh:
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

        with open(output[0], "w") as out:
            for gene_id, transcripts in sorted(gene_transcripts.items()):
                if len(transcripts) == 1:
                    transcript_id = next(iter(transcripts))
                    gene_name = gene_names[gene_id]
                    out.write(f"{gene_id} {gene_name} {transcript_id}\n")
