configfile: "config.yaml"

SAMPLES = config["samples"]
DATA_DIR = config["data_dir"]
OUT_DIR = config["out_dir"]
TSO = config["tso"]

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
    shell:
        """
        cutadapt \
            -g {params.tso} \
            --discard-untrimmed \
            -e {params.error_rate} \
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
        bc1_whitelist=config["bc1_whitelist"],
        bc2_whitelist=config["bc2_whitelist"],
        bc3_whitelist=config["bc3_whitelist"],
        prefix=OUT_DIR + "/{sample}/star/{sample}_",
    threads: 14
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
        """
