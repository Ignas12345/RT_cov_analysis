import argparse, os , re, sys
import numpy as np

sys.path.append('/tmp/Mazutislab-out/Ignas/RT_comparison/code')
from get_max_poly_x_run_in_bins import get_bin_ends

def get_exon_positions(exon_gtf_path):
    exons_per_transcript = {}
    for i, line in enumerate(open(exon_gtf_path)):
        line = line.split('\t')
        if 'havana_transcript "' in line[8]:
            meta_data = line[8].split('; ')
            for item in meta_data:
                if 'havana_transcript "' in item:
                    key = item.strip()
                    key = re.search(r'havana_transcript\s*"([^"]+)"', key).group(1)
                    if key not in exons_per_transcript:
                        exons_per_transcript[key] = {}
                        exons_per_transcript[key]['exon_start_positions'] = []
                        exons_per_transcript[key]['exon_end_positions'] = []
                    break 
            exons_per_transcript[key]['exon_start_positions'].append(int(line[3]))
            exons_per_transcript[key]['exon_end_positions'].append(int(line[4]))
    return exons_per_transcript

def characterize_exons_by_bins(exons_per_transcript, gene_bed_file, n_bins):

    for i, line in enumerate(open(gene_bed_file)):
        line = line.split('\t')
        gene_start = int(line[1])
        gene_end = int(line[2])
        strand = line[5]
        gene_length = gene_end - gene_start

    n_bins = 100
    bin_ends = get_bin_ends(gene_length, n_bins)
    #if the gene is on the negative strand, then the 5' end will be corresponding to the last bin, as the exon start positions are counted from the 5' end of the positive strand
    if strand == '-':
            bin_ends = bin_ends[::-1]

    binned_transcript_exon_maps = {}
    for key in exons_per_transcript:
        binned_transcript_exon_maps[key] = []
        shifted_starts = np.array(exons_per_transcript[key]['exon_start_positions']) - gene_start
        shifted_ends = np.array(exons_per_transcript[key]['exon_end_positions']) - gene_start
        #print(shifted_starts)
        #print(shifted_ends)
        #print(bin_ends)
        for i, bin_end in enumerate(bin_ends):
            in_exon = np.where(
                (shifted_starts <= bin_end) & (shifted_ends >= bin_end)
            )
            if len(in_exon[0]) > 0:
                binned_transcript_exon_maps[key].append(1)
            else:
                binned_transcript_exon_maps[key].append(0)
    return binned_transcript_exon_maps

if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage='This code takes the exon lines in a gtf file and produces binned coverage of a gene, where 1 indicates an exonic region and 0 non-exonic. Only transcripts that have a havana trascript id in the gtf file are used.')
    parser.add_argument('-e', '--exon_gtf', type=str, help="The input gtf file with the exons. Exons will be grouped together with havana transcript ID used as key. Also, better to ensure that all transcrpits come from the same strand (although maybe that's always the case).")
    parser.add_argument('-g', '--gene_bed', type=str, help='Gene bed file that contains the strand, start and end positions of the gene.')
    parser.add_argument('-n', '--n_bins', type=int, default = 100, help='The number of bins that will be used for representing the exonic coverage.')

    args = parser.parse_args()
    exons_per_transcript = get_exon_positions(args.exon_gtf)
    exon_in_bins = characterize_exons_by_bins(exons_per_transcript=exons_per_transcript, gene_bed_file=args.gene_bed, n_bins=args.n_bins)

    for key in exon_in_bins:
        print(key)
        print(exon_in_bins[key])