import os
import csv
import numpy as np


def filter_gene_hits(input, p_percent, hits_ratio2max, beta_percent):
    # Read original gene hits table
    header = []
    tmp_tbl = []
    with open(input, 'r') as i_csv:
        reader = csv.reader(i_csv)
        for index, line in enumerate(reader):
            if index == 0:
                header = line
                i_p = line.index('maxp')
                i_hits = line.index('hits')
                i_beta = -1

            elif len(line[1].split('_')[0]) != 8:
                tmp_tbl.append(line)

    # Calculate cutoff
    p_cutoff = np.percentile([float(line[i_p]) for line in tmp_tbl], p_percent)
    hits_cutoff = max(int(line[i_hits]) for line in tmp_tbl)*hits_ratio2max
    beta_cutoff = np.percentile([float(line[i_beta]) for line in tmp_tbl], beta_percent)

    filtered_hits = []

    # Apply filters
    for line in tmp_tbl:
        if float(line[i_p]) >= p_cutoff and float(line[i_hits]) >= hits_cutoff and float(line[i_beta]) >= beta_cutoff:
            filtered_hits.append(line)
        else:
            pass

    return header, filtered_hits


for root, dirs, files in os.walk('./all_gene_hits'):
    gene_hits_lst = [file for file in files]

if not os.path.exists('./filtered_gene_hits/'):
    os.mkdir('./filtered_gene_hits/')

for in_csv in gene_hits_lst:
    header = filter_gene_hits('./all_gene_hits/'+in_csv, 90, 0.1, 10)[0]
    filtered_tbl = filter_gene_hits('./all_gene_hits/'+in_csv, 90, 0.1, 10)[1]
    out_csv = in_csv.split('.')[0]+'_filtered.csv'

    with open('./filtered_gene_hits/'+out_csv, 'w') as o_csv:
        writer = csv.writer(o_csv)
        writer.writerow(header)
        for line in filtered_tbl:
            writer.writerow(line)



