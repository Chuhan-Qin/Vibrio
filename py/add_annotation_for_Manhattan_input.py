import pandas as pd
import numpy as np
import itertools
import csv
import os

# Annotate input for Manhattan plot with gene names
drug = ['DOX', 'CIP', 'AZI']
oxygen_levels = ['A', 'AN']
plot_suffix = '_combined_to_EN1443_kmers.txt'
plot_all_suffix = '_combined_to_EN1443_kmers_all.txt'
plot_annotated_suffix = '_combined_to_EN1443_kmers_annotated.txt'
plot_all_annotated_suffix = '_combined_to_EN1443_kmers_all_annotated.txt'
hits_anno_suffix = '_SNP_hits_annotated.csv'
discrad_range = []

'''
for d in drug:
    plot1 = pd.read_table('manhattan_input/newAsem_snps/'+d+plot_suffix, sep='\t')
    ae_plot = plot1[plot1['oxygen'] == 'A']
    anae_plot = plot1[plot1['oxygen'] == 'AN']
    com_lst = list(itertools.product(oxygen_levels, [d]))
    prefix_lst = []
    for c in com_lst:
        prefix_lst.append(''.join(map(str, c)))
    #print(prefix_lst)

    for p in prefix_lst:
        tmp_df = pd.read_csv('SNP_hits/reads_from_newAsem_medakaTree_matched_wGene_AC10/'+
                             p+hits_anno_suffix)
        tmp_df['oxygen'] = p.replace(d, '')
        tmp_df = tmp_df.loc[:, ['snploc', 'oxygen', 'concentration', 'gene']].rename(columns={'snploc': 'minLOG10(P)'})

        if p.replace(d, '') == 'A':
            tmp_ae_df = pd.merge(ae_plot, tmp_df, how='left', on=['minLOG10(P)', 'oxygen', 'concentration'])
        elif p.replace(d, '') == 'AN':
            tmp_anae_df = pd.merge(anae_plot, tmp_df, how='left', on=['minLOG10(P)', 'oxygen', 'concentration'])

    combined_df = pd.concat([tmp_ae_df, tmp_anae_df], ignore_index=True)
    combined_df.to_csv('manhattan_input/newAsem_snps/'+d+plot_suffix.strip('.txt')+'_annotated.txt',
                       sep='\t', index=False)
'''

# Annotate kmer hits with corresponding SNP hits
'''
for d in drug:
    plot1 = pd.read_table('manhattan_input/newAsem_kmers/'+d+plot_suffix, sep='\t')
    #loc = plot1['minLOG10(P)'].values
    plot1_nm = plot1.columns.values.tolist()
    overlap_nm = plot1_nm + ['gene']
    plot1_lst = plot1.values.tolist()

    ae_lst = [i for i in plot1_lst if i[5] == 'A']
    anae_lst = [i for i in plot1_lst if i[5] == 'AN']

    com_lst = list(itertools.product(oxygen_levels, [d]))
    prefix_lst = []
    for c in com_lst:
        prefix_lst.append(''.join(map(str, c)))
    #print(prefix_lst)

    overlap_hits = []

    for p in prefix_lst:
        tmp_df = pd.read_csv('SNP_hits/reads_from_newAsem_medakaTree_matched_wGene_AC10/'+
                             p+hits_anno_suffix)

        tmp_df['oxygen'] = p.replace(d, '')
        tmp_df = tmp_df.loc[:, ['start', 'end', 'oxygen', 'concentration', 'gene']]

        tmp_nm = tmp_df.columns.values.tolist()
        tmp_lst = tmp_df.values.tolist()

        derep_lst = []
        for line in tmp_lst:
            if line not in derep_lst:
                derep_lst.append(line)
        #print(len(derep_lst))

        if p.replace(d, '') == 'A':
            tmp_hits = []
            for line1 in ae_lst:
                counter = 0

                if len(derep_lst) == 0:
                    line1.append('NA')
                    tmp_hits.append(line1)

                else:
                    for line2 in derep_lst:
                        if counter < len(derep_lst) - 1:
                            if int(line2[0]) <= int(line1[2]) <= int(line2[1]) and line1[5] == line2[2] and line1[6] == line2[3]:
                                line1.append(line2[4])
                                tmp_hits.append(line1)
                                break
                            else:
                                counter += 1

                        elif counter == len(derep_lst) - 1:
                            if int(line2[0]) <= int(line1[2]) <= int(line2[1]) and line1[5] == line2[2] and line1[6] == line2[3]:
                                line1.append(line2[4])
                                tmp_hits.append(line1)
                            else:
                                #print(line1)
                                line1.append('NA')
                                tmp_hits.append(line1)
                                break

            for line in tmp_hits:
                if len(line) == 10:
                    line.pop(-2)

            overlap_hits += tmp_hits


        elif p.replace(d, '') == 'AN':
            tmp_hits = []
            for line1 in anae_lst:
                counter = 0

                if len(derep_lst) == 0:
                    line1.append('NA')
                    tmp_hits.append(line1)

                else:
                    for line2 in derep_lst:
                        if counter < len(derep_lst) - 1:
                            if int(line2[0]) <= int(line1[2]) <= int(line2[1]) and line1[5] == line2[2] and line1[6] == \
                                    line2[3]:
                                line1.append(line2[4])
                                tmp_hits.append(line1)
                                break
                            else:
                                counter += 1

                        elif counter == len(derep_lst) - 1:
                            if int(line2[0]) <= int(line1[2]) <= int(line2[1]) and line1[5] == line2[2] and line1[6] == \
                                    line2[3]:
                                line1.append(line2[4])
                                tmp_hits.append(line1)
                            else:
                                # print(line1)
                                line1.append('NA')
                                tmp_hits.append(line1)
                                break

            for line in tmp_hits:
                if len(line) == 10:
                    line.pop(-2)

            overlap_hits += tmp_hits

    overlap_df = pd.DataFrame(overlap_hits, columns=overlap_nm)
    print(len(overlap_df))
    overlap_df.to_csv('manhattan_input/newAsem_kmers/annotated/'+d+'_combined_to_EN1443_kmers_annotated.txt',
                      sep='\t', index=False)
'''

# Annotate input for all kmers
for d in drug:
    plot1 = pd.read_table('manhattan_input/newAsem_kmers/annotated/'+d+plot_annotated_suffix, sep='\t')
    plot1.drop_duplicates(inplace=True, ignore_index=True)

    plot1_all = pd.read_table('manhattan_input/newAsem_kmers/'+d+plot_all_suffix, sep='\t')
    plot1_all.drop_duplicates(inplace=True, ignore_index=True)
    #print(len(plot1_all))

    col_nm_all = plot1_all.columns.values.tolist()
    overlap_plot = pd.merge(plot1_all, plot1, how='left', on=col_nm_all)
    #print(len(overlap_plot))

    derep_plot = overlap_plot[~overlap_plot.duplicated(subset=['oxygen', 'gene']) | (overlap_plot['gene'].isna())]
    #print(len(derep_plot))
    derep_plot.loc[(2912215 < derep_plot['minLOG10(P)']) & (derep_plot['minLOG10(P)'] < 3022483), 'gene'] = ''

    derep_plot.to_csv('manhattan_input/newAsem_kmers/annotated/'+d+plot_all_annotated_suffix,
                        sep='\t', index=False)
