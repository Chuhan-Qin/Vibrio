import csv

homedir = '/Users/chuhanqin/Desktop/my_pc/Vibrio/Bangladesh_assemblies/'
outdir = homedir+'R_files/plot_GWAS_gene_overlap/filtered_hits/'

group_lst = ['AAZI', 'ANAZI', 'ACIP', 'ANCIP', 'ADOX', 'ANDOX']


def extract_uniq_gene_lst(gene_lst_path):
    with open(gene_lst_path, 'r') as i_csv:
        reader = csv.reader(i_csv)
        unique_gene_lst = []
        for index, line in enumerate(reader):
            if index == 0:
                pass
            else:
                if '_' in line[1]:
                    if line[1].split('_')[0] not in unique_gene_lst:
                        unique_gene_lst.append(line[1].split('_')[0])
                else:
                    if line[1] not in unique_gene_lst:
                        unique_gene_lst.append(line[1])

    return unique_gene_lst


# Write comma-separated lists of unique gene hits
for g in group_lst:
    temp_lst = extract_uniq_gene_lst('filtered_gene_hits/'+g+'_gene_hits_EN1443_filtered.csv')
    with open(outdir+g+'_genes.txt', 'w') as o_txt:
        for ind, i in enumerate(temp_lst):
            if ind != len(temp_lst)-1:
                o_txt.writelines(i+',')
            else:
                o_txt.writelines(i)


