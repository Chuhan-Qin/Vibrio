import os
import csv
import gffutils
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# List all gene name and function annotations in the gff files
gffdir = '/Users/chuhanqin/Desktop/my_pc/Vibrio/Bangladesh_assemblies/gff_files_newAsem/'
fnadir = '/Users/chuhanqin/Desktop/my_pc/Vibrio/Bangladesh_assemblies/all_fna_reoriented/'
refgff = '/Users/chuhanqin/Desktop/my_pc/Vibrio/Bangladesh_assemblies/Reference/GCF_900205735.1.gff'
reffna = '/Users/chuhanqin/Desktop/my_pc/Vibrio/Bangladesh_assemblies/Reference/GCF_900205735.1.fna'
outdir = '/Users/chuhanqin/Desktop/my_pc/Vibrio/Bangladesh_assemblies/pyseer/continuous_AUC_phenotype/blastx_fix_annotation/'

refgff2 = '/Users/chuhanqin/Desktop/my_pc/Vibrio/Bangladesh_assemblies/gff_files_newAsem/EN_1443.gff'

# List of gene hits for annotating
genehitsdir = '/Users/chuhanqin/Desktop/my_pc/Vibrio/Bangladesh_assemblies/py_files/gwas_gene_hits/filtered_gene_hits_newAsem/'


def create_gene_func_ls(path_ref, path_else):
    # Tabulate gene names and corresponding functional annotations
    # Reference genome
    dict = {}
    with open(path_ref, 'r') as i_txt:
        for line in i_txt.readlines():
            if len(line.split('\t')) == 9 and 'gene' in line.split('\t')[8]:   # Skip hypothetical proteins
                for index, item in enumerate(line.split('\t')[8].split(';')):
                    if item.split('=')[0] == 'gene':
                        temp_gen_nm = item.split('=')[1]
                        dict[temp_gen_nm] = ''
                    if item.split('=')[0] == 'product':
                        dict[temp_gen_nm] = item.split('=')[1].strip('\n')
    #print(dict['nanH'])
    #print(len(gene_func_dict))

    # Iterate through all gff files of the isolates
    for root, dirs, files in os.walk(path_else):
        for file in files:
            with open(path_else+file, 'r') as i_txt:
                for line in i_txt.readlines():
                    if len(line.split('\t')) == 9:
                        fields = line.split('\t')[8].split(';')
                        for idx, field in enumerate(fields):
                            if field.startswith('gene='):
                                temp_key = field.split('=')[1]
                            elif field.startswith('product='):
                                temp_val = field.split('=')[1].strip('\n')
                            else:
                                pass
                        if temp_key not in dict.keys():
                            dict[temp_key] = temp_val
    return dict


def get_hit_function(group_name, gene_func_dict):
    for g in group_name:
        new_hits_matrix = []
        with open(genehitsdir+g+'_gene_hits_EN1443_filtered.csv', 'r') as i_csv:
            reader = csv.reader(i_csv)
            for index, line in enumerate(reader):
                if index == 0:
                    header = line[0:2]+['function']+line[2:6]+[line[6].strip('\n')]
                else:
                    if line[1] in gene_func_dict.keys():
                        new_line = line[0:2] + [gene_func_dict[line[1]]] + line[2:7]
                        new_hits_matrix.append(new_line)
                    else:
                        if '|' in line[1]:
                            for i in line[1].split('|'):
                                #print(i)
                                if i in gene_func_dict.keys():
                                    new_line = line[0:2] + [gene_func_dict[i]] + line[2:7]
                                    new_hits_matrix.append(new_line)
                        else:
                            new_line = line[0:2] + ['hypothetical'] + line[2:7]
                            new_hits_matrix.append(new_line)

        with open(genehitsdir+g+'_gene_hits_wFunc_EN1443_filtered.csv', 'w') as o_csv:
            writer = csv.writer(o_csv)
            writer.writerow(header)
            for line in new_hits_matrix:
                writer.writerow(line)


# Run - get function and coordinates for gene hits
# Get gene name function pairs
'''
gene_func_dict = create_gene_func_ls(refgff, gffdir)
print(gene_func_dict['umuC_1'])
analysis_group = ['AAZI', 'ACIP', 'ADOX', 'ANAZI', 'ANCIP', 'ANDOX']
get_hit_function(analysis_group, gene_func_dict)
'''

'''
# Get corresponding gff files and IDs
gff_id_pair_dict = {'BEDIDOIH': 'GCF_900205735.1.gff'}
for root, dirs, files in os.walk(gffdir):
    for file in files:
        with open(gffdir+file, 'r') as i_txt:
            for line in i_txt.readlines():
                if len(line.split('\t')) == 9:
                    for index, item in enumerate(line.split('\t')[8].split(';')):
                        if item.split('=')[0] == 'ID':
                            gff_id_pair_dict[item.split('=')[1].split('_')[0]] = file
                    break


# Extract sequences of unnamed genes using coordinates
def extract_unnamed_gene_seq(group_ls):
    for g in group_ls:
        seq_coord = []
        with open(genehitsdir+g+'_gene_hits_wFunc.csv', 'r') as i_csv:
            reader = csv.reader(i_csv)
            for index, line in enumerate(reader):
                if index == 0:
                    pass
                elif line[2] == 'hypothetical':
                    if line[1].split('_')[0] == 'BEDIDOIH':
                        with open(refgff, 'r') as i_gff:
                            for anno in i_gff.readlines():
                                if len(anno.split('\t')) == 9 and anno.split('\t')[8].split(';')[0].split('=')[1] == line[1]:
                                    seq_coord.append(['GCF_900205735.1.gff', line[1], anno.split('\t')[0],
                                                      anno.split('\t')[3], anno.split('\t')[4]])
                    else:
                        if line[1].split('_')[0] in gff_id_pair_dict.keys():
                            with open (gffdir+gff_id_pair_dict[line[1].split('_')[0]], 'r') as i_gff:
                                for anno in i_gff.readlines():
                                    if len(anno.split('\t')) == 9 and anno.split('\t')[8].split(';')[0].split('=')[1] == line[1]:
                                        seq_coord.append([gff_id_pair_dict[line[1].split('_')[0]], line[1], anno.split('\t')[0],
                                                          anno.split('\t')[3], anno.split('\t')[4]])
        #print(seq_coord, '\n', len(seq_coord))

        geneSeq_dict = {}
        for line in seq_coord:
            chrs_dict = {}
            geneSeq_dict[line[1]] = ''
            if line[0] == 'GCF_900205735.1.gff':
                for seq in SeqIO.parse(reffna, 'fasta'):
                    chrs_dict[seq.id] = seq.seq
                try:
                    geneSeq_dict[line[1]] = SeqRecord(chrs_dict[line[2]][int(line[3]):int(line[4])], id=line[1])
                except KeyError:
                    print(line[1])
            else:
                for seq in SeqIO.parse(fnadir+line[0].strip('.gff')+'.fasta', 'fasta'):
                    chrs_dict[seq.id] = seq.seq
                try:
                    geneSeq_dict[line[1]] = SeqRecord(chrs_dict[line[2]][int(line[3]):int(line[4])], id=line[1])
                except KeyError:
                    print(line[1])
        #print(len(geneSeq_dict['BEDIDOIH_01742'].seq))

        # Write extracted gene sequences
        if not os.path.exists(outdir+g):
            os.mkdir(outdir+g)

        for k, v in geneSeq_dict.items():
            SeqIO.write(v, outdir + g + '/'+k+'.fasta', 'fasta')
'''

# Run function for extracting gene sequences
#extract_unnamed_gene_seq(analysis_group)

# Test
'''
seq_coord = []
with open(genehitsdir+'AAZI_gene_hits_wFunc.csv', 'r') as i_csv:
    reader = csv.reader(i_csv)
    for index, line in enumerate(reader):
        if index == 0:
            pass
        elif line[2] == 'hypothetical':
            if line[1].split('_')[0] == 'BEDIDOIH':
                with open(refgff, 'r') as i_gff:
                    for anno in i_gff.readlines():
                        if len(anno.split('\t')) == 9 and anno.split('\t')[8].split(';')[0].split('=')[1] == line[1]:
                            seq_coord.append(['GCF_900205735.1.gff', line[1], anno.split('\t')[0],
                                              anno.split('\t')[3], anno.split('\t')[4]])
            else:
                if line[1].split('_')[0] in gff_id_pair_dict.keys():
                    with open (gffdir+gff_id_pair_dict[line[1].split('_')[0]], 'r') as i_gff:
                        for anno in i_gff.readlines():
                            if len(anno.split('\t')) == 9 and anno.split('\t')[8].split(';')[0].split('=')[1] == line[1]:
                                seq_coord.append([gff_id_pair_dict[line[1].split('_')[0]], line[1], anno.split('\t')[0],
                                                    anno.split('\t')[3], anno.split('\t')[4]])
print(seq_coord, '\n', len(seq_coord))

chrs_dict = {}
geneSeq_dict = {}
for line in seq_coord:
    chrs_dict[line[1]] = ''
    geneSeq_dict[line[1]] = ''
    if line[0] == 'GCF_900205735.1.gff':
        #print(line[1])
        for seq in SeqIO.parse(reffna, 'fasta'):
            chrs_dict[seq.id] = seq.seq
            geneSeq_dict[line[1]] = SeqRecord(chrs_dict[line[2]][int(line[3]):int(line[4])], id=line[1])
    else:
        #print(line[1])
        for seq in SeqIO.parse(fnadir+line[0].strip('.gff')+'.fasta', 'fasta'):
            chrs_dict[seq.id] = seq.seq
            geneSeq_dict[line[1]] = SeqRecord(chrs_dict[line[2]][int(line[3]):int(line[4])], id=line[1])
print(geneSeq_dict['BEDIDOIH_03266'].seq)
'''
'''
# Run blastX on extracted unknown genes
for root, dirs, files in os.walk(outdir+'AAZI/'):
    for file in files:
        os.system('conda run -n blast_env blastx -db nr -query '+outdir+'AAZI/'+file+' -out '+
                  file+'.blastx.txt -remote')
'''


# Create reference library for annotation
def gffuntils_create_func_lib (path_ref, dir_else):
    ref_func_dict = {}

    ref_fn = gffutils.example_filename(path_ref)
    ref_db = gffutils.create_db(ref_fn, dbfn=':memory:')
    for feature in ref_db.all_features():
        if 'gene' in feature.attributes.keys():
            if feature.attributes['gene'][0] not in ref_func_dict.keys():
                ref_func_dict[feature.attributes['gene'][0]] = feature.attributes['product'][0]
        else:
            if feature.attributes['ID'][0] not in ref_func_dict.keys():
                ref_func_dict[feature.attributes['ID'][0]] = feature.attributes['product'][0]

    for root, dirs, files in os.walk(dir_else):
        for file in files:
            temp_fn = gffutils.example_filename(root+file)
            temp_db = gffutils.create_db(temp_fn, dbfn=':memory:')

            for feature in temp_db.all_features():
                if 'gene' in feature.attributes.keys():
                    if feature.attributes['gene'][0] not in ref_func_dict.keys():
                        ref_func_dict[feature.attributes['gene'][0]] = feature.attributes['product'][0]
                else:
                    if feature.attributes['ID'][0] not in ref_func_dict.keys():
                        ref_func_dict[feature.attributes['ID'][0]] = feature.attributes['product'][0]

    return ref_func_dict


# Get positions of significant hits on the refseq
def gffutils_create_loc_lab(path_ref, chrom1, chrom2, chrom1_len):
    ref_pos_dict = {}

    ref_fn = gffutils.example_filename(path_ref)
    ref_db = gffutils.create_db(ref_fn, dbfn=':memory:')

    for feature in ref_db.all_features():
        if 'gene' in feature.attributes.keys():
            if feature.attributes['gene'][0] not in ref_pos_dict.keys():
                if feature.seqid == chrom1:
                    ref_pos_dict[feature.attributes['gene'][0]] = [feature.start, feature.end]
                if feature.seqid == chrom2:
                    ref_pos_dict[feature.attributes['gene'][0]] = [feature.start+chrom1_len, feature.end+chrom1_len]
        else:
            if feature.attributes['ID'][0] not in ref_pos_dict.keys():
                if feature.seqid == chrom1:
                    ref_pos_dict[feature.attributes['ID'][0]] = [feature.start, feature.end]
                if feature.seqid == chrom2:
                    ref_pos_dict[feature.attributes['ID'][0]] = [feature.start+chrom1_len, feature.end+chrom1_len]

    return ref_pos_dict


# Filter and split gene hits table
'''
p_cutoff = 5.37
for root, dirs, files in os.walk('SNP_hits/reads_from_newAsem_medakaTree_all_SNPs/'):
    for file in files:
        drug = file.split('_')[0]
        snp_hit_df = pd.read_table('SNP_hits/reads_from_newAsem_medakaTree_all_SNPs/'+file, sep='\t')

        asnp_hit_df = snp_hit_df.loc[(snp_hit_df['log10(p)'] > p_cutoff) & (snp_hit_df['oxygen'] == 'A'),
                                    ['minLOG10(P)', 'log10(p)', 'concentration']]

        ansnp_hit_df = snp_hit_df.loc[(snp_hit_df['log10(p)'] > p_cutoff) & (snp_hit_df['oxygen'] == 'AN'),
                                     ['minLOG10(P)', 'log10(p)', 'concentration']]

        asnp_hit_df = asnp_hit_df.rename(columns={'minLOG10(P)': 'snploc'})
        ansnp_hit_df = ansnp_hit_df.rename(columns={'minLOG10(P)': 'snploc'})

        asnp_hit_df.to_csv('SNP_hits/reads_from_newAsem_medakaTree_filtered_SNPs/'+'A'+drug+'_SNP_hits.csv', index=False)
        ansnp_hit_df.to_csv('SNP_hits/reads_from_newAsem_medakaTree_filtered_SNPs/'+'AN'+drug+'_SNP_hits.csv', index=False)
'''

# Run
'''
# Create reference library for products and positions
pd.options.mode.chained_assignment = None

ref_func_dict = gffuntils_create_func_lib(refgff, gffdir)
func_df = pd.DataFrame(ref_func_dict.items(), columns=['gene', 'product'])

ref_pos_dict = gffutils_create_loc_lab(refgff2, 'contig_1', 'contig_2', 3081133)

ref_pos_lst = []
for k, v in ref_pos_dict.items():
    ref_pos_lst.append([k, v[0], v[1]])

pos_df = pd.DataFrame(ref_pos_lst, columns=['gene', 'start', 'end'])

# Process GWAS hits results
for root, dirs, files in os.walk('all_gene_hits/new_asem/'):
    for file in files:
        temp_df = pd.read_csv('all_gene_hits/new_asem/'+file)

        for i, v in temp_df['gene'].iteritems():
            if '|' in v:
                temp_df['gene'][i] = v.split('|')[0]

        combined_df1 = pd.merge(temp_df, func_df, how='left', on='gene')
        combined_df2 = pd.merge(combined_df1, pos_df, how='left', on='gene')

        combined_df2.to_csv('all_gene_hits/new_asem_wFunc_loc/'+file, index=False)
'''

# Match SNP hits with kmer hits by coordinates
analysing_group = ['AAZI', 'ANAZI', 'ACIP', 'ANCIP', 'ADOX', 'ANDOX']
out_path = 'SNP_hits/reads_from_newAsem_medakaTree_matched_wGene_AC10/'
o_suffix = '_SNP_hits_annotated.csv'

for i in analysing_group:
    snp_df = pd.read_csv('SNP_hits/reads_from_newAsem_medakaTree_filtered_SNPs/'+i+'_SNP_hits.csv')
    snp_header = snp_df.columns.tolist()
    snp_lst = snp_df.values.tolist()

    gene_df = pd.read_csv('all_gene_hits/new_asem_wFunc_loc/'+i+'_gene_hits_EN1443.csv').sort_values(by=['start'])
    gene_df = gene_df[(~gene_df['start'].isnull()) & (gene_df['hits'] > 10)]
    gene_header = gene_df.columns.tolist()
    gene_lst = gene_df.values.tolist()

    combined_header = snp_header + [i.strip('\n') for i in gene_header[1:]]

    matched_hit_lst = []
    for snp in snp_lst:
        for gene in gene_lst:
            if snp[0] < gene[8]:
                pass
            elif snp[0] < gene[9]:
                matched_hit_lst.append(snp+gene[1:])
            else:
                pass

    if not os.path.exists(out_path):
        os.mkdir(out_path)

    with open(out_path+i+o_suffix, 'w') as o_csv:
        writer = csv.writer(o_csv)
        writer.writerow(combined_header)

        for line in matched_hit_lst:
            writer.writerow(line)

















