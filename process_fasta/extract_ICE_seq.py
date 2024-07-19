from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import csv
import pandas as pd
import os

# Extract ICE sequences from FASTAs of whole genomes
# 3' coordinate of ICE corresponds to the start of prfC
# 5' coordinate of ICE corresponds to the end of gmr

# Process csv file containing blast hits
work_dir = '/Users/chuhanqin/Desktop/my_pc/blast_searches/fasta_sets/BanP1P2P3_fna_files/'

'''
# Filter hits by alignment size
with open(work_dir+'prfC/prfC.csv', 'r') as i_csv:
    reader = csv.reader(i_csv)
    prfC_hits = [line for index, line in enumerate(reader) if index >= 1 and int(line[3]) >= 1000]
prfC_cods = [[line[1], min(line[8], line[9]), max(line[8], line[9])] for line in prfC_hits]
prfC_codDf = pd.DataFrame(prfC_cods, columns=['Isolate', 'prfC_s', 'prfC_e'])

with open(work_dir+'VchInd6_5endhit/VchInd6_5endhit.csv', 'r') as i_csv:
    reader = csv.reader(i_csv)
    vch5e_hits = [line for index, line in enumerate(reader) if index >= 1 and int(line[3]) >= 10000]
vch5e_cods = []+[[line[1], min(line[8], line[9]), max(line[8], line[9])] for line in vch5e_hits]
vch5e_codDf = pd.DataFrame(vch5e_cods, columns=['Isolate', 'vch5e_s', 'vch5e_e'])

with open(work_dir+'gmr/gmr.csv', 'r') as i_csv:
    reader = csv.reader(i_csv)
    gmr_hits = [line for index, line in enumerate(reader) if index >= 1]
gmr_cods = [[line[1], min(line[8], line[9]), max(line[8], line[9])] for line in gmr_hits]
gmr_codDf = pd.DataFrame(gmr_cods, columns=['Isolate', 'gmr_s', 'gmr_e'])

# Quality control - compare coordinates on isolates genomes
# Combine coordinates for extracting ICEs
ICE_cods = gmr_codDf.set_index('Isolate').join(vch5e_codDf.set_index('Isolate')).join(prfC_codDf.set_index('Isolate'))
#result = pd.merge(gmr_codDf, vch5e_codDf, on='Isolate')
print(ICE_cods)

ICE_cods.to_csv(work_dir+'ICE_coordinate_combined/ICE_coordinate_summarisation.csv')
'''

fasta_dir = '/Users/chuhanqin/Desktop/my_pc/Vibrio/Bangladesh_assemblies/BanP1P2P3_fna_files/'


def extract_seq_by_cod(id_cod_list):  # Input list in the format: ID, feature start coordinate, feature end coordinate
    chrs_dir = {}
    ices_dir = {}

    for line in id_cod_list:
        chrs_dir[line[0]] = ''
        ices_dir[line[0]] = ''
        for seq in SeqIO.parse(fasta_dir+line[0]+'.fasta', 'fasta'):
            chrs_dir[line[0]] += seq.seq
            ices_dir[line[0]] = SeqRecord(chrs_dir[line[0]][int(line[1]):int(line[2])], id=line[0]+'_ICE',
                                          description='SXT-R391 family of ICE')

    return ices_dir

# Read coordinates of markers
cods_matrix = []
with open(work_dir+'ICE_coordinate_combined/ICE_coordinate_summarisation.csv', 'r') as i_csv:
    reader = csv.reader(i_csv)
    for index, line in enumerate(reader):
        if index >= 1:
            cods_matrix.append(line)

# Running
ICE_range_list = []
for line in cods_matrix:
    if line[3] == '':
        pass
    else:
        if int(line[5])>int(line[1]):
            ICE_range_list.append([line[0].split('.')[0], line[2], line[5]])
        if int(line[5])<int(line[1]):
            ICE_range_list.append([line[0].split('.')[0], line[6], line[1]])

ice_dict = extract_seq_by_cod(ICE_range_list)
#print(ice_dict)

#chrs = {'1318700': ''}
#for seq in SeqIO.parse('../Comparison/SXT_R391/1318700.fasta', 'fasta'):
    #chrs['1318700'] += seq.seq
#print(chrs['1318700'][4028377:4029923])

# Write extracted ICE sequences
out_dir = '/Users/chuhanqin/Desktop/my_pc/Vibrio/Bangladesh_assemblies/SXT_ICE_fasta/'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

for k, v in ice_dict.items():
    SeqIO.write(v, out_dir+k+'_ICE.fasta', 'fasta')

