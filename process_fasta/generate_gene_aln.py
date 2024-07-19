from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os


def extract_gene_alignment(isolate_lst, input_aln, gene_start, gene_end, gene_nm, extra_len):
    seq_dict = {}
    plus_seq_dict = {}

    for record in SeqIO.parse(work_dir+input_aln, 'fasta'):
        if record.id in isolate_lst:
            seq_dict[record.id] = ''
            plus_seq_dict[record.id] = ''
            seq_dict[record.id] = SeqRecord(record.seq[gene_start-1:gene_end].upper(),
                                            id=record.id,
                                            description=gene_nm)

            plus_seq_dict[record.id] = SeqRecord(record.seq[gene_start-1-extra_len:gene_end+extra_len].upper(),
                                                      id=record.id,
                                                      description=gene_nm+' plus '+str(extra_len)+\
                                                                  ' bps up- and down-stream')

    return seq_dict, plus_seq_dict


work_dir = '/Users/chuhanqin/Desktop/my_pc/Vibrio/Bangladesh_assemblies/py_files/process_fasta/'
id_lst = ['1383600', 'D02182158',
          'D02182361', 'D15180683',
          '1385750', '1296600',
          'D02182360', 'D10182089',
          'D02182354', '1283200']

# Store rpsG seqs to a dictionary
rpsg_seq_dict, rpsg_plus_seq_dict = extract_gene_alignment(id_lst, 'genome_alignment/all_newAsem_for_GWAS_EN1443.aln',
                                                           2578366, 2578836, 'rpsG', 2000)

#for k, v in rpsg_seq_dict.items():
    #print(k, v[79])

# Write alignment
#SeqIO.write(rpsg_seq_dict.values(), 'output/rpsg/rpsg.aln', 'fasta')
#SeqIO.write(rpsg_plus_seq_dict.values(), 'output/rpsg/rpsg_plus_flanking.aln', 'fasta')

# Write sequences of the refseq
refseq_rpsg_dict = {}

for record in SeqIO.parse('EN_1443.fasta', 'fasta'):
    if record.id == 'contig_1':
        refseq_rpsg_dict['refseq'] = SeqRecord(record.seq[2576365:2580836].upper(), id='rpsG plus 2000bp flanking')
        break

#print(refseq_rpsg_dict.values())

SeqIO.write(refseq_rpsg_dict.values(), 'refseq_rpsg_plus_flanking.fasta', 'fasta')

'''
# Extract nanH [1933330:1935676]
# plus up- down-stream [1931330:1937676]
nanh_seq_dict = {}
nanh_plus_seq_dict = {}

for record in SeqIO.parse('./all_mapped_to_N16961_wB1B2.aln', 'fasta'):
    #print(record.seq[1935115])

    if record.id in id_lst:
        nanh_seq_dict[record.id] = ''
        nanh_plus_seq_dict[record.id] = ''
        nanh_seq_dict[record.id] = SeqRecord(record.seq[1933330:1935676].upper(), id=record.id, description='nanH')
        nanh_plus_seq_dict[record.id] = SeqRecord(record.seq[1931330:1937676].upper(), id=record.id,
                                             description='nanH plus 2000 bps up- and down-stream')

#print(nanh_seq_dict)
for k, v in nanh_plus_seq_dict.items():
    print(k, v.seq[3785])

#SeqIO.write(nanh_seq_dict.values(), './nanh.aln', 'fasta')
#SeqIO.write(nanh_plus_seq_dict.values(), './nanh_plus_flanking.aln', 'fasta')

#for record in SeqIO.parse('./nanh_pair_oldSeq.aln', 'fasta'):
    #print(record.id, record.seq[1935115])
'''
'''
file_lst = []
for root, dirs, files in os.walk('../../BanP1P2P3_fna_files'):
    for file in files:
        if 'fasta' in file and 'fasta.' not in file and 'EN' not in file:
            file_lst.append(file)
            #os.system('cp '+root+'/'+file+' ../../all_fna_wOldSeq/')

print(file_lst, '\n', len(file_lst))

for root, dirs, files in os.walk('../../all_fna_reoriented'):
    for file in files:
        if 'fasta' in file and 'fasta.' not in file and not file.startswith('DRR') and not file.startswith('SRR') and not file.startswith('EN') and file not in file_lst:
            os.system('cp ' + root + '/' + file + ' ../../all_fna_wOldSeq/')
'''
