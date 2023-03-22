#from distutils.log import warn
#from itertools import product
import re

#import matplotlib.pyplot as plt

import pandas as pd
import numpy as np 

import re
#from collections import Counter
#from enum import Enum

from Bio.Seq import Seq
#from Bio import SeqIO
from itertools import repeat
import warnings


def revCompl(seq):
    _bseq = Seq(seq)
    _rev_comp = _bseq.reverse_complement()
    # complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    # reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))

    return str(_rev_comp)

# copied from Mohammed 2019

def get_genetic_code():
    
    aa_keys = {
        'I' : ['ATA', 'ATC', 'ATT'],
        'M' : ['ATG'],
        'T' : ['ACA', 'ACC', 'ACG', 'ACT'],
        'N' : ['AAC', 'AAT'],
        'K' : ['AAA', 'AAG'],
        'S' : ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
        'R' : ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
        'L' : ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
        'P' : ['CCA', 'CCC', 'CCG', 'CCT'],
        'H' : ['CAC', 'CAT'],
        'Q' : ['CAA', 'CAG'],
        'V' : ['GTA', 'GTC', 'GTG', 'GTT'],
        'A' : ['GCA', 'GCC', 'GCG', 'GCT'],
        'D' : ['GAC', 'GAT'],
        'E' : ['GAA', 'GAG'],
        'G' : ['GGA', 'GGC', 'GGG', 'GGT'],
        'F' : ['TTC', 'TTT'],
        'Y' : ['TAC', 'TAT'],
        'C' : ['TGC', 'TGT'],
        'W' : ['TGG'],
        '_' : ['TAA', 'TAG', 'TGA']
        }
        
        
        
    codon_keys = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }
    
    return aa_keys, codon_keys


## function to load gene definition file from gff
def create_gene_definition_file(gene_defintion_filename, nr_headers=7):

    # % read gene definitions
    gene_def = pd.read_csv(gene_defintion_filename, sep='\t', skiprows=nr_headers, header=None)
    columns = ['id', 'or', 'tp', 'start', 'stop', 'nn1', 'strand', 'nn3', 'info']
    gene_def.columns = columns
    # remove last row
    gene_def.drop(gene_def.index[-1], inplace=True)

    # %
    def gname(_row):
        # print(_row)
        try:
            lt = _row.split('locus_tag=')[1].split(';')[0]
        except:
            lt = ''
        try:
            gn = _row.split('gene=')[1].split(';')[0]
        except:
            gn = ''

        return [gn, lt]

    # do not select the whole gene or gene regions.. they overlap completely with the CDS definitions
    # df_tmp = gene_def.loc[(gene_def.tp != 'region') & (gene_def.tp != 'gene') & (gene_def.tp !='sequence_feature'), ['start', 'stop','strand','info']].copy()
    df_tmp = gene_def.loc[(gene_def.tp == 'CDS'), ['start', 'stop','strand','info']].copy()
    # df_tmp = gene_def.loc[(gene_def.tp == 'gene'), ['start', 'stop','strand','info']].copy()
    df_tmp.index = range(df_tmp.shape[0])

    A = df_tmp['info'].apply(lambda x: pd.Series(gname(x)))
    A.columns = ['gene', 'locus_tag']
    df_tmp = df_tmp.merge(A, left_index=True, right_index=True)
    # because we want to map on genes -> remove the definitions that do not have such a tag
    df_genes = df_tmp[df_tmp.gene.apply(lambda x: len(x)) > 0].copy()

    #
    df_genes = df_genes.drop(['info'], axis=1).copy()
    # start and stop are always on the reference genome .. i.e. for all i start(i) < stop(i)
    
    df_genes.start = df_genes.start.astype('int32')
    df_genes.stop = df_genes.stop.astype('int32')
    df_genes.index = range(df_genes.shape[0])
    df_genes['len'] = df_genes.stop - df_genes.start + 1 # +1 added on 29/06/2022
    

    return df_genes

# this routine finds all occurences of the pattern.. also overlapping ... so 
# "CC" in "CCCC" returns 3 counts
def CntSubstr(pattern, string):
    a = [m.start() for m in re.finditer(
        '(?={0})'.format(re.escape(pattern)), string)]
    return a

def locate_amino_acids():
    return None        



def augment_gff_file(f_gff_file,nr_headers=7,rng=30,funGetFastaStr=None):
    ''' adjust the ORF regions by rng (default 30)
    '''
    gns = create_gene_definition_file(f_gff_file,nr_headers=nr_headers).copy()
    # print the last rows of the gns dataframe
    gns.tail()
    gns['ROI_START']=(gns.strand=='+')*(gns['start']-rng) + (gns.strand=='-')*(gns['stop']+rng)
    gns['ROI_STOP']=(gns.strand=='+')*(gns['start']+rng) + (gns.strand=='-')*(gns['stop']-rng)

    _aa_keys,_=get_genetic_code()

    columns2add = [_cd for _cd in _aa_keys.keys()]
    # add columns
    for _aa in columns2add:
        gns[_aa]=[[] for i in repeat(None, gns.shape[0])]
        gns[_aa]=gns[_aa].astype('object')

    for _gene in gns.index:
        _start_stop_strand = gns.loc[_gene,['start','stop','strand','gene']]

        _str = funGetFastaStr(_start_stop_strand.start, _start_stop_strand.stop)
        if _start_stop_strand.strand=='-':
            _str = revCompl(_str)

        # make dictionary of occuring aminoacids on gene
        _seq = Seq(_str)
        
        with warnings.catch_warnings(record=True) as w:            
            #warnings.simplefilter("always")
            _aa_seq = _seq.translate()
            # replace the stop codon definition from * to _
            _aa_seq = _aa_seq.replace('*','_')
            if len(w)==1:            
                print('{0} {1}'.format(_start_stop_strand.gene,w[-1].message))
                    
        _aa_cnt = 0
        _aa_pos = {}

        for _aa in _aa_seq:
            if _aa in _aa_pos:
                _aa_pos[_aa].append(_aa_cnt)
            else:
                _aa_pos[_aa]=[_aa_cnt]            
            _aa_cnt+=1        
        
        for _aa in _aa_pos.keys():    
            # make them nucleotide positions, so multiply by 3                                      
            if _start_stop_strand.strand == '+':
                gns.at[_gene,_aa]=[_start_stop_strand.start+_p*3  for _p in _aa_pos[_aa]]
            else: 
                gns.at[_gene,_aa]=[_start_stop_strand.stop-_p*3 for _p in _aa_pos[_aa]]
                      
    return gns



# %%
# make sure that your gene definitions are loaded first.. cannot pass it as an argument because of map call in process_sam_data

def getAsiteData53_V0(readStart:int,geneStart:int,strand:int,offset=14):
    '''
    returns the A site position determined from the 5'->3' direction with an offset of 15 nt
    readStart: position of the start of the read strand independent, i.e. if strand == -, start > end
    '''
    
    if strand==16: # reverse read
        return (readStart-offset)-geneStart

    return geneStart - (readStart+offset) # update 17/7/2022


def getAsiteData35_V0(readStop:int,geneStart:int,strand,offset=12):
    '''
    returns the A site position determined from the 3'->5' direction with a default offset of (-)12 nt
    readStop: position of the end of the read (=strand independent, i.e. if strand == -, end < start)        
    '''
    # 385, 410, 0, 9 -> (385-410)-15 -> -40
    # 385, 410, 16,9 -> (410-385)+15 -> +4

    if strand==16: # reverse read
        return (readStop+offset)-geneStart 

    return geneStart - (readStop-offset)
    
def getAsiteData53_V5(poi:int,refPoint:int,strand:int,offset=14):    
    if strand==16: # reverse read
        return -((poi-offset)-refPoint)
    return (poi + offset) - refPoint # update 19/7/2022

def getAsiteData35_V5(poi:int,refPoint:int,strand:int,offset=12):    
    if strand==16: # reverse read
        return -((poi+offset)-refPoint)
    return (poi-offset)-refPoint # update 19/7/2022


getAsiteData53 = getAsiteData53_V5
getAsiteData35 = getAsiteData35_V5


def get_gene_AAs(agene, aAA):
    global gns_def
    gene = gns_def[agene]        
    gene_aas = gene[aAA]
    return gene_aas.values[0]


def process_sam_data(df_sam, fasta_str, gns_def, index=0,offset35:int=12,offset53:int=14):    
    '''
    process the sam file to determine the A-site positions, relative to start of ORF
    '''

    def ongene(readPosition,strand,readLength):
        '''
        return gene if match, its position from start, its codon position from start,
        readPosition is strand sensitive (i.e. on reverse strand it start position > stop position)

        '''

        # readPosition is 1-based indexed
        
        _read_begin = readPosition 
        _read_end = readPosition + readLength - 1

        if strand==16:
            _read_end = readPosition - (readLength - 1)

        str_map = {16:'-',0:'+'}
        _on_roi = False

        if strand==0: # if on normal strand
            on_roi = (_read_begin >= gns_def.ROI_START) & (_read_end <= gns_def.ROI_STOP)
            _gns_start_chk = gns_def.start<=(readPosition)
            _gns_stop_chk = (_read_end)<=gns_def.stop
        if strand==16:
            on_roi = (gns_def.ROI_START>=_read_begin) & (gns_def.ROI_STOP <= _read_end)            
            _gns_start_chk = gns_def.start<=(_read_end)
            _gns_stop_chk = gns_def.stop>=(readPosition)

        # gns_def.start & stop are not strand sensitive.. i.e. stop > start for all genes 

        # _gene = (gns_def.start<=(readPosition-1)) & (gns_def.stop>(readPosition+readLength)) & (gns_def.strand==str_map[strand])    
        _gene = _gns_start_chk & _gns_stop_chk & (gns_def.strand==str_map[strand])    
        
        if (sum(on_roi)==1) | (sum(_gene)==1):

            if sum(_gene)>1:
                # if the hit count is not unique for a single gene then continue
                # return ["",np.nan,np.nan,np.nan,np.nan,False,[]]
                return ["",np.nan,np.nan,np.nan,np.nan,False]


            _gene_roi = on_roi & (gns_def.strand==str_map[strand])
            #(gns_def.start<=(aposition-1)) & (gns_def.stop>aposition) & (gns_def.strand==str_map[strand])    
            if(sum(_gene)==0):
                gene = gns_def[_gene_roi]
            else:
                gene = gns_def[_gene]        
            

            if strand==0:            
                _strt = gene.start.values[0]
                _codon = (readPosition-_strt)//3
                pos_from_start = readPosition- (_strt-1)                     
                _gstart = gene.start.values[0]
                _gstop = gene.stop.values[0]

            else: # strand == 16, reversed
                _strt = gene.stop.values[0]
                _codon = (_strt-readPosition)//3
                pos_from_start = _strt - readPosition
                _gstop = gene.start.values[0]
                _gstart = gene.stop.values[0]
        
            return [gene.iloc[0].gene,_gstart,_gstop,pos_from_start, _codon,sum(on_roi)==1]
        
        return ["",np.nan,np.nan,np.nan,np.nan,False]

    
    # declare local variable
    aa_to_extract = 'I'   
    
    def get_fasta_part(l,r):        
        _fs = fasta_str[l-1:r]        
        return _fs

    def get_revCompl(aStr,rev):        
        _fs = aStr
        if rev:
            _fs = revCompl(_fs)        
        return _fs

    
    print('processing sam data (block {0}) ...'.format(index+1))   
    df_select = df_sam.copy()
    df_select['read_len'] = df_select.sam_str.str.len()
    df_select['right_pos'] = df_select.left_pos + df_select.read_len - 1    

    print('determining reverse complements data (block {0}) ...'.format(index+1))
    
    df_select['read_str'] = list(map(get_revCompl, df_select.sam_str,df_select.dir))
    df_select['fasta_match'] = list(map(get_fasta_part,df_select.left_pos,df_select.right_pos))
    
    df_select['begin_read'] = (df_select.dir==0)*df_select.left_pos + (df_select.dir==16)*df_select.right_pos
    df_select['end_read'] = (df_select.dir==0)*df_select.right_pos + (df_select.dir==16)*df_select.left_pos
    
    print('determining hits in area of ORF (block {0}) ...'.format(index+1))

    df_select['gene'],df_select['gene_start'],df_select['gene_stop'],df_select['gene_pos'], \
        df_select['gene_codon'],df_select['on_roi'] = \
            zip(*map(ongene,df_select.begin_read,df_select.dir,df_select.read_len))
    
    print('determining A site positions (block {0}) ...'.format(index+1))
    
    df_select['posAsite35']=df_select.apply(lambda x:getAsiteData35(x.end_read,x.gene_start,x.dir,offset35),axis=1)
    #list(map(getAsiteData35,df_select.end_read,df_select.gene_start,df_select.dir))
    df_select['posAsite53']=df_select.apply(lambda x:getAsiteData53(x.begin_read,x.gene_start,x.dir,offset53),axis=1)
    #list(map(getAsiteData53,df_select.begin_read,df_select.gene_start,df_select.dir))         
    
    cols_to_order = ['read_str', 'sam_str','fasta_match']
    new_columns =  (df_select.columns.drop(cols_to_order).tolist())+cols_to_order 
    df_select = df_select[new_columns]

    # here we add the columns containing the positions for the individual amino acids (from the gene definition)
    svec={0:'+',16:'-'}
    def get_gene_AAs(agene, strand, aatoe):
        #print(aatoe)
        if len(agene)>0:        
            gene = gns_def[(gns_def.gene==agene) & (gns_def.strand==svec[strand])]        
            #gene_aas = gene[aa_to_extract]
            gene_aas = gene[aatoe]
            return gene_aas.values[0]
        return []

    AA_dict, _ = get_genetic_code()
    
    print('create A site dictionary (block {0}) ...'.format(index+1))

    for aa in AA_dict.keys():        
        df_select[aa]=df_select.apply(lambda x:get_gene_AAs(x.gene,x.dir,aa),axis=1)


    df_reads_normal = np.zeros((len(fasta_str),))
    df_reads_reverse = np.zeros((len(fasta_str),))
    
    def update_read_counts(strt,stp,strand):
        if strand==0:
            df_reads_normal[range(strt-1,stp)]+=1
        if strand==16:
            df_reads_reverse[range(strt-1,stp)]+=1

    print('summarizing counts per amino acid (block {0}) ...'.format(index+1))

    df_select.apply(lambda x:update_read_counts(x.left_pos,x.right_pos,x.dir),axis=1);


    return df_select,df_reads_normal,df_reads_reverse
# %%
