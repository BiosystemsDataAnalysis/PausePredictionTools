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

# https://builtin.com/data-science/how-to-speed-up-pandas

# https://stackoverflow.com/questions/68547194/apply-numpy-where-along-one-of-axes



# add ORF 10% information
def addORFinfo(df):
    
    orf_range = 0.05 # percentage 
    df['gene_len']=np.abs(df['gene_stop']-df['gene_start'])
    df['gene_orf_start_tenperc']=(df.dir==0)*(df.gene_start+np.round(df.gene_len*orf_range,0))+(df.dir==16)*(df.gene_start-np.round(df.gene_len*orf_range,0))
    df['gene_orf_stop_tenperc']=(df.dir==0)*(df.gene_stop-np.round(df.gene_len*orf_range,0))+(df.dir==16)*(df.gene_stop+np.round(df.gene_len*orf_range,0))

    df['within_90_ORF']=((df.dir==0) & ((df.begin_read >= df.gene_orf_start_tenperc) & (df.end_read<= df.gene_orf_stop_tenperc))) | \
        ((df.dir==16) & ((df.begin_read <= df.gene_orf_start_tenperc) & (df.end_read>= df.gene_orf_stop_tenperc)))

    return df


def detpos(vecdata, offset35:int=12,offset53:int=14):
      
    #/offset53=_offset53
    #offset35=_offset35

    # print(offset53,offset35)

    strand = vecdata.iloc[0]
    start_read = vecdata.iloc[1]
    read_len = vecdata.iloc[2]
    position_I = vecdata.iloc[3]
    
    stop_read = start_read + (read_len-1)
    if strand==16:
        stop_read = start_read - (read_len - 1)  

    result_53 = [getAsiteData53(start_read,p,strand,offset53) for p in position_I]
    result_35 = [getAsiteData35(stop_read,p,strand,offset35) for p in position_I]
  
    return [result_53,result_35]


def add_aa_scores(dfin:pd.DataFrame,aa:str='I',offset35:int=12,offset53:int=14):
    
    # copy I data vectors and determine distances to reads and I position, in detpos the offsets are determined 
    # they can be set by the global variables _offset53 and _offset35
    ccc = dfin[['dir','begin_read','read_len',aa]].apply(lambda x:detpos(x),axis=1,result_type="expand")
        
    col35_, col53_= aa+"35_"+str(offset35), aa+"35_"+str(offset53)
    ccc.columns = [col53_,col35_]

    # find minimum distance between I and A-site 
    myvec = []
    for x in ccc[col53_].values.tolist():
        if(len(x)>0):
            myvec.append(x[np.argmin(np.abs(x))])
        else:        
            myvec.append(None)

    dfin["min"+ aa + "53"] = myvec

    myvec = []
    for x in ccc[col35_].values.tolist():
        if(len(x)>0):
            myvec.append(x[np.argmin(np.abs(x))])
        else:        
            myvec.append(None)

    dfin["min"+ aa + "35"] = myvec

    return dfin


def where(arr):
    # check for those items that are empty
    adj = np.sum(arr, axis=1) == 0
    cs = np.cumsum(arr,axis=1)
    x = np.argmax(cs, axis=1) - adj
    # set those items that have multiple genes to false
    x[cs.max(axis=1)>1]=-1
    return x

def process_sam_vectorized(blk_in, fasta_str, gns_info, index=0,offset35:int=12,offset53:int=14,verbose=False):
    
    try:
    
        gns_def = gns_info['DEF']
        ROISTART = gns_info['ROISTART']
        ROISTOP = gns_info['ROISTOP']
        GENESTART=gns_info['GENESTART']
        GENESTOP=gns_info['GENESTOP']
        GENESTRANDS=gns_info['GENESTRANDS']
        NP_AA=gns_info['NP_AA']
        AA_dict=gns_info['AA_DICT']

        def get_fasta_part(l,r):        
            _fs = fasta_str[l-1:r]        
            return _fs

        def get_revCompl(aStr,rev):        
            _fs = aStr
            if rev:
                _fs = revCompl(_fs)        
            return _fs

        print('processing sam data (block {0}) ...'.format(index+1))   
        
        if verbose:
            print("copy data ..")
        blk = blk_in.copy()
        gs = np.array(GENESTRANDS)

        if verbose:
            print("determine the (stranded) ends of the reads .. {0}".format(index+1))
        
        blk['read_len'] = blk.sam_str.str.len()
        blk['right_pos'] = blk.left_pos + blk.read_len - 1    

        if verbose:
            print('determining reverse complements data (block {0}) ...'.format(index+1))
        
        blk['read_str'] = list(map(get_revCompl, blk.sam_str,blk.dir))
        blk['fasta_match'] = list(map(get_fasta_part,blk.left_pos,blk.right_pos))
        
        blk['begin_read'] = (blk.dir==0)*blk.left_pos + (blk.dir==16)*blk.right_pos
        blk['end_read'] = (blk.dir==0)*blk.right_pos + (blk.dir==16)*blk.left_pos
        
        # blk.apply(
        #     lambda row: (row.begin_read + row.read_len - 1) if row.dir == 0 
        #     else (row.begin_read - (row.read_len - 1)), axis=1)

        f_start_0  = lambda row: (row.begin_read  >= GENESTART)
        f_end_0    = lambda row: (row.end_read    <= GENESTOP) 
        
        f_start_16 = lambda row: (row.end_read   >= GENESTART) 
        f_end_16   = lambda row: (row.begin_read <= GENESTOP) 
            
        if verbose:
            print("match regions of interest {0} .. ".format(index+1))   

        df_roi_0  = blk.apply(lambda row: (row.begin_read >= ROISTART) & ((row.end_read) <= ROISTOP),axis=1)
        df_roi_16 = blk.apply(lambda row: (ROISTART >= row.begin_read) & (ROISTOP <= (row.end_read)),axis=1)        

        np_roi_0 = where(np.array(df_roi_0.to_list()))
        np_roi_16 = where(np.array(df_roi_16.to_list()))
        # map ROI hits to their associated genes (indices)
        gns_roi_all = np.where(blk.dir==0,np_roi_0,np_roi_16)

        if verbose:
            print("match genes {0} ..".format(index+1))
        # find the start and stop of the genes
        np_start_0 = np.array(blk.apply(f_start_0, axis=1).to_list())
        np_end_0 = np.array(blk.apply(f_end_0, axis=1).to_list())
        np_start_16 = np.array(blk.apply(f_start_16, axis=1).to_list())
        np_end_16 = np.array(blk.apply(f_end_16, axis=1).to_list())
        # 
        rp = np.repeat(gs,blk.shape[0],axis=0).reshape(4178,blk.shape[0]).transpose()
        # find genes matches per strand
        gns_0 = where( (np_start_0 & np_end_0) & (rp == 0))
        gns_16 = where((np_start_16 & np_end_16) & (rp == 16))

        # select the genes per strand
        gns_all = np.where(blk.dir==0,gns_0,gns_16)
        
        # use gene definitions from ROI if no genes found
        gns_all = np.where(gns_all>=0,gns_all,gns_roi_all)
        
        gns_match = gns_all>=0
        
        # copy the gene strings to frame

        blk["gene"]=""
        blk.loc[gns_match,"gene"]=gns_def.iloc[gns_all[gns_match]].gene.values
        blk["gene_start"]=np.NaN
        blk.loc[gns_match,"gene_start"]=GENESTART[gns_all[gns_match]]
        blk["gene_stop"]=np.NaN
        blk.loc[gns_match,"gene_stop"]=GENESTOP[gns_all[gns_match]]
        
        
        def st2codon(row):
            if row.dir==0:
                _strt = row.gene_start
                _codon = (row.begin_read - _strt) // 3
                pos_from_start = row.begin_read - (_strt - 1)            
            else:
                _strt = row.gene_stop
                _codon = (_strt - row.begin_read) // 3
                pos_from_start = _strt -  row.begin_read
                
            return [pos_from_start,_codon]        

        blk["gene_pos"]=np.NaN
        blk['gene_codon']=np.NaN

        if verbose:
            print("determine position relative to start of gene {0} ..".format(index+1))    
        
        _aa =blk.loc[gns_match].apply(lambda x:st2codon(x),axis=1,result_type='expand')
        # check if there are any genes matched in this block (happens with small blocks)
        if _aa.shape[0]>0:
            blk.loc[gns_match,"gene_pos"]=_aa.iloc[:,0]
            blk.loc[gns_match,"gene_codon"]=_aa.iloc[:,1]
        
        # for full datablock
        blk["on_roi"]=gns_roi_all!=-1    

        # reverse start/stop 
        _tmp_start = np.where(blk.dir==0,blk.gene_start,blk.gene_stop)
        _tmp_stop = np.where(blk.dir==0,blk.gene_stop,blk.gene_start)
        
        blk["gene_start"]=_tmp_start.copy()
        blk["gene_stop"]=_tmp_stop.copy()

        if verbose:
            print('determining A site positions {0} ...'.format(index+1))
        
        blk['posAsite35']=blk.apply(lambda x:getAsiteData35(x.end_read,x.gene_start,x.dir,offset35),axis=1)
        blk['posAsite53']=blk.apply(lambda x:getAsiteData53(x.begin_read,x.gene_start,x.dir,offset53),axis=1)
            
        if verbose:
            print("augment with AA gene positions {0} ..".format(index+1))
    
        _arrAA = {}
        for ii,aa in zip(range(len(AA_dict.keys())),AA_dict.keys()):
            _arrAA[aa] = [NP_AA[g][ii] if g!=-1 else [] for g in gns_all]

        df_AA = pd.DataFrame.from_dict(_arrAA,orient='columns')
        # set the index otherwise empty rows
        df_AA.index = blk.index

        blk = pd.concat([blk,df_AA],axis=1)

        np_coverage_normal = np.zeros((len(fasta_str),))
        np_coverage_reverse = np.zeros((len(fasta_str),))
        
        def update_read_counts(strt:int,stp:int,strand:int):
            if strand==0:
                np_coverage_normal[range(strt-1,stp)]+=1
            if strand==16:
                np_coverage_reverse[range(strt-1,stp)]+=1

        if verbose:
            print('determining overall genome coverage {0} ..'.format(index+1))        

        blk.apply(lambda x:update_read_counts(x.left_pos,x.right_pos,x.dir),axis=1)

        blk = addORFinfo(blk)
        # find lastkey (for printing)
        # lstkey = list(AA_dict.keys())[-1]
        
        if verbose:
            print("adding distances for the different aminoacids {0} .. ".format(index+1))
        for k_ in AA_dict.keys():       
            # if verbose:
            #     if k_!=lstkey:
            #         print("{0}, ".format(k_),end="", flush=True)     
            #     else:
            #         print("{0}".format(k_))
            blk = add_aa_scores(blk,k_)
    

        return blk, np_coverage_normal, np_coverage_reverse
    
    except Exception as me:
        print("Memory error (block {0})".format(index), flush=True)
        return -1, [], []    



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
