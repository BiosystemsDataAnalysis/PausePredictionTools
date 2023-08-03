#%%
from fileinput import filename
from os.path import exists
import sys

# from sys import platform
interactive_mode = hasattr(sys, 'ps1')

from io import StringIO

from hashlib import md5
from time import localtime
# check for missing packages

import subprocess
import pkg_resources

required = {'biopython', 'pandas','typed-argument-parser', 'numpy'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    python = sys.executable
    subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)

from Bio.Seq import Seq
import numpy as np
import pandas as pd
from tap import Tap # typed-argument-parser

# import psutil
import tempfile
import os

#%% check for arguments and installed packages
if interactive_mode:
    sys.argv=['file.py','--input_file','Mup_FK.sam','--fasta_file','WT-Prspb-amyI.fa']


class myargs(Tap):
    
    input_file: str  # the input file (only sam files).
    fasta_file: str # the fasta file to generate sequences for motif detection.
    gff_file: str="" # optional, if supplied the peaks are annnotated with the location on the genome e.g. gene, ORF of gene
    nt_upstream: int = 50 # the number of positions upstream (towards 5') from the first nt of the A-site .
    nt_downstream:int = 48 # the number of positions downstream (towards 3') from the first nt of the A-site.
    min_hits:int = 1 # pauses should have at least min_hits reads.
    down_offset:int = 12 # the offset when read from the 3' end.
    up_offset:int = 14 # the offset when read from the 5' end.
    min_qual:int = 42 # the minimum mapping quality of the reads.
    nlines:int = 1000000 # the number of lines to be processed in a single block.
    

args = myargs().parse_args()

def check_arguments():
    result = True
    
    if not exists(args.input_file):    
        print('The input file cannot be found')
        result = False

    if not exists(args.fasta_file):        
        print('The fasta file cannot be found')
        result = False

    if len(args.gff_file)>0 and (not exists(args.gff_file)):
        print('The gff file cannot be found')
        result = False    

    return result

if not check_arguments():
    exit(1)


        
## function to load gene definition file from gff
def create_gene_definition_file(gene_defintion_filename):

    f_ = open(gene_defintion_filename,'r')
    headers = []
    for l in f_:
        if l.startswith("#"):
            headers.append(l)
        else:
            break
    f_.close()

    # % read gene definitions
    gene_def = pd.read_csv(gene_defintion_filename, sep='\t', skiprows=len(headers), header=None)
    columns = ['id', 'or', 'tp', 'start', 'stop', 'nn1', 'strand', 'nn3', 'info']
    gene_def.columns = columns
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
    # df_tmp = gene_def.loc[(gene_def.tp == 'CDS'), ['start', 'stop','strand','info']].copy()
    # NB! This could be different for different GFF files.. perhaps we need some optional settings here
    
    df_tmp = gene_def.loc[(gene_def.tp == 'gene'), ['start', 'stop','strand','info']].copy()
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
    df_genes['len'] = df_genes.stop - df_genes.start
    

    orf_range = 0.05 # percentage 
    # df['gene_len']=np.abs(df['gene_stop']-df['gene_start'])
    df_genes['gene_orf_start_tenperc']=(df_genes.strand=='+')*(df_genes.start+np.round(df_genes.len*orf_range,0))+\
        (df_genes.strand=='-')*(df_genes.start-np.round(df_genes.len*orf_range,0))  
    df_genes['gene_orf_stop_tenperc']=(df_genes.strand=='+')*(df_genes.stop-np.round(df_genes.len*orf_range,0))+\
        (df_genes.strand=='-')*(df_genes.stop+np.round(df_genes.len*orf_range,0))


    return df_genes


def get_fasta_str(fastastr:str,astart:int,aend:int,strand:int=0):
    '''
    return part of reference genome based on 1-based indices
    '''
    str_ = fastastr[astart-1:aend]
    if strand==16:
        return Seq(str_).reverse_complement()

    return str_

def load_fasta_file(filename:str=""):
    fasta_str = ''
    with open(filename, "r") as myfile:
        fasta_ref = myfile.readlines()    
    for _l in fasta_ref:
            if not _l.startswith('>'):
                fasta_str = fasta_str + _l.strip()            

    return fasta_str


def write_header(hdr,thefile):
    for k in hdr.keys():
        _kv = hdr[k][0]
        thefile.write('@{0}\t'.format(k))
        for _k in _kv.keys():            
            _v = _kv[_k]
            thefile.write('{0}:{1}\t'.format(_k,_v))
        thefile.write('\n')


def filter_read_file(file_name):

    ext_ = extFile(args.input_file)
    
    if False: #can_read_file(args.input_file):
        print('map file to memory ...')
        ff = open(':samfile:',"w")
        file_name = ':samfile:'
    else:
        print('write filtered reads to temporary file ...')
        file_name = next(tempfile._get_candidate_names())
        ff = open(file_name,"w")

    # write reads into memory mapped sam file and filters on the fly... crucial for memory issues        
    if ext_ == 'sam':
        a = open(args.input_file,'r')
        for line in a:            
            if line.startswith('@'):
                ff.write(line)
            else: # this is an important line .. apparently low quality reads have additional fields that make the import via csv in pandas very complicated
                l = line.rstrip().split('\t')
                qual = int(l[4])
                if qual>=args.min_qual:            
                    ff.write(line)
        a.close()
        ff.close()            
        
    print('determining sam data columns ... ')  
    headers = []    
    f = open(file_name,'r')    
    # determine number of header lines
    for l in f:        
        if l.startswith('@'):
            headers.append(l)
        else: # assume this is no blank line
            break
    f.close()        
    if l is None: 
        print('No lines found in the file')
        return

    nrfields = l.split()    
    cols = ['C{0}'.format(i + 1) for i in range(len(nrfields))]


    return file_name,headers,cols


def extFile(filename:str):
    return filename.split(".")[-1].lower()
    

def revCompl(seq):
    _bseq = Seq(seq)
    _rev_comp = _bseq.reverse_complement()    
    return str(_rev_comp)


def getAsiteData53(readStart:int,strand:int,offset=args.up_offset):    
    if strand==16: # antisense read
        return readStart-offset
    return readStart+offset


def getAsiteData35(readStop:int,strand:int,offset=args.down_offset):    
    if strand==16: # antisense read
        return readStop+offset
    return readStop-offset


def augment_sam_data(dfin:pd.DataFrame):

    dfin['len'] = dfin['read_str'].str.len()
    dfin['right_pos']=dfin.left_pos+dfin.len-1

    dfin['begin_read']=(dfin.dir==16)*dfin.right_pos + (dfin.dir==0)*dfin.left_pos
    dfin['end_read']=(dfin.dir==16)*dfin.left_pos + (dfin.dir==0)*dfin.right_pos

    dfin['3->5'] = dfin.apply(lambda x:getAsiteData35(x.end_read,x.dir),axis=1)
    dfin['5->3'] = dfin.apply(lambda x:getAsiteData53(x.begin_read,x.dir),axis=1)

    return dfin

fasta_data = load_fasta_file(args.fasta_file)

def add_fasta_sequence(offset:int,strand:int=0,ntup:int=args.nt_upstream,ntdown:int=args.nt_downstream):
    # offset is zero based here     
    # assume position = 80 (=offset) -> (80-nt_upstream) : (80 + (ndown+2))
    # e.g. (80-50) : (80+2+48) = 30:130 => with 0 based = (30,31,32 ... 129) -> element[50] == offset
    # antisense
    # (offset-2-48) : (offset+50) -> 30:130 => range(30,130) -> element[50] == offset
    
    minr_ = offset-ntup
    maxr_ = offset + 2 + ntdown

    if strand==16:
        minr_ = offset - 2 - ntdown
        maxr_ = offset + ntup
    
    return get_fasta_str(fastastr=fasta_data,astart=minr_,aend=maxr_,strand=strand)

def check_access(file, mode):
   ''' Check if file is already open so it cannot be written to. If open, return False else True.
   '''
   try:
      open(file, mode)
      return True
   except PermissionError:
      return False


def add_prefix(filename):
    ''' Add time stamp as a prefix to a ame of a file that needs to be written to.
    '''
    prefix = md5(str(localtime()).encode('utf-8')).hexdigest()
    return f"{prefix}_{filename}"


if args.gff_file!="":
    gns_def = create_gene_definition_file(args.gff_file)

def augment_with_gene_info(dfin:pd.DataFrame):
    
    str_map = {0:'+',16:'-'}
    
    def find_gene(pos:int,strand:int):
        
        gene_ = (str_map[strand]==gns_def.strand) & (pos >= gns_def.start) & (pos<= gns_def.stop)    
        if sum(gene_)==1:
            entry_ = gns_def[gene_]
            gene_name = entry_.gene.values[0]            
            locus_tag = entry_.locus_tag.values[0]
            
            from_beg = pos - entry_.start.values[0]
            if strand==16:
                from_beg = entry_.stop.values[0]-pos

            gene_length = entry_.len
            in_orf = (pos>=entry_.gene_orf_start_tenperc) & (pos<=entry_.gene_orf_stop_tenperc)
            return gene_name, locus_tag, from_beg, in_orf.values[0], gene_length.values[0]

        return "","",None,False,None
        
    dfin['gene'],dfin['locus_tag'],dfin['offset'],dfin['in_orf_90'], dfin['gene_length']= zip(*map(find_gene,dfin.position,dfin.strand))
    dfin = dfin.replace(np.NaN,"")
    dfin = dfin[['genome','position','strand','gene','locus_tag','gene_length','offset','in_orf_90','count','sequence']]
    return dfin

def process_block(dfIn:pd.DataFrame, blocknr:int):

    global df_35, df_53
    
    # copy organism from first row
    org_ = dfIn.iloc[0][2] #hdrs[1].split("\t")[1].split(":")[1]

    dfIn.rename(columns={'C5': 'qual', 'C10': 'read_str', 'C2': 'dir', 'C4': 'left_pos', 'C11': 'phred_scores'}, inplace=True)
    
    def concat_frames(dfTot,dfIn,direction='3->5'):
        if dfTot.shape[0]>0:
            dfIn.rename(columns={'Count':'Count2'},inplace=True)        
            _out=pd.merge(dfTot,dfIn,left_on=['org','C2',direction,'Strand'],right_on=['org','C2',direction,'Strand'],how='outer')
            _out.replace(np.nan,0,inplace=True)
            _out.Count = _out.Count + _out.Count2
            _out.drop(columns='Count2',inplace=True)
            _out.Count = _out.Count.astype("int")
            return _out
        else:
            return dfIn


    print("processing block {0}...".format(blocknr))
    dfaug = augment_sam_data(dfIn)    
    
    print("summarizing 3->5 data  ...")
    # also aggregate on direction to include the strand info
    _df_35 = dfaug.groupby('3->5').agg(
            Count=pd.NamedAgg(column='C1',aggfunc='count'),
            Strand=pd.NamedAgg(column='dir',aggfunc='max')).reset_index()

    _df_35['org'] = org_
    _df_35['C2'] = _df_35['3->5']-1
    _df_35 = _df_35[['org','C2','3->5','Count','Strand']]
    
    _df_35 = _df_35[_df_35.Count>=args.min_hits]

    df_35 = concat_frames(df_35,_df_35)
    
    # repeat for 5->3 direction
    # also aggregate on direction to include the strand info
    print("summarizing 5->3 data  ...")

    _df_53 = dfaug.groupby('5->3').agg(
        Count=pd.NamedAgg(column='C1',aggfunc='count'),
        Strand=pd.NamedAgg(column='dir',aggfunc='max')).reset_index()
    _df_53['org'] = org_
    _df_53['C2'] = _df_53['5->3']-1
    _df_53 = _df_53[['org','C2','5->3','Count','Strand']]
    _df_53 = _df_53[_df_53.Count>=args.min_hits]

    df_53 = concat_frames(df_53,_df_53,direction='5->3')


def process_blocks(filename,hdrs,cols, cs=args.nlines):
    block_nr = 1
    with pd.read_csv(filename, sep='\t', skiprows=len(hdrs), header=None, names=cols,chunksize=cs) as reader:        
        for chunk in reader:
            process_block(chunk.copy(),block_nr)
            block_nr += 1

#%% the main routine
global df_35, df_53

df_35 = df_53 = pd.DataFrame()

def main():
    global df_35, df_53
    
    outputfiles = {}

    if check_arguments():        
        print("Analyzing {0}".format(args.input_file))
        print("\tOffset 3->5: +/-{0}".format(args.down_offset))
        print("\tOffset 5->3: +/-{0}".format(args.up_offset))
        print("\tCut-off reporting: minimal {0} reads".format(args.min_hits))
        print("\tReporting from -{0} to {1} nt\n\n".format(args.nt_upstream,args.nt_downstream))
    else:
        exit(1)
    
    _filename,_hdrs,_cols = filter_read_file(args.input_file)        
    

    # generate the summaries df_35 and df_53
    process_blocks(_filename,_hdrs,_cols)

    # sort values according to position
    df_35 = df_35.sort_values("C2")
    df_53 = df_53.sort_values("C2")
    
    # dir write to bed file(s)
    ext_ = args.input_file.split('.')[-1]
    output_prefix = args.input_file.strip("."+ext_)

    file_ = output_prefix+"_35.bed"    
    if not check_access(file_,"w"):
        file_ = add_prefix("35.bed")
    df_35.to_csv(file_,header=None,index=False,sep='\t')
    outputfiles['35bed']=file_

    file_ = output_prefix+"_53.bed"
    if not check_access(file_,"w"):
        file_ = add_prefix("53.bed")
    df_53.to_csv(file_,header=None,index=False,sep='\t')
    outputfiles['53bed']=file_


    # add fasta parts
    print("add fasta sequences to results")
    if df_35.shape[0]>0:
        df_35['seq']=df_35.apply(lambda x:add_fasta_sequence(x['3->5'],x.Strand),axis=1)
        df_35.drop(columns=['C2'],inplace=True)        
        df_35.columns = ['genome','position','count','strand','sequence']                

    if df_53.shape[0]>0:
        df_53['seq']=df_53.apply(lambda x:add_fasta_sequence(x['5->3'],x.Strand),axis=1)
        df_53.drop(columns=['C2'],inplace=True)        
        df_53.columns = ['genome','position','count','strand','sequence']        

    # augment with gene info if specified
    if (args.gff_file!=""):
        print("augmenting with gene definition data")
        df_35 = augment_with_gene_info(df_35)
        df_53 = augment_with_gene_info(df_53)

    # saving files to csv, and create temporary file names if they're opened already

    file_ = output_prefix+"_35_full.csv"
    if not check_access(file_,"w"):
        file_ = add_prefix("35_full.csv")            
    df_35.to_csv(file_,index=False)
    outputfiles['35csv']=file_

    file_ = output_prefix+"_53_full.csv"
    if not check_access(file_,"w"):
        file_ = add_prefix("35_full.csv")                    
    df_53.to_csv(file_,index=False)
    outputfiles['53csv']=file_

    print("finished\n\n")
    print("output written to:")
    for k_ in outputfiles.keys():
        print("\t"+outputfiles[k_])    

    print('\nremove temporary file ...')
    os.remove(_filename)
# %%
#run main program .. 
main()
