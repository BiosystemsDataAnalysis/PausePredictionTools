#%% import/install some modules

## from https://stackoverflow.com/questions/44210656/how-to-check-if-a-module-is-installed-in-python-and-if-not-install-it-within-t

import sys
import subprocess
import pkg_resources
import tempfile   
from os.path import exists
from map_functions import *
import multiprocessing
from collections import Counter

required = {"ipython", "pandas","numpy","biopython","plotly","typed-argument-parser"}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    python = sys.executable
    subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)

from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

import plotly.express as px


#%% define the functions
def load_fasta_file(fasta_file_loc:str=""):  
    ''' loads the fasta file with sequence data .. skips the lines starting with ">"
    '''
   
    print("loading the fasta file {0}".format(fasta_file_loc))
    fasta_file = fasta_file_loc
    fasta_str = ''
    with open(fasta_file, "r") as myfile:
        fasta_ref = myfile.readlines()    

    for _l in fasta_ref:
            if not _l.startswith('>'):
                fasta_str = fasta_str + _l.strip()
            
    return fasta_str

def load_gene_defs(fasta_str:str="", gff_file:str=""):
    ''' loads the gene definitions from the gff file '''
    
    def get_fasta_str(astart,aend):
        '''
        return part of reference genome based on 1-based indices
        '''
        return fasta_str[astart-1:aend]
        
    # augment the gene definitions and extend the ORF range with site +/- 30 nucleotides

    print("loading the gene definitions file {0}".format(gff_file))
    gns = augment_gff_file(gff_file,rng=30,funGetFastaStr=get_fasta_str)
    
    return gns


def load_sam_file(sam_file_name:str="", saveFormatted:bool=False)->pd.DataFrame:
    '''
    Load the sam file to pandas dataframe 
    saveFormatted:bool - save the intermediate data to a pkl file for faster processing at a later stage, default False
    '''

    if not(exists(sam_file_name)):
        print("{0} not found, please enter valid file name/path".format(sam_file_name))

    sam_file = sam_file_name
    pkl_file = sam_file.replace(".sam",".pkl")
    load_sam_data_from_pkl_file = exists(pkl_file)

    val = ["is not","is"]
    print("The **sam** data **{0}** loaded from previous results".format(val[load_sam_data_from_pkl_file]))

    # load data from sam file ..
    if not(load_sam_data_from_pkl_file):        

        print('determining headers sam file ... ')
        headers = []
        f = open(sam_file,'r')
        # determine number of header lines
        for l in f:
            if l.startswith('@'):
                headers.append(l)
            else:
                break
        f.close()
        nr_headers = len(headers)

        cols = ['C{0}'.format(i + 1) for i in range(20)]    
        
        df_sam = pd.read_csv(sam_file, sep='\t', skiprows=nr_headers, header=None, names=cols, engine='python')
        df_sam.rename(columns={'C5': 'qual', 'C10': 'sam_str', 'C2': 'dir', 'C4': 'left_pos', 'C11': 'phred_scores'},inplace=True)
        df_sam = df_sam.drop(['C6','C7','C8','C9','C12','C13','C14','C15','C16','C17','C18','C19','C20'],axis=1)
        
        df_sam.left_pos = df_sam.left_pos.astype('int')
        df_sam.qual = df_sam.qual.astype('int')

        if saveFormatted:
            df_sam.to_pickle(pkl_file)
    else:
        df_sam = pd.read_pickle(pkl_file)

    return df_sam



def mythreadfunc(index,df_work,multi_dict_,d_normal,d_reverse,blocks,fasta_str,gns,offset35:int=12,offset53:int=14):
        rng = range(blocks[index][0],blocks[index][1])
        reg_out,nrm_out,rev_out = process_sam_data(df_work.iloc[rng],fasta_str,gns, index, offset35,offset53)
        multi_dict_[index] = reg_out
        d_normal[index]=nrm_out
        d_reverse[index]=rev_out

#%% define the plot functions

def detpos(vecdata, offset35:int=12,offset53:int=14):
      
    #/offset53=_offset53
    #offset35=_offset35

    # print(offset53,offset35)

    strand = vecdata[0]
    start_read = vecdata[1]
    read_len = vecdata[2]
    position_I = vecdata[3]
    
    stop_read = start_read + (read_len-1)
    if strand==16:
        stop_read = start_read - (read_len - 1)  

    result_53 = [getAsiteData53(start_read,p,strand,offset53) for p in position_I]
    result_35 = [getAsiteData35(stop_read,p,strand,offset35) for p in position_I]
  
    return [result_53,result_35]


# add ORF 10% information
def addORFinfo(df):
    
    orf_range = 0.05 # percentage 
    df['gene_len']=np.abs(df['gene_stop']-df['gene_start'])
    df['gene_orf_start_tenperc']=(df.dir==0)*(df.gene_start+np.round(df.gene_len*orf_range,0))+(df.dir==16)*(df.gene_start-np.round(df.gene_len*orf_range,0))
    df['gene_orf_stop_tenperc']=(df.dir==0)*(df.gene_stop-np.round(df.gene_len*orf_range,0))+(df.dir==16)*(df.gene_stop+np.round(df.gene_len*orf_range,0))

    df['within_90_ORF']=((df.dir==0) & ((df.begin_read >= df.gene_orf_start_tenperc) & (df.end_read<= df.gene_orf_stop_tenperc))) | \
        ((df.dir==16) & ((df.begin_read <= df.gene_orf_start_tenperc) & (df.end_read>= df.gene_orf_stop_tenperc)))

    return df


def createPlotData(df_data:pd.DataFrame, postfix="", cutoff_dist=20):
    # i53vec_15_11 = list(chain(*df_data.minI53))
    # i35vec_15_11 = list(chain(*df_data.minI35))

    aa_found = df_data.columns[df_data.columns.str.contains("min")]
    aa_found = list(aa_found.str.strip("min53").unique())
        
    df_wrk_plt = None
    
    if len(aa_found)==0:
        print("no distances to know AAs were found")
        return None

    for AA in aa_found:
        col53_,col35_ = "min" + AA + "53", "min" + AA + "35"
        aa53vec_ = df_data[col53_]
        aa35vec_ = df_data[col35_]

        if postfix=="":
            colName = "5->3"
        else:
            colName = "5->3 ({0})".format(postfix)

        aa53_ = Counter(aa53vec_)
        df_aa53 = pd.DataFrame.from_dict(aa53_,orient='index')
        df_aa53.columns = [colName]

        aa35_ = Counter(aa35vec_)
        df_aa35 = pd.DataFrame.from_dict(aa35_,orient='index')
        if postfix=="":
            colName = "3->5"
        else:
            colName = "3->5 ({0})".format(postfix)        

        df_aa35.columns = [colName]

        df_wrk_aa53 =df_aa53[(df_aa53.index>-cutoff_dist) & (df_aa53.index<cutoff_dist)]
        df_wrk_aa35 =df_aa35[(df_aa35.index>-cutoff_dist) & (df_aa35.index<cutoff_dist)]

        df_wrk = df_wrk_aa35.merge(df_wrk_aa53,left_index=True,right_index=True,how='outer')

        df_wrk_plt_= df_wrk.stack().reset_index()    
        df_wrk_plt_.columns=['position','serie','counts']
        df_wrk_plt_['amino acid']=AA

        if df_wrk_plt is None:
            df_wrk_plt = df_wrk_plt_.copy()
        else:
            df_wrk_plt = pd.concat([df_wrk_plt,df_wrk_plt_],axis=0)

    combined_ = df_wrk_plt.groupby(['position','serie']).sum().reset_index()
    combined_['amino acid']='combined'
    df_wrk_plt = pd.concat([df_wrk_plt,combined_],axis=0)
    return df_wrk_plt[['position','amino acid','serie','counts']]


def createReadDistData(df_data):
    df_ret = df_data.groupby('read_len').count()['C1']    
    return pd.DataFrame(df_ret)
    # return pd.DataFrame(df_ret,columns=['frequency'])



def createAsiteDistributionPlotData(df_in:pd.DataFrame):

    # create filters
    genes_summary_ = df_in.groupby('gene').agg(['count','min'])['gene_len']      
    
    # calculate relative frequency
    genes_rel_=genes_summary_['count']/genes_summary_['min']
    genes_rel_.name='relfreq'
    
    top5_ = genes_rel_.sort_values(ascending=False).iloc[0:5]
    top10_ = genes_rel_.sort_values(ascending=False).iloc[0:10]
    morethan10 = genes_summary_[genes_summary_['count']>10]

    all_reads = createPlotData(df_in)
    all_reads['condition'] = 'all'

    df_reads_ = all_reads.copy()

    reads_ORF_ = createPlotData(df_in[df_in.within_90_ORF])
    reads_ORF_['condition']= 'within 90% ORF'    
    df_reads_ = pd.concat([df_reads_,reads_ORF_],axis=0)

    filter = (df_in.gene.isin(morethan10.index.values)) & (df_in.within_90_ORF)
    read_ORF_and_10 = createPlotData(df_in[filter])
    read_ORF_and_10['condition'] = 'ORF + >10'
    df_reads_ = pd.concat([df_reads_,read_ORF_and_10],axis=0)

    filter = ((df_in.gene.isin(morethan10.index.values))) & (df_in.within_90_ORF) & (~(df_in.gene.isin(top5_.index.values)))
    if (filter.sum()>0):
        read_ORF_and_10_minTop5 = createPlotData(df_in[filter])
        read_ORF_and_10_minTop5['condition'] = 'ORF + >10 - top5'
        df_reads_ = pd.concat([df_reads_,read_ORF_and_10_minTop5],axis=0)
    
    
    filter = ((df_in.gene.isin(morethan10.index.values))) & (df_in.within_90_ORF) & (~(df_in.gene.isin(top10_.index.values)))
    if (filter.sum()>0):
        read_ORF_and_10_minTop10 = createPlotData(df_in[filter])
        read_ORF_and_10_minTop10['condition'] = 'ORF + >10 - top10'
        df_reads_ = pd.concat([df_reads_,read_ORF_and_10_minTop10],axis=0)
    
    return df_reads_


def add_aa_scores(dfin:pd.DataFrame,aa:str='I',offset35:int=12,offset53:int=14):
    
    # copy I data vectors and determine distances to reads and I position, in detpos the offsets are determined 
    # they can be set by the global variables _offset53 and _offset35
    cc = dfin[['dir','begin_read','read_len',aa]].apply(lambda x:detpos(x),axis=1)
    ccc=cc.apply(pd.Series)
    
    col35_, col53_= aa+"35_"+str(offset35), aa+"35_"+str(offset53)
    ccc.columns = [col53_,col35_]

    # find minimum distance between I and A-site ?
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

#%% the main plot routine
def make_plots(dfw:pd.DataFrame, pkl_file_name:str, sample_type:str="",offset35:int=12,offset53:int=14):
    kys_,cdns_ = get_genetic_code()
    
    output_file_name = pkl_file_name.replace(".pkl","_with_AAs.pkl")
    
    if exists(output_file_name):
        print("reading input data from previously generated data {0}".format(output_file_name))
        dfm = pd.read_pickle(output_file_name)
    else:
        
        # focus on mapped to known genes only
        dfmapped_on_genes = dfw[~(dfw.gene=="")].copy()
        dfm = addORFinfo(dfmapped_on_genes)
        
        print("adding distances for the different aminoacids ")
        for k_ in kys_:            
            dfm = add_aa_scores(dfm,k_)

        print("storing data in {0}".format(output_file_name))
        dfm.to_pickle(output_file_name)


    file_name_pickle = output_file_name.replace(".pkl","plot_data_{0}_{1}.pkl".format(offset53,offset35))

    if exists(file_name_pickle):
        print("reading data from previously generated file {0}".format(file_name_pickle))
        df_n_m = pd.read_pickle(file_name_pickle)
    else:
        df_n_m = createAsiteDistributionPlotData(dfm)
        df_n_m.to_pickle(file_name_pickle)

    work_data = df_n_m
    plot_data = work_data[(work_data.position>=-offset35) & (work_data.position<=offset53)]

    # plot 1
    _title = r"A-site position relative to first codon of different amino acids " + \
    "in {0} sample determined with +{1}/-{2} offsets".format(sample_type,offset53,offset35)

    fig = px.bar(data_frame=plot_data,x='position',y='counts',color='amino acid',facet_col="serie",animation_frame='condition', barmode='group')
    fig["layout"].pop("updatemenus") # optional, drop animation buttons
    fig.update_traces(dict(marker_line_width=0))
    xaxis_dict = dict(tickmode = 'linear',tick0 = 0,dtick = 3)
    fig.update_layout(title_text=_title,xaxis=xaxis_dict,xaxis2=xaxis_dict)

    fig.show();
    file_name_html = output_file_name.strip(".pkl")+"_all_{0}_{1}_wo_sel.html".format(offset53,offset35)
    fig.write_html(file_name_html)



#%% the main routine for processing
from tap import Tap # typed-argument-parser
import sys
interactive_mode = hasattr(sys, 'ps1')

#gff_file = "./reference_files/wt-prspb-amyI-GFF3.gff"
#fasta_file = "./reference_files/WT-Prspb-amyI-FASTA.fa"

if interactive_mode:
 sys.argv = ['script.py', '--sam_file', './input/1X_PEG-Spin_filtered_SAM.sam', 'output_file','./output/1X_PEG-Spin_filtered.pkl']
#  sys.argv = ['script.py', '--input_file', 'test_out.sam','--fasta_file','WT-Prspb-amyI-FASTA.fa']

class myargs(Tap):
    
    sam_file: str  # the input sam(or bam) file
    gff_file: str # the file with gene definitions
    fasta_file: str # the file with nucleotide sequences
    #output_file:str # the output pkl file
    nr_cores: int=1 # the number of cores to use
    mq:int=41 # the minimum mapping quality of a read
    si:bool=False # save the intermediate results to (pkl) file
    o53:int=14 # the offset from 5->3 direction
    o35:int=12 # the offset from 3->5 direction
    log:bool=False # create a log file
    


# add this line for processing in multiple threads
if __name__ ==  '__main__': 
    
    args = myargs().parse_args()        
    sam_file_name = args.sam_file
    output_file_name = args.sam_file.lower().replace(".sam","_ASITE.pkl")

    print('processing file {0}'.format(sam_file_name))

    if args.log:
        (fd, filename) = tempfile.mkstemp()
        temp_stdout = sys.stdout
        print("logging to file {0}".format(filename))
        sys.stdout = open("./"+filename,'wt')


    if not(exists(output_file_name)):

        fasta_str = load_fasta_file(args.fasta_file)
        gns = load_gene_defs(fasta_str, args.gff_file)
        df_sam = load_sam_file(sam_file_name, args.si)
        #     #  filter the reads with low mapping quality
        df_work_data = df_sam[df_sam.qual >= args.mq].copy()

        maxcpu = args.nr_cores             
        maxcpu = min(multiprocessing.cpu_count(),maxcpu)
        print('running on {0} threads'.format(maxcpu))

        blocksize = df_work_data.shape[0]//maxcpu
        blocks = {}
        for b in range(maxcpu):
            blocks[b]=[(b)*blocksize,(b+1)*blocksize]

        blocks[b]=[blocks[b][0],df_work_data.shape[0]]

        print("nr cpus = {0}".format(maxcpu))
            
        multiprocessing_dict = multiprocessing.Manager().dict()
        dct_normal = multiprocessing.Manager().dict()
        dct_reverse = multiprocessing.Manager().dict()
        multiprocessing_loop_  = []

        for c in range(maxcpu):    
            _process = multiprocessing.Process(target=mythreadfunc, args=(c,df_work_data, multiprocessing_dict,dct_normal,dct_reverse,blocks,fasta_str,gns,args.o35,args.o53))
            multiprocessing_loop_.append(_process)
        
        for process_  in multiprocessing_loop_ :
            process_.start()

        for process_  in multiprocessing_loop_ :
            process_.join()

        df = multiprocessing_dict[0].copy()
        normal_reads = dct_normal[0].copy()
        reverse_read = dct_reverse[0].copy()
        
        for i in range(1,maxcpu):
            if not(multiprocessing_dict[i] is None):
                df = pd.concat([df,multiprocessing_dict[i]],axis=0)
                normal_reads += dct_normal[i]
                reverse_read += dct_reverse[i]
        
        df.to_pickle(output_file_name)
    
    else:        
        print("reading from previous results {0}".format(output_file_name))
        df = pd.read_pickle(output_file_name)


    make_plots(df, output_file_name, offset35=args.o35, offset53=args.o53)

    if args.log:
        sys.stdout = temp_stdout
    
    print("finished..")

