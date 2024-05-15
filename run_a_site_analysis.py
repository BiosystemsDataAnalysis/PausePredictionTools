#%% import/install some modules

## from https://stackoverflow.com/questions/44210656/how-to-check-if-a-module-is-installed-in-python-and-if-not-install-it-within-t

import sys
import subprocess
from importlib import metadata
#import pkg_resources
import tempfile   
from os.path import exists, dirname, realpath, basename, join
from map_functions import *
import multiprocessing
from collections import Counter
import psutil

if sys.platform == "linux" or sys.platform == "linux2":    
    import resource

import h5py

import datatable as dt
import threading

import time
import os
required = {"ipython", "resource","pandas","numpy","biopython","plotly","datatable","typed-argument-parser","psutil"}

def checkrequired():
    missing = []
    for r in required:
        try:
            metadata.version(r)
        except:
            missing.append(r)
    
    return missing

missing = checkrequired()

lock = threading.Lock()

if missing:
    python = sys.executable
    subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)

from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

import plotly.express as px
import plotly.graph_objects as go

from tap import Tap # typed-argument-parser
import sys
interactive_mode = hasattr(sys, 'ps1')
from enum import IntEnum

import pandas as pd 
import warnings
warnings.filterwarnings('ignore',category=pd.io.pytables.PerformanceWarning)


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

def count_lines_in_file(file_name:str="")->int:
    ''' count the number of lines in a text file
    '''
    count = 0
    with open(file_name,'r') as f:
      for _ in f:
        count += 1
    return count


def len_headers_sam_file(sam_file):    
    f = open(sam_file,'r')
    print('determining headers sam file ... ')
    headers = []
    # determine number of header lines
    for l in f:
        if l.startswith('@'):
            headers.append(l)
        else:
            break
    f.close()
    return len(headers)

    
#%% define the plot functions

def createPlotData_HDF(pDict:dict, postfix="", cutoff_dist=20):
    
    # i53vec_15_11 = list(chain(*df_data.minI53))
    # i35vec_15_11 = list(chain(*df_data.minI35))

    print("create plotting data")
    
    with pd.HDFStore(pDict[PROJECT_INFO.STORE]) as hdf_store:
                
        init_df =  hdf_store.get("/augmented/block_0")
        aa_found = init_df.columns[init_df.columns.str.contains("min")]
        # remove characters m,i,n,5 and 3
        aa_found_id = list(aa_found.str.strip("min53").unique())
    
        if len(aa_found)==0:
            print("no distances to know AAs were found")
            return None

        #if "/raw_aa_plot_data" in hdf_store.keys():
        try:
            print("loading raw plot data from file")
            df_data = hdf_store.get("/raw_aa_plot_data")            
        except:
            # load all data in array .. not actually a block function.. for later
            print("processing combined plot data")
            df_data = pd.DataFrame()
            _keys = [k for k in hdf_store.keys() if k.startswith("/augmented/block_")]
            for k in _keys:
                df = hdf_store.get(k)                
                df = df[['gene','within_90_ORF'] + aa_found.values.tolist()]
                # remove entries that have no gene mapping
                df = df[df.gene!=""]
                if(df_data.shape[0]==0):
                    df_data = df.copy()
                else:            
                    if not(df.empty):                           
                        df_data = pd.concat([df_data,df],axis=0).copy()
           
            hdf_store.put("/raw_aa_plot_data",df_data)

        
        df_data_work = df_data.copy()

        for filter in range(5):
            # read filter data from store            

            if filter==1:
                filter_ = df_data.within_90_ORF
            if filter==2:
                morethan10 = hdf_store.get("/filter_morethan10")
                filter_ = (df_data.gene.isin(morethan10.index.values)) & df_data.within_90_ORF            
            if filter==3:    
                morethan10 = hdf_store.get("/filter_morethan10")
                top5_ = hdf_store.get("/filter_top5")
                filter_ = ((df_data.gene.isin(morethan10.index.values))) & (df_data.within_90_ORF) & (~(df_data.gene.isin(top5_.index.values)))
            if filter==4:
                morethan10 = hdf_store.get("/filter_morethan10")
                top10_ = hdf_store.get("/filter_top10")
                filter_ = ((df_data.gene.isin(morethan10.index.values))) & (df_data.within_90_ORF) & (~(df_data.gene.isin(top10_.index.values)))
        
            if (filter>0):
                df_data_work = df_data[filter_].copy()
                            
            df_wrk_plt_export = None
            for AA in aa_found_id:
                col53_,col35_ = "min" + AA + "53", "min" + AA + "35"
                
                # aa53vec_ = df_data_work[col53_]
                # aa35vec_ = df_data_work[col35_]                

                if postfix=="":
                    colName53 = "5->3"
                    colName35 = "3->5"
                else:
                    colName53 = "5->3 ({0})".format(postfix)
                    colName35 = "3->5 ({0})".format(postfix)
                    
                # make a frequency tables
                df_aa53_new  = df_data_work[col53_].value_counts().to_frame()
                df_aa53_new.columns = [colName53]
                                
                df_aa35_new  =  df_data_work[col35_].value_counts().to_frame()                
                df_aa35_new.columns = [colName35]
                               
                # initialize full vector
                df_wrk_aa35_new = pd.DataFrame(index=range(-cutoff_dist,cutoff_dist+1),data=np.zeros((cutoff_dist*2+1,1)),columns=['3->5'])
                # filter for cutoff distance
                df_wrk_aa35_ = df_aa35_new.loc[ (df_aa35_new.index>-cutoff_dist) & (df_aa35_new.index<cutoff_dist)]
                # combine both
                df_wrk35_new = df_wrk_aa35_new.add(df_wrk_aa35_,fill_value=0)
                
                # initialize full vector
                df_wrk_aa53_new = pd.DataFrame(index=range(-cutoff_dist,cutoff_dist+1),data=np.zeros((cutoff_dist*2+1,1)),columns=['5->3'])
                # filter for cutoff distance
                df_wrk_aa53_ = df_aa53_new.loc[ (df_aa53_new.index>-cutoff_dist) & (df_aa53_new.index<cutoff_dist)]
                # combine both
                df_wrk53_new = df_wrk_aa53_new.add(df_wrk_aa53_,fill_value=0)
                      
                # combine 53 and 35 series for the AA
                df_wrk_new = df_wrk35_new.merge(df_wrk53_new,left_index=True,right_index=True,how='outer')
                # stack results to appropriate for plotting
                df_wrt_plt_new = df_wrk_new.stack().reset_index()
                # rename columns for plotting options
                df_wrt_plt_new.columns=['position','serie','counts']
                # assign for the current Amino Acid
                df_wrt_plt_new['amino acid']=AA
                
                # concatenate plot results                    
                if df_wrk_plt_export is None:
                    df_wrk_plt_export = df_wrt_plt_new.copy()
                else:
                    df_wrk_plt_export = pd.concat([df_wrk_plt_export,df_wrt_plt_new],axis=0)                                     
            
            combined_new = df_wrk_plt_export.groupby(['position','serie']).sum(numeric_only=True).reset_index()
            combined_new['amino acid']='combined'
            df_wrk_plt_export = pd.concat([df_wrk_plt_export,combined_new],axis=0)
            
            # store tables
            hdf_store.put("plotdata/filter_{0}".format(filter),df_wrk_plt_export[['position','amino acid','serie','counts']])


summary_dict = multiprocessing.Manager().dict()
def summary_block_thread(index, store_ptr,sema,verbose):            
    if verbose:
        print("processing block {0} for summary".format(index))
    df = store_ptr.get(index)      
    genes_summary__ = df.groupby('gene').agg(['count','min'])['gene_len']            
    sema.release()    
    with lock:
        summary_dict[index]=genes_summary__
        

def createAsiteDistributionPlotData_HDF(pDict:dict):

    with pd.HDFStore(pDict[PROJECT_INFO.STORE]) as hdf_store:
    # load raw filter data only if necessary
        if ("/genes_summary_" in hdf_store.keys()) and  ("/filter_morethan10" in hdf_store.keys()):            
            print("load distribution summary data")
            genes_summary_ = hdf_store.get("/genes_summary_")

        else: # create filters and summary
            print("process distrubtion data per block")
            
            multiprocessing_loop_  = []
            sema = multiprocessing.Semaphore(pDict[PROJECT_INFO.CPU]);
            _keys = [k for k in hdf_store.keys() if k.startswith("/augmented/block_")]
            for k in _keys:                
                # countdown available semaphores
                sema.acquire();
                _process = multiprocessing.Process(target=summary_block_thread, args=(k,hdf_store,sema, pDict[PROJECT_INFO.ARGS].verbose))            
                multiprocessing_loop_.append(_process)
                _process.start()
            
            # wait for all processes to finish    
            for process_  in multiprocessing_loop_ :
                process_.join()                                
                        
            print("joining summary blocks")
            genes_summary_= summary_dict[_keys[0]]
           
            for k in _keys[1:]:
                _df = summary_dict[k]
                _df = genes_summary_.merge(_df,left_on=['gene'],right_on=['gene'],how='outer')
                _df['count'] = _df[['count_x','count_y']].sum(axis=1).astype('int')
                _df['min'] = _df[['min_x','min_y']].min(axis=1)            
                _df = _df.drop(['count_x','count_y','min_x','min_y'],axis=1).copy()
                genes_summary_ = _df.copy()
            
            hdf_store.put("/genes_summary_",genes_summary_)

    # close hdf_store

    with pd.HDFStore(pDict[PROJECT_INFO.STORE]) as hdf_store:
        
        if not("/filter_top5" in hdf_store.keys()): # assume the others are also created                    
            print("calculate plotting filters")
            # calculate relative frequency
            genes_rel_=genes_summary_['count']/genes_summary_['min']
            genes_rel_.name='relfreq'
        
            top5_ = genes_rel_.sort_values(ascending=False).iloc[0:5]
            top10_ = genes_rel_.sort_values(ascending=False).iloc[0:10]
            morethan10 = genes_summary_[genes_summary_['count']>10]        
            
            # store filters
            hdf_store.put("/filter_top5",top5_)
            hdf_store.put("/filter_top10",top10_)
            hdf_store.put("/filter_morethan10",morethan10)
        
    # close hdf_store, recreate plot data from here
    
    createPlotData_HDF(pDict)        

    with pd.HDFStore(pDict[PROJECT_INFO.STORE]) as hdf_store:
        all_reads = hdf_store.get("plotdata/filter_0")    
        all_reads['condition'] = 'all'
        df_reads_ = all_reads.copy()

        reads_ORF_ = hdf_store.get("plotdata/filter_1")    
        reads_ORF_['condition']= 'within 90% ORF'    
        df_reads_ = pd.concat([df_reads_,reads_ORF_],axis=0)

        read_ORF_and_10 =  hdf_store.get("plotdata/filter_2")
        read_ORF_and_10['condition'] = 'ORF + >10'
        df_reads_ = pd.concat([df_reads_,read_ORF_and_10],axis=0)

        read_ORF_and_10_minTop5 =  hdf_store.get("plotdata/filter_3")
        read_ORF_and_10_minTop5['condition'] = 'ORF + >10 - top5'
        df_reads_ = pd.concat([df_reads_,read_ORF_and_10_minTop5],axis=0)

        read_ORF_and_10_minTop10 =  hdf_store.get("plotdata/filter_4")
        read_ORF_and_10_minTop10['condition'] = 'ORF + >10 - top10'
        df_reads_ = pd.concat([df_reads_,read_ORF_and_10_minTop10],axis=0)
    
    return df_reads_
   

class ORF_FILTER(IntEnum):
    NONE = 0
    MTOP10 = 1
    MTOP20 = 2
    MTOP50 = 3
    L10 = 4


def make_ORF_plot_data(dataIn:pd.DataFrame,dir53:bool=True):

    if dir53:
        _lbl = '5->3'
        #posAsite = 'posAsite53'
        dataIn = dataIn[dataIn.direction==53].copy()

    if not dir53:
        _lbl = '3->5'
        dataIn = dataIn[dataIn.direction==35].copy()
        #posAsite = 'posAsite35'

    print("reformating to ORF plot data direction {0}".format(_lbl))

    # dfm_wrk = dataIn.groupby([posAsite,'dir']).count()['C1'].reset_index()
    dfm_wrk = dataIn.groupby(['posAsite','dir']).count()['count'].reset_index()
    
    df_plot = dfm_wrk[(dfm_wrk['posAsite']>-50) & (dfm_wrk['posAsite']<50) ].copy()
    df_plot.columns = ['position','dir','count']
    total_plot = df_plot.groupby('position').sum()['count'].reset_index()
    total_plot['series']=_lbl + ' (total)'

    df_plot['series']=_lbl
    df_plot.loc[df_plot.dir==16,'series']=_lbl+'(*)'
    #df_plot['series'].loc[df_plot.dir==16]=_lbl+'(*)'

    dataIn = dataIn[dataIn.direction=='35'].copy().drop(columns=['dir'],axis=1,inplace=True)
    
    return df_plot,total_plot


orf_dict = multiprocessing.Manager().dict()

def orf_block_thread(index, store_ptr,sema,verbose):        
    
    if verbose:
        print("processing block {0} for ORF".format(index))
    df = store_ptr.get(index)
    dfm_35 = df.groupby(['posAsite35','dir','gene']).agg(['count','min'])['gene_len'].reset_index()
    dfm_53 = df.groupby(['posAsite53','dir','gene']).agg(['count','min'])['gene_len'].reset_index()
    
    sema.release()
    
    with lock:
        orf_dict[index]=[dfm_35,dfm_53]
        


def prepare_ORF_plot_data_HDF(pDict:dict, filter:ORF_FILTER=ORF_FILTER.NONE):
    #
    # output_file_name = join(_op,file_prefix.replace("_ASITE.pkl","_ORF_unfiltered_data_HDF.pkl"))
    
    with pd.HDFStore(pDict[PROJECT_INFO.STORE]) as hdf_store:
    
        if "/orf_plot_data" in hdf_store.keys():    
            print("reading ORF data from previously generated results")
            dfm_tot = hdf_store.get("orf_plot_data")
    
        else:
            print("preparing data for ORF plots... ")

            dfm_35 = dfm_53 = pd.DataFrame()

            _keys = [k for k in hdf_store.keys() if k.startswith("/augmented/block_")]
            _nblocks = pDict[PROJECT_INFO.NBLOCKS]
            
            assert len(_keys)==_nblocks, "file corruption error, not enough data blocks found"
            multiprocessing_loop_  = []
            sema = multiprocessing.Semaphore(pDict[PROJECT_INFO.CPU]);
            for k in _keys:                
                # countdown available semaphores
                sema.acquire();
                _process = multiprocessing.Process(target=orf_block_thread, args=(k,hdf_store,sema, True))            
                multiprocessing_loop_.append(_process)
                _process.start()
            
            # wait for all processes to finish    
            for process_  in multiprocessing_loop_ :
                process_.join()

            print("combining ORF blocks")                             
            dfm_35 = orf_dict[_keys[0]][0]
            dfm_53 = orf_dict[_keys[0]][1]                          
            for k in _keys[1:]:
                _dfm_35 = orf_dict[k][0]
                _dfm_53 = orf_dict[k][1]                  
                _df = dfm_35.merge(_dfm_35,left_on=['posAsite35','dir','gene'],right_on=['posAsite35','dir','gene'],how='outer')
                _df['count'] = _df[['count_x','count_y']].sum(axis=1).astype('int')
                _df['min'] = _df[['min_x','min_y']].min(axis=1)
                dfm_35 = _df.drop(['count_x','count_y','min_x','min_y'],axis=1).copy()
                _df = dfm_53.merge(_dfm_53,left_on=['posAsite53','dir','gene'],right_on=['posAsite53','dir','gene'],how='outer')
                _df['count'] = _df[['count_x','count_y']].sum(axis=1).astype('int')
                _df['min'] = _df[['min_x','min_y']].min(axis=1)
                dfm_53 = _df.drop(['count_x','count_y','min_x','min_y'],axis=1).copy()
                    
                # dfm_35 = orf_dict[_keys[k]][0].drop(['count_x','count_y','min_x','min_y'],axis=1).copy()
                # dfm_53 = orf_dict[_keys[k]][1].drop(['count_x','count_y','min_x','min_y'],axis=1).copy()
            
           
            # for k in _keys:
            #     df = hdf_store.get(k)
            #     _dfm_35 = df.groupby(['posAsite35','dir','gene']).agg(['count','min'])['gene_len'].reset_index()
            #     _dfm_53 = df.groupby(['posAsite53','dir','gene']).agg(['count','min'])['gene_len'].reset_index()
            #     if(dfm_35.shape[0]==0):
            #         dfm_35 = _dfm_35.copy()
            #         dfm_53 = _dfm_53.copy()
            #     else:
            #         _df = dfm_35.merge(_dfm_35,left_on=['posAsite35','dir','gene'],right_on=['posAsite35','dir','gene'],how='outer')
            #         _df['count'] = _df[['count_x','count_y']].sum(axis=1).astype('int')
            #         _df['min'] = _df[['min_x','min_y']].min(axis=1)
            #         dfm_35 = _df.drop(['count_x','count_y','min_x','min_y'],axis=1).copy()
            #         _df = dfm_53.merge(_dfm_53,left_on=['posAsite53','dir','gene'],right_on=['posAsite53','dir','gene'],how='outer')
            #         _df['count'] = _df[['count_x','count_y']].sum(axis=1).astype('int')
            #         _df['min'] = _df[['min_x','min_y']].min(axis=1)
            #         dfm_53 = _df.drop(['count_x','count_y','min_x','min_y'],axis=1).copy()
        
            dfm_35['relfreq']=dfm_35['count']/dfm_35['min']
            dfm_53['relfreq']=dfm_53['count']/dfm_53['min']
        
            dfm_tot = pd.concat([dfm_35,dfm_53],axis=0)
            dfm_tot['direction']=35
            dfm_tot.loc[np.isnan(dfm_tot.posAsite35),'direction']=53
            dfm_tot.direction=dfm_tot.direction.astype('int')

            dfm_tot['posAsite'] = dfm_tot.posAsite53.fillna(0)+dfm_tot.posAsite35.fillna(0)
            dfm_tot = dfm_tot.drop(columns=['posAsite35','posAsite53'],axis=1)

            print("storing intermediate results to file.")
            hdf_store.put("/orf_plot_data",dfm_tot)
        
    
    #create datasets based on filter        

    # all data
    dfm_35 = df_pd_35 = dfm_tot[dfm_tot.direction==35].copy()
    dfm_53 = df_pd_53 = dfm_tot[dfm_tot.direction==53].copy()

    # excluding top 10
    if filter==ORF_FILTER.MTOP10:
        top10_35 = dfm_35.groupby('gene').count().sort_values('relfreq',ascending=False).iloc[0:10].index
        top10_53 = dfm_53.groupby('gene').count().sort_values('relfreq',ascending=False).iloc[0:10].index
        df_pd_35 = dfm_35[~dfm_35.gene.isin(top10_35)]
        df_pd_53 = dfm_53[~dfm_53.gene.isin(top10_53)]

    # excluding top 20
    if filter==ORF_FILTER.MTOP20:
        top20_35 = dfm_35.groupby('gene').count().sort_values('relfreq',ascending=False).iloc[0:20].index
        top20_53 = dfm_53.groupby('gene').count().sort_values('relfreq',ascending=False).iloc[0:20].index
        df_pd_35 = dfm_35[~dfm_35.gene.isin(top20_35)]
        df_pd_53 = dfm_53[~dfm_53.gene.isin(top20_53)]

    # excluding top 50
    if filter==ORF_FILTER.MTOP50:
        top50_35 = dfm_35.groupby('gene').count().sort_values('relfreq',ascending=False).iloc[0:50].index
        top50_53 = dfm_53.groupby('gene').count().sort_values('relfreq',ascending=False).iloc[0:50].index
        df_pd_35 = dfm_35[~dfm_35.gene.isin(top50_35)]
        df_pd_53 = dfm_53[~dfm_53.gene.isin(top50_53)]

    if filter==ORF_FILTER.L10:
        # less than 10 reads
        grouped35 = dfm_35.groupby('gene')
        more_than_10_35 = grouped35.filter(lambda x:x['relfreq'].count()>10).gene.unique()
        grouped53 = dfm_53.groupby('gene')
        more_than_10_53 = grouped53.filter(lambda x:x['relfreq'].count()>10).gene.unique()
        df_pd_35 =  dfm_35[~dfm_35.gene.isin(more_than_10_35)]
        df_pd_53 =  dfm_53[~dfm_53.gene.isin(more_than_10_53)]

    # create plot data, min top10
    df_plot35, df_total35 = make_ORF_plot_data(df_pd_35,False)
    df_plot53, df_total53 = make_ORF_plot_data(df_pd_53)

    # combine the 2 datasets to 1 figure, perhaps skip the detailed bars (df_plot35 & df_plot53)
    df_plot_total= pd.concat([df_total35,df_total53],axis=0)
    df_plot_detail = pd.concat([df_plot35,df_plot53],axis=0)

    return df_plot_total, df_plot_detail


def make_ORF_plot_HDF(pDict:dict, filter:ORF_FILTER=ORF_FILTER.NONE):
    
    df_total, df_detail = prepare_ORF_plot_data_HDF(pDict, filter)

    fig = px.bar(data_frame=df_detail,x='position',y='count',color="series", barmode='overlay')
    _flt = df_total.series == '3->5 (total)'
    fig.add_traces(go.Scatter(x=df_total[_flt].position, y=df_total[_flt]['count'], mode = 'lines',name='3->5 (total)'))

    _flt = df_total.series == '5->3 (total)'
    fig.add_traces(go.Scatter(x=df_total[_flt].position, y=df_total[_flt]['count'], mode = 'lines',name='5->3 (total)'))

    _title = "distribution plot around ORF"
    #_title = 'chloramphenicol (16h) sample min top 10 genes'
    # _title = 'mupirocin sample min top 10 genes'

    fig["layout"].pop("updatemenus") # optional, drop animation buttons
    fig.update_traces(dict(marker_line_width=0))
    # xaxis_dict = dict(tickmode = 'linear',tick0 = 0,dtick = 3)
    fig.update_layout(title_text=_title)    

    fname=join(pDict[PROJECT_INFO.OUTPATH],  pDict[PROJECT_INFO.BASENAME].replace(".sam","_sCDS_HDF.html"))    
    fig.write_html(fname)
    
    print("output file written to {0}".format(fname))

#%% the main plot routine
def make_plots_HDF(pDict:dict):
                
    ofile = os.path.join(pDict[PROJECT_INFO.OUTPATH],pDict[PROJECT_INFO.BASENAME])
    
    if pDict[PROJECT_INFO.ARGS].orf:    
        make_ORF_plot_HDF(pDict,ORF_FILTER.NONE)
    
    o53 = pDict[PROJECT_INFO.ARGS].o53
    o35 = pDict[PROJECT_INFO.ARGS].o35
    file_name_csv = ofile.replace(".sam","_plot_data_{0}_{1}.csv".format(o53,o35))

    df_n_m = createAsiteDistributionPlotData_HDF(pDict)

    work_data = df_n_m
    plot_data = work_data[(work_data.position>=-o35) & (work_data.position<=o53)]
   

    if pDict[PROJECT_INFO.ARGS].csv:
        plot_data.to_csv(file_name_csv)

    _title = pDict[PROJECT_INFO.ARGS].title
    if _title == "" :
        _title = pDict[PROJECT_INFO.BASENAME].replace(".sam","")
                
    _title = r"A-site position relative to first codon of different amino acids " + \
    "in sample ({0}) determined with +{1}/-{2} offsets".format(_title,o53,o35)

    fig = px.bar(data_frame=plot_data,x='position',y='counts',color='amino acid',facet_col="serie",animation_frame='condition', barmode='group')
    fig["layout"].pop("updatemenus") # optional, drop animation buttons
    fig.update_traces(dict(marker_line_width=0))
    xaxis_dict = dict(tickmode = 'linear',tick0 = 0,dtick = 3)
    fig.update_layout(title_text=_title,xaxis=xaxis_dict,xaxis2=xaxis_dict)

    file_name_html = ofile.replace(".sam","_all_{0}_{1}.html".format(o53,o35))
    fig.write_html(file_name_html)
    print("output written to {0}".format(file_name_html))        
               
    separate = False
    # if separate plots are needed
    if separate:
        # from 5->3
        p53 = plot_data.loc[plot_data.serie=="5->3"] 
        p53_all_aa = p53.loc[p53.condition=="all"]
        _title53 = r"A-site position relative to first codon of different amino acids " + \
        "in sample ({0}) determined from 5->3 with +{1} offset".format(_title,o53)
        
        fig = px.bar(data_frame=p53,x='position',y='counts',color='amino acid',animation_frame='condition', barmode='group')
        fig["layout"].pop("updatemenus") # optional, drop animation buttons
        fig.update_traces(dict(marker_line_width=0))
        xaxis_dict = dict(tickmode = 'linear',tick0 = 0,dtick = 3)
        fig.update_layout(title_text=_title53,xaxis=xaxis_dict,xaxis2=xaxis_dict)

        file_name_html = ofile.replace(".sam","_all_53_{0}.html".format(o53,o35))
        fig.write_html(file_name_html)
        
        print("output written to {0}".format(file_name_html))
        
        # from 3->5    
        
        p35 = plot_data.loc[plot_data.serie=="3->5"]
        _title35 = r"A-site position relative to first codon of different amino acids " + \
        "in sample ({0}) determined from 3->5 with -{1} offset".format(_title,o35)
    
        fig = px.bar(data_frame=p35,x='position',y='counts',color='amino acid',animation_frame='condition', barmode='group')
        fig["layout"].pop("updatemenus") # optional, drop animation buttons
        fig.update_traces(dict(marker_line_width=0))
        xaxis_dict = dict(tickmode = 'linear',tick0 = 0,dtick = 3)
        fig.update_layout(title_text=_title35,xaxis=xaxis_dict,xaxis2=xaxis_dict)

        file_name_html = ofile.replace(".sam","_all_35_{0}.html".format(o35))
        fig.write_html(file_name_html)
        
        print("output written to {0}".format(file_name_html))


def refvariables(gns):
     # prepare arrays for faster access 

    AA_dict, _ = get_genetic_code()    

    NP_AA = np.array([])
    
    for aa in AA_dict.keys():    
        aa_list = np.array([l for l in gns[aa]],dtype=object)
    
        if(NP_AA.shape[0]==0):
            NP_AA = aa_list
        else:
            NP_AA = np.column_stack((NP_AA,aa_list))

    # create dictionary with numpy arrays for faster processing
    gene_info = {'ROISTART':gns.ROI_START.to_numpy('int64'),'ROISTOP': gns.ROI_STOP.to_numpy('int64'),
                'GENESTART':gns.start.to_numpy('int64'),'GENESTOP':gns.stop.to_numpy('int64'),
                'GENESTRANDS': gns.strand.apply(lambda x: 0 if x == '+' else 16),'DEF':gns,
                'NP_AA':NP_AA,'AA_DICT':AA_dict}
    
    return gene_info


def limit_memory(maxperc=0.8): 
    m = psutil.virtual_memory()
    maxsize = int(m.total*maxperc)
    soft, hard = resource.getrlimit(resource.RLIMIT_AS) 
    resource.setrlimit(resource.RLIMIT_AS, (maxsize, hard)) 


#%%
# https://stackoverflow.com/questions/41231678/obtaining-a-exclusive-lock-when-writing-to-an-hdf5-file

class SafeHDF5Store(pd.HDFStore):
    def __init__(self, *args, **kwargs):
        interval   = kwargs.pop('probe_interval', 1)
        # create lock file name
        self._lock = "%s.lock" % args[0]
        while True:
            try:
                # try to open lock file
                self._flock = os.open(self._lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                # continue if suceeded
                break
            except (IOError, OSError):
                # wait for a second                
                time.sleep(interval)
        pd.HDFStore.__init__(self, *args, **kwargs)

    def __exit__(self, *args, **kwargs):
        pd.HDFStore.__exit__(self, *args, **kwargs)
        # close lock file
        os.close(self._flock)
        # remove lock file
        os.remove(self._lock)        

    def write_hdf(f, key, df):    
        with SafeHDF5Store(f) as store:
            # use the put option here instead of df.to_hdf because of type case issues
            store.put(key,df)

#%%
def mythreadfunc_block_hdf(index, store_ptr, df_work, multi_dict_, d_normal, d_reverse, fasta_str, gns_dict, sema, offset35:int=12, offset53:int=14, verbose=True):        
    
    reg_out,nrm_out,rev_out = process_sam_vectorized(df_work,fasta_str,gns_dict, index, offset35,offset53,verbose)    
    if verbose:
        print("storing results for block {0}".format(index+1))    
    
    sema.release();
        
    with lock:

        if verbose:
            print("storing block to file {0}{1}".format(index,store_ptr))            
        _block ="/augmented/block_{0}".format(index)
        SafeHDF5Store.write_hdf(store_ptr,_block,reg_out)
        multi_dict_[index] = index
        d_normal[index]=nrm_out
        d_reverse[index]=rev_out        
        

def store_hdf(data, name:str, group:str, filename):
    # https://stackoverflow.com/questions/22922584/how-to-overwrite-array-inside-h5-file-using-h5py
    # Read/write, file must exist
    print("storing block data {0},{1},{2}".format(name,group,filename))
    with h5py.File(filename, "r+") as f:
        if not(group in f.keys()):
            gr = f.create_group(group)
        else:
            gr = f.get(group)        
            # create if not exising            
            if not(name in gr.keys()):
                gr.create_dataset(name,data)
            _f = gr[name]
            _f[...] = data
            

def get_hdf_items(filename,mainfld=""):    
    folders = []
    if exists(filename):
        try:
            with h5py.File(filename, "r") as f:
                if mainfld!="":
                    folders = [k for k in f.get(mainfld).keys() ]
                else:
                    folders = [k for k in f.keys() ]
        except:
            return []
    return folders

def get_hdf_item(filename,item=""):        
    if exists(filename):
        try:
            with h5py.File(filename, "r") as f:
                return f[item][()]
        except:
            return []
    return []

def get_hdf_attribute(filename,attr=""):
    if(exists(filename)):
        try:
            with h5py.File(filename,"r") as f:
                return f.attrs[attr]
        except:
            return "" 
    return "" 

def set_hdf_attribute(filename,val,attr=""):
    if(exists(filename)):
        try:
            with h5py.File(filename,"r+") as f:
                f.attrs[attr]=val
        except:
            pass
        

class PROJECT_INFO(IntEnum):
    TITLE = 0
    NBLOCKS = 1
    ARGS = 2
    STORE = 3
    OUTPATH = 4
    BASENAME = 5
    CPU = 6

#%% the main routine for processing
from tap import Tap # typed-argument-parser
import sys
interactive_mode = hasattr(sys, 'ps1')

if interactive_mode:
    #sys.argv = ['script','--sam','demo/filtered_Histag_Standard.sam','--gff', 'reference_files/wt-prspb-amyI-GFF3.gff','--fa',r'reference_files/WT-Prspb-amyI-FASTA.fa','--nc','10','--nb','100']
    #sys.argv = ['script','--sam','demo/demo_file.sam','--gff', 'reference_files/wt-prspb-amyI-GFF3.gff','--fa',r'reference_files/WT-Prspb-amyI-FASTA.fa','--nb','10']
    sys.argv = ['script','--sam','demo/demo_1.sam','--gff', 'reference_files/wt-prspb-amyI-GFF3.gff','--fa',r'reference_files/WT-Prspb-amyI-FASTA.fa','--orf','--nc','10','--verbose']

class myargs(Tap):
    
    sam: str  # the input sam(or bam) file
    gff: str # the file with gene definitions
    fa: str # the file with nucleotide sequences    
    nc: int=1 # the number of cores to use
    nb: int=1 # in how many parts the sam file is split
    mq:int=41 # the minimum mapping quality of a read
    #si:bool=False # save the intermediate results to (pkl) file
    o53:int=14 # the offset from 5->3 direction
    o35:int=12 # the offset from 3->5 direction
    title:str="" # the sample title
    log:bool=False # create a log file
    op:str="" # output folder, if not specified same as folder where sam resides
    ow:bool=False # overwrite existing output
    orf:bool=False # include ORF plots
    csv:bool=False # export to csv
    verbose:bool=False # show extra logging information
    lim:float=0.8 # set maximum memory limit, throws error if exceeded    
    
#%%
# semaphore from https://stackoverflow.com/questions/20886565/using-multiprocessing-process-with-a-maximum-number-of-simultaneous-processes
# add this line for processing in multiple threads
if __name__ ==  '__main__': 

    # define project dictionary for all settings
    project_data = {}
    
    args = myargs().parse_args()       
    
    if sys.platform == "linux" or sys.platform == "linux2":    
        limit_memory(args.lim)
 
    sam_file_name = args.sam    
    output_file_name = sam_file_name.lower().replace(".sam","_results.h5")
    sam_file_name_base = basename(args.sam)
    nr_blocks = args.nb


    dir_path = dirname(realpath(output_file_name))
    if args.op != "":
        dir_path = args.op 

        if not(exists(dir_path)):
            print("Error, the output path {0} does not exist, alter or create".format(dir_path))
            exit(1)

    output_file_name = basename(output_file_name)
    output_file_full = join(dir_path,output_file_name)

    print('processing file {0}'.format(sam_file_name))

    if args.log:    
        temp_stdout = sys.stdout
        logfile = join(dir_path,sam_file_name_base.lower().replace(".sam",".log"))
        print("logging to file {0}".format(logfile))
        sys.stdout = open(logfile,'wt')

    # remove the output file, otherwise it only grows
    if args.ow:
        print("Overwriting exsiting output")        
        os.remove(output_file_full)

    
    folders = get_hdf_items(output_file_full,"augmented")
    nblocks = get_hdf_attribute(output_file_full,"nblocks")

    maxcpu = args.nc             
    maxcpu = min(multiprocessing.cpu_count(),maxcpu)
    print('{0} cores/threads detected, running on {1} threads'.format(multiprocessing.cpu_count(),maxcpu))

    project_data[PROJECT_INFO.NBLOCKS]=nblocks
    project_data[PROJECT_INFO.STORE]=output_file_full
    project_data[PROJECT_INFO.OUTPATH]=dir_path
    project_data[PROJECT_INFO.BASENAME]=sam_file_name_base.lower()
    project_data[PROJECT_INFO.ARGS]=args    
    project_data[PROJECT_INFO.CPU]=maxcpu
    

    if not(exists(output_file_full)) or (len(folders)==0) or (args.ow) or (nblocks!=len(folders)):
        
            
        fasta_str = load_fasta_file(args.fa)
        gns = load_gene_defs(fasta_str, args.gff)                       
        
        concurrency = maxcpu
        total_task_num = nr_blocks

        if maxcpu>nr_blocks:
            print("nr of blocks {0} is smaller than number of threads {1} reserved, decreasing number of threads to number of blocks".format(nr_blocks,maxcpu))
            maxcpu = nr_blocks
                                    
        # create semaphore object 
        sema = multiprocessing.Semaphore(concurrency);
            
        multiprocessing_dict = multiprocessing.Manager().dict()
        dct_normal = multiprocessing.Manager().dict()
        dct_reverse = multiprocessing.Manager().dict()
        multiprocessing_loop_  = []

        # determine block sizes based on number of blocks argument

        file_lines = count_lines_in_file(sam_file_name)
        chunksize = (file_lines//nr_blocks)
        b = [i for i in range(0,file_lines,chunksize)]
        bb = b[1:]+[file_lines]
        
        nr_blocks = len(bb)

        # store number of blocks 
        project_data[PROJECT_INFO.NBLOCKS]=nr_blocks

        assert len(b) == nr_blocks, "Error defining block start positions"
        assert len(bb) == nr_blocks, "Error defining block end positions"
        assert(max(bb)==file_lines), "Error defining last block"

        nrheaders = len_headers_sam_file(sam_file_name)                
        dtSAM = dt.fread(sam_file_name,skip_to_line=nrheaders+1)
        
        # select variables of interest
        colindex = [0,1,2,3,4,9,10]
        colnames = ['x','dir','xx','left_pos','qual','sam_str','phred_scores']
                                                                        
        # make dictionary with reference variables for faster access 
        gene_info = refvariables(gns)       
            
        for i in range(nr_blocks):
            _df_sam = dtSAM[b[i]:bb[i],colindex].to_pandas()
            _df_sam.columns = colnames
            _df_sam = _df_sam[_df_sam.qual >= args.mq].copy()                                             
            # countdown available semaphores
            sema.acquire();
            _process = multiprocessing.Process(target=mythreadfunc_block_hdf, args=(i,output_file_full,_df_sam, multiprocessing_dict,dct_normal,dct_reverse,fasta_str,gene_info,sema,args.o35,args.o53,args.verbose))            
            multiprocessing_loop_.append(_process)
            _process.start()
        
        for process_  in multiprocessing_loop_ :
            process_.join()
                                                         
        blks2rerun = list(set([i for i in range(nr_blocks)]).difference(set(multiprocessing_dict.keys())))
              
        # rerun missed blocks in a single thread 
        
        for i in blks2rerun:
            _df_sam = dtSAM[b[i]:bb[i],colindex].to_pandas()
            _df_sam.columns = colnames
            _df_sam = _df_sam[_df_sam.qual >= args.mq].copy()                                             
            # create only a single semaphore
            sema = multiprocessing.Semaphore(1)
            sema.acquire();
            _process = multiprocessing.Process(target=mythreadfunc_block_hdf, args=(i,output_file_full,_df_sam, multiprocessing_dict,dct_normal,dct_reverse,fasta_str,gene_info,sema,args.o35,args.o53,args.verbose))                        
            _process.start()
            # wait for process to finish
            _process.join()
        
            if type(multiprocessing_dict[i]) is int:
                if multiprocessing_dict[i] == -1:
                    print("Not enough memory to process block {0}, exiting ... ".format(i))
                    sys.exit(1)                                                      
                
        # store nrblocks to attribute in hdf file
        set_hdf_attribute(output_file_full,nr_blocks,"nblocks")
        
        normal_reads = dct_normal[0].copy()
        reverse_read = dct_reverse[0].copy()
        
        for k in range(1,nr_blocks):            
            normal_reads += dct_normal[k]
            reverse_read += dct_reverse[k]            
    
    else:        
        print("reading from previous results {0}".format(output_file_full))

    make_plots_HDF(project_data)
                                       
    if args.log:
        sys.stdout = temp_stdout
    
    print("finished..")


# %%
