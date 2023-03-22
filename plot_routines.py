#%% import some modules


# from weblogo import *
import io
# from numpy import concatenate, result_type
from os.path import exists

# import seaborn as sns

from IPython.display import Markdown as md
from IPython.display import *

#import ipywidgets as widgets
#from ipywidgets import *

from collections import Counter
#import random
#from matplotlib.lines import Line2D


import plotly.express as px
import plotly.graph_objects as go
#import plotly

from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

#from ipywidgets import interact, interactive, fixed, interact_manual

# import timeit
from map_functions import getAsiteData53,getAsiteData35,get_genetic_code
# from itertools import chain

import pandas as pd 
import numpy as np
from tap import Tap 

# _offset53 = 14
# _offset35 = 12 # will be made negative in the subfunction

def detpos(vecdata):
      
    offset53=_offset53
    offset35=_offset35

    # print(offset53,offset35)

    strand = vecdata[0]
    start_read = vecdata[1]
    read_len = vecdata[2]
    position_I = vecdata[3]
    
    stop_read = start_read + (read_len-1)
    if strand==16:
        stop_read = start_read - (read_len - 1)  

    # print(vecdata)   
    result_53 = [getAsiteData53(start_read,p,strand,offset53) for p in position_I]
    result_35 = [getAsiteData35(stop_read,p,strand,offset35) for p in position_I]

    # result_53 = [getAsiteData53(p,start_read,strand,offset53) for p in position_I]
    # result_35 = [getAsiteData35(p,stop_read,strand,offset35) for p in position_I]
  

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
    return pd.DataFrame(df_ret,columns=['frequency'])



def createAsiteDistributionPlotData(df_in:pd.DataFrame):

    # create filters
    genes_summary_ = df_in.groupby('gene').agg(['count','min'])['gene_len']      
    
    
    #genes_.name = 'frequency'
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




def add_aa_scores(dfin:pd.DataFrame,aa:str='I'):
    
    # copy I data vectors and determine distances to reads and I position, in detpos the offsets are determined 
    # they can be set by the global variables _offset53 and _offset35
    cc = dfin[['dir','begin_read','read_len',aa]].apply(lambda x:detpos(x),axis=1)
    ccc=cc.apply(pd.Series)
    
    col35_, col53_= aa+"35_"+str(_offset35), aa+"35_"+str(_offset53)
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




#%%

import sys
interactive_mode = hasattr(sys, 'ps1')


if interactive_mode:
 sys.argv = ['script.py', '--input_file', 'output\\1X_PEG-Spin_filtered_SAM.pkl', '--sample_type','4']


class  myargs(Tap):
    input_file:str # input file name
    sample_type:int # sample type 1,2 or 3

sample_types = {1:'16 h chloramphenicol',2:'mupirocin',3:'in silico',4:'other'}

_args = myargs().parse_args()


# load the data with calculated distances from reads to isoleucine first nucleotide


kys_,cdns_ = get_genetic_code()


file_name = _args.input_file
# file_name = 'is_update_v2_1807.pkl'
if _args.sample_type==4:
    _sample = file_name.strip("/output").strip(".pkl")
else:
    _sample = sample_types[_args.sample_type]

#file_name = 'MUP_FK_2906_amyI_25_wIdist_13to15_9to12.pkl'


output_file_name = file_name.strip(".pkl")+"_with_AAs.pkl"

if exists(output_file_name):
    dfm = pd.read_pickle(output_file_name)
else:
    dfw = pd.read_pickle(file_name)
    # focus on mapped to known genes only
    dfmapped_on_genes = dfw[~(dfw.gene=="")].copy()
    dfm = addORFinfo(dfmapped_on_genes)

    for k_ in kys_:
        print("adding distances for AA "+k_)
        dfm = add_aa_scores(dfm,k_)

    dfm.to_pickle(output_file_name)


# to_remove = ['rrn','trn','ssrA']

# _flt = ~(dfm.gene.str.startswith('rrn')) & ~(dfm.gene.str.startswith('trn')) & ~(dfm.gene.str.startswith('ssrA'))
# dfm_flt = dfm[_flt]


# %% create read distributions for different settings
file_name_pickle = file_name.strip('.pkl')+"plot_data_14_12.pkl"

if exists(file_name_pickle):
    df_14_12 = pd.read_pickle(file_name_pickle)
else:
    df_14_12 = createAsiteDistributionPlotData(dfm)
    df_14_12.to_pickle(file_name_pickle)

#df_14_12_flt = createAsiteDistributionPlotData(dfm_flt)

#%%

work_data = df_14_12

plot_data = work_data[(work_data.position>=-12) & (work_data.position<=14)]


#%%



# _sample = 'mupirocin'
# _sample = '16 h chloramphenicol'
# _sample = 'in silico'

_title = r"A-site position relative to first codon of different amino acids " + \
    "in {0} sample determined with +14/-12 offsets".format(_sample)

fig = px.bar(data_frame=plot_data,x='position',y='counts',color='amino acid',facet_col="serie",animation_frame='condition', barmode='group')
fig["layout"].pop("updatemenus") # optional, drop animation buttons
fig.update_traces(dict(marker_line_width=0))
xaxis_dict = dict(tickmode = 'linear',tick0 = 0,dtick = 3)
fig.update_layout(title_text=_title,xaxis=xaxis_dict,xaxis2=xaxis_dict)

fig.show();
file_name_html = file_name.strip(".pkl")+"_all_14_12_wo_sel.html"
fig.write_html(file_name_html)

#%%
df_14_12_I = plot_data[plot_data['amino acid'].isin(['I'])]


_title = r"A-site position relative to first nucleotide of nearest Isoleucine codon " + \
    "in {0} sample determined with +14/-12 offsets".format(_sample)

fig = px.bar(data_frame=df_14_12_I,x='position',y='counts',facet_col="serie",animation_frame='condition', barmode='group')
fig["layout"].pop("updatemenus") # optional, drop animation buttons
fig.update_traces(dict(marker_line_width=0))
xaxis_dict = dict(tickmode = 'linear',tick0 = 0,dtick = 3)
fig.update_layout(title_text=_title,xaxis=xaxis_dict,xaxis2=xaxis_dict)



fig.show();
file_name_html = file_name.strip(".pkl")+"_I_14_12_wo_sel.html"
fig.write_html(file_name_html)

#%%
df_14_12_total = plot_data[plot_data['amino acid'].isin(['combined'])]


_title = r"A-site position relative to first nucleotide of nearest amino acid (any) codon " + \
    "in {0} sample determined with +14/-12 offsets".format(_sample)

fig = px.bar(data_frame=df_14_12_total,x='position',y='counts',facet_col="serie",animation_frame='condition', barmode='group')
fig["layout"].pop("updatemenus") # optional, drop animation buttons
fig.update_traces(dict(marker_line_width=0))
xaxis_dict = dict(tickmode = 'linear',tick0 = 0,dtick = 3)
fig.update_layout(title_text=_title,xaxis=xaxis_dict,xaxis2=xaxis_dict)

fig.show();
file_name_html = file_name.strip(".pkl")+"_total_14_12_wo_sel.html"
fig.write_html(file_name_html)
#%%

# #_sample = 'mupirocin'
# _sample = '16 h chloramphenicol'
# # _sample = 'in silico'


# _title = r"A-site position relative to first codon of amino acid " + \
#     "on {0} sample determined with +14/-12 offsets".format(_sample)

# fig = px.bar(data_frame=df_14_12,x='position',y='counts',color='amino acid',facet_col="serie",animation_frame='condition', barmode='group')
# fig["layout"].pop("updatemenus") # optional, drop animation buttons
# fig.update_layout(title_text=_title)
# fig.show();
# file_name_html = file_name.strip(".pkl")+"_all_14_12.html"
# fig.write_html(file_name_html)




# # %% combine results of CM and MUP

# df_14_12_mup = pd.read_pickle('Mup_0407plot_data_14_12.pkl')
# df_14_12_cm = pd.read_pickle('16_CM_0407plot_data_14_12.pkl')

# df_14_12_cm['sample']='chloramphenicol'
# df_14_12_mup['sample']='mupirocin'

# df_plot = pd.concat([df_14_12_cm,df_14_12_mup],axis=0)
# # %%

# _title = r"A-site position from first nt of iso-leucine codon over subset of genes (depending on condition) <br>" + \
#     " on mupirocin and chloramphenicol sample determined with +14(5'->3')/-12(3'->5') offsets"

# fig = px.bar(data_frame=df_plot,x='position',y='counts',color='sample',facet_col="serie",animation_frame='condition', barmode='group')
# fig["layout"].pop("updatemenus") # optional, drop animation buttons
# fig.update_layout(
#     title=dict(text=_title,x=0.5,y=0.95, xanchor='center',yanchor= 'top'),
#     xaxis_title="position",
#     yaxis_title="counts",
#     legend_title="Sample",
#     font=dict(      
#         size=24,        
#     )
# )

# fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
# fig.update_yaxes(showline=False)
# #, linewidth=2, linecolor='black')

# # fig.update_layout(title_text=_title)
# fig.show();


# fig.update_layout(margin=dict(l=50, r=50, t=200, b=50),
#     paper_bgcolor="White",    
#     plot_bgcolor='rgb(255,255,255)')

# # file_name_html = file_name.strip(".pkl")+"14_11.html"
# fig.write_html("MUP_CM_14_12.html")
# # %%
