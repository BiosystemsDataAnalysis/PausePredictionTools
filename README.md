# A python script to visualize the amino acid distributions of reads around the open reading frame

### genome = genome-build-accession NCBI_Assembly:GCF_000009045.1

## required packages ipython, pandas, numpy, biopython, plotly and typed-argument-parser

usage: run_a_site_analysis.py --sam SAM --gff GFF --fa FA [--nc NC] [--nb NB] [--mq MQ] [--o53 O53] [--o35 O35] [--title TITLE] [--log] [--op OP] [--orf] [-h]  

options:  
  --sam SAM      (str, required) the input sam(or bam) file  
  --gff GFF      (str, required) the file with gene definitions  
  --fa FA        (str, required) the file with nucleotide sequences  
  --nc NC        (int, default=1) the number of cores to use  
  --nb NB        (int, default=1) the number of blocks to use (NC<=NB)
  --mq MQ        (int, default=41) the minimum mapping quality of a read  
  --o53 O53      (int, default=14) the offset from 5->3 direction  
  --o35 O35      (int, default=12) the offset from 3->5 direction  
  --title TITLE  (str, default=) the sample title  
  --log          (bool, default=False) write output to a log file  
  --op OP        (str, default=) output folder, if not specified same as folder where sam resides  
  --orf          (bool, default=False) generate ORF plots  
  -h, --help     show this help message and exit  