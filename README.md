# Python scripts to analyse a-site distributions in ribosomal profile data

These scripts are not optimized for speed and are likely to be improved over time. 


### genome = genome-build-accession NCBI_Assembly:GCF_000009045.1

## required packages ipython, pandas, numpy, biopython, plotly and typed-argument-parser

### RUN_A_SITE_ANALYSIS.PY, visualize the amino acid distributions of reads around the open reading frame of genes

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

### ALT_PREDICT_V2.PY, alternative pause prediction tool to Ribogalaxy

usage: alt_predict_v2.py --input_file INPUT_FILE --fasta_file FASTA_FILE [--gff_file GFF_FILE]
                         [--nt_upstream NT_UPSTREAM] [--nt_downstream NT_DOWNSTREAM] [--min_hits MIN_HITS]
                         [--down_offset DOWN_OFFSET] [--up_offset UP_OFFSET] [--min_qual MIN_QUAL] [--nlines NLINES]
                         [-h]

options:
  --input_file INPUT_FILE       (str, required) the input file (only sam files).  
  --fasta_file FASTA_FILE       (str, required) the fasta file to generate sequences for motif detection.  
  --gff_file GFF_FILE           (str, default=) optional, if supplied the peaks are annnotated with the location on the genome e.g. gene, ORF of gene.  
  --nt_upstream NT_UPSTREAM     (int, default=50) the number of positions upstream (towards 5') from the first nt of the A-site.  
  --nt_downstream NT_DOWNSTREAM (int, default=48) the number of positions downstream (towards 3') from the first nt of the A-site.  
  --min_hits MIN_HITS           (int, default=1) pauses should have at least min_hits reads.  
  --down_offset DOWN_OFFSET     (int, default=12) the offset when read from the 3' end.  
  --up_offset UP_OFFSET         (int, default=14) the offset when read from the 5' end.  
  --min_qual MIN_QUAL           (int, default=42) the minimum mapping quality of the reads.  
  --nlines NLINES               (int, default=1_000_000) the number of lines to be processed in a single block.  
  -h, --help                    show this help message and exit.  