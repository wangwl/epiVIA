## Description of epiVIA
Vector Integration Analysis with epigenomic assays.

### (1) Identify Chimeric reads from bam file

EpiVIA parses each read and identifies the chimeric fragments, which are further classified into three different categories based on how the read pair is aligned to the combined reference genome: (1) “pair-chimeric” are the fragments with one read mapped to the host genome and the other mapped to the provirus (or vector) genome; (2) “host-chimeric” are the fragments that both ends mapped to the host genome, with a small soft-clipped fragment at one end which can exactly match the start or end of LTR sequence of provirus; (3) “viral-chimeric” are the fragments that are properly pair mapped to either end of the provirus genome, with a soft-clipped fragment at the end which can be mapped to the host genome with a non-zero mapping quality. Since the real chimeric fragment won’t be soft-clipped at 3’ end of 5’ LTR and 5’ end of 3’ LTR, while these cases can present in the alignment results, we corrected the aligned position for these fragments. 

### (2) Identify integration site with chimeric reads

Although there is pair rescue during the alignment of BWA, we still found some pair-chimeric fragments which have soft-clipped ends from the other genome. Hence, we intend to identify the integration site with the soft-clipped read. If the read mapped to host genome is a soft-clipped alignment, we can determine the precise integration site if the clipped oligonucleotide exactly match either end of the LTR of provirus. In the other case, if the viral read that is aligned to the end of LTR is soft-clipped, we try to search the clipped fragment in the context of up or down 200 bp where the host read aligned, and an exact match can tell us the precise integration site. For host-chimeric and viral-chimeric fragments, we identify the integration site while classifying these fragments. These integration sites are further annotated with genes, transposable elements and enhancer annotation from UCSC genome browser database.
 

## Install

### install with pip
pip install epiVIA

### intall from github page
git clone https://github.com/wangwl/epiVIA

cd epiVIA

python setup.py install

## Usage
usage: epiVIA [-h] [--Provirus PROVIRUS] [--Vector VECTOR]
              [--ltr5_start LTR5_START] --ltr5_end LTR5_END
              [--ltr3_start LTR3_START] [--ltr3_end LTR3_END] --HostIndex
              HOSTINDEX [--LTRClipLen LTRCLIPLEN] [--HostClipLen HOSTCLIPLEN]
              --Host2bit HOST2BIT [--tempdir TEMPDIR] [--quality QUALITY]
              [--candidate_bam CANDIDATE_BAM] [--chimera_bam CHIMERA_BAM]
              [--gbdb GBDB] [--genome GENOME]
              bamfile outdir

Find chimera in alignment result with combined host and vector sequence as
reference.

positional arguments:
  bamfile
  outdir                A directory where the result files (stat on cells and
                        integrations) written

optional arguments:
  -h, --help            show this help message and exit
  
  --Provirus PROVIRUS   Fasta file of provirus sequence that is integrated in the host genome
                        
  --Vector VECTOR       Fasta file of vector sequence
  
  --ltr5_start LTR5_START
                        Start position of LTR5, required if using vector
                        
  --ltr5_end LTR5_END   End position of LTR5
  
  --ltr3_start LTR3_START
                        Start position of LTR3, required if using vector
                        
  --ltr3_end LTR3_END   End position of LTR3, required if using vector
  
  --HostIndex HOSTINDEX
                        The BWA Index prefix of Host genome
                        
  --LTRClipLen LTRCLIPLEN
                        The soft clipped read length to search in LTR
                        
  --HostClipLen HOSTCLIPLEN
                        The soft clipped length to search in host genome
                        
  --Host2bit HOST2BIT   The 2bit (UCSC) file of host genome sequence
  
  --tempdir TEMPDIR     A directory for temporary files generated in the
                        pipeline
                        
  --quality QUALITY     Set the MAPQ for the host read that locate the
                        integration site
                        
  --candidate_bam CANDIDATE_BAM
                        A bamfile to write the suspicious reads used for
                        identifying chimeras, saving time if rerun the
                        pipeline
                        
  --chimera_bam CHIMERA_BAM
                        A bamfile to write the records that are classified as
                        chimeras
                        
  --gbdb GBDB           gbdb directory of UCSC, default:
                        http://hgdownload.soe.ucsc.edu/gbdb/.
                        
  --genome GENOME       Genome name in UCSC, used for annotation of the
                        integration sites

## Output

EpiVIA reports the results of integration events both from every integration site and every single cell. In this case, we are able to identify cells that bear the same integration site which could possibly result from clonal expansion, and cells that have multiple integration sites, which might be caused by doublets in the single-cell ATAC experiment. Specifically, we found there are multiple chimeric reads aligned to the poly(G) region in ‘LINC00486’, which additionally has many chimeric reads with other chromosomes of the host genome (data not shown), in both bulk and single-cell ATAC-seq data. Therefore, the identified integration sites in this region were excluded from the results.

We also calculated the coverage of the provirus sequence in each cell in EpiVIA. Specifically, for the fragments that are properly mapped to the LTR region of the provirus genome, which can be alternatively aligned to the other LTR, we shifted the aligned position to the 3’ LTR if the pair is reported to aligned to the 5’ LTR in the bam file. Therefore, we are able to remove the PCR duplicates in LTR region. While calculating the coverage of the LTR region, we divided the coverage to both ends.

### example
