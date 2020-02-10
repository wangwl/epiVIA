## Description of epiVIA
Vector Integration Analysis with epigenomic assays.

![epiVIA pipeline](images/Figure1.pdf)

### (1) Identify Chimeric reads from bam file

EpiVIA parses each read and identifies the chimeric fragments, which are further classified into three different categories based on how the read pair is aligned to the combined reference genome: (1) “pair-chimeric” are the fragments with one read mapped to the host genome and the other mapped to the provirus (or vector) genome; (2) “host-chimeric” are the fragments that both ends mapped to the host genome, with a small soft-clipped fragment at one end which can exactly match the start or end of LTR sequence of provirus; (3) “viral-chimeric” are the fragments that are properly pair mapped to either end of the provirus genome, with a soft-clipped fragment at the end which can be mapped to the host genome with a non-zero mapping quality. Since the real chimeric fragment won’t be soft-clipped at 3’ end of 5’ LTR and 5’ end of 3’ LTR, while these cases can present in the alignment results, we corrected the aligned position for these fragments. 

### (2) Identify integration site with chimeric reads

Although there is pair rescue during the alignment of BWA, we still found some pair-chimeric fragments which have soft-clipped ends from the other genome. Hence, we intend to identify the integration site with the soft-clipped read. If the read mapped to host genome is a soft-clipped alignment, we can determine the precise integration site if the clipped oligonucleotide exactly match either end of the LTR of provirus. In the other case, if the viral read that is aligned to the end of LTR is soft-clipped, we try to search the clipped fragment in the context of up or down 200 bp where the host read aligned, and an exact match can tell us the precise integration site. For host-chimeric and viral-chimeric fragments, we identify the integration site while classifying these fragments. These integration sites are further annotated with genes, transposable elements and enhancer annotation from UCSC genome browser database.
 

## Install

### install with pip
> pip install epiVIA

### intall from github page
> git clone https://github.com/wangwl/epiVIA
>
> cd epiVIA
>
> python setup.py install


## Output
This pipeline has multiple output files for the downstream analysis of integration sites, cells with integration sites and with proviral reads. The following tables are detailed descriptions of the output files and their file formats.

### Result files
| filename                 | description                                                          |
| ------------------------ |:-------------------------------------------------------------------- |
| Integration_sites.list   | Identified integration sites                                         |
| barcode_integration.list | Ingretation sites sorted by cell barcoeds                            |
| Cell_summary.txt         | Number of chimeric and viral reads in each cell                      |
| Vector_coverage.csv      | Coverage of viral reads at each positioin along the proviral genome  |

#### Integration_sites.list

| Column name  |  Description                                                                                    |
| ------------ | ----------------------------------------------------------------------------------------------- |
| Genome       | Reference genome version (eg: hg19, hg38, mm9, mm10)                                            |
| Chrom        | Chromosome of the integration site                                                              |
| ChrStart     | Start position in host genome where the chimeric reads aligned, or the integration site         |
| ChrEnd       | End Position in host genome where the chimeric reads aligned, or the integration site           |
| VectorStart  | Start posision in proviral genome where the chimeric reads aligned                              |
| VectorEnd    | End position in proviral genome where the chimeric reads aligned                                |
| LTR          | The chimeric reads is from 5' LTR or 3' LTR                                                     |
| InsertOri    | Integration orietation of the proviral genome based on the reference sequence                   |
| InsertGene   | The annotated gene at the integration site of the host genome                                   |
| GeneOri      | Orientation of the gene                                                                         |
| ExonIntron   | The integration is in intron or exon region of the gene                                         |
| Enhancer     | Enhancer or Promoter information at the integration site                                        |
| Elite        | If the integration site locates in an enhancer based on UCSC annotation, whether it's Elite     |
| TEFamily     | Family of transposable elements annotation at the integration site                              |
| TEClass      | Class of transposable elements annotation at the integration site                               |
| CellNumber   | Number of cells have this integration site                                                      |
| CellBarcodes | Cell barcodes that have this integration site                                                   |
| ReadInCell   | Number of reads in each cell that support this integration site                                 |
| ReadNumber   | Total number of reads that support this integration site                                        |
| Quality      | The quality of the integration site, adapted from the mapping quality of BWA                    |

#### barcode_integration.list
Each cell that have integration sites. If a cell have multiple integration sites, ".1", ".2", ".3", ..., were appended to the cell barcode.

| Column name  |  Description                                                                                    |
| ------------ | ----------------------------------------------------------------------------------------------- |
| CellBarcode  | Cell Barcode of the cell                                                                        |
| Chrom        | Chromosome of the integration site                                                              |
| ChrStart     | Start position in host genome where the chimeric reads aligned, or the integration site         |
| ChrEnd       | End Position in host genome where the chimeric reads aligned, or the integration site           |
| VectorStart  | Start posision in proviral genome where the chimeric reads aligned                              |
| VectorEnd    | End position in proviral genome where the chimeric reads aligned                                |
| InsertGene   | The annotated gene at the integration site of the host genome                                   |
| ExonIntron   | The integration is in intron or exon region of the gene                                         |
| Enhancer     | Enhancer or Promoter information at the integration site                                        |
| TEFamily     | Family of transposable elements annotation at the integration site                              |
| TEClass      | Class of transposable elements annotation at the integration site                               |
| Quality      | The quality of the integration site, adapted from the mapping quality of BWA                    |

#### Cell_summary.txt

| Column name  |  Description                                                                                    |
| ------------ | ----------------------------------------------------------------------------------------------- |
| CellBarcode  | Cell Barcode of the cell                                                                        |
| ChimeraSum   | Total number of chimeric reads from the cell                                                    |
| Pair         | Number of 'Pair' chimeric reads (case a) from the cell                                          |
| HostRead     | Number of 'HostRead' chimeric reads (case b) from the cell                                      |
| VectorRead   | Number of 'VectorRead' chimeric reads (case c) from the cell                                    |
| integration  | Number of integration sites in the cell                                                         |
| VectorSum    | Total number of reads from the integrated proviral genome                                       |
| noLTR        | Number of proviral reads that is not from LTR regions                                           |
| LTR          | Number of proviral reads that is from LTR regions                                               |



### Temporary files

epiVIA has options of "--candidate_bam" and "--chimera" for writing the bam records used for identifying chimeric reads and those that are classified as chimeric between host genome and viral genome.


## Example
