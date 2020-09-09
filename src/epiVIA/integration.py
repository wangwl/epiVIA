#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-05-10 14:41:41
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import sys
import pybedtools
from epiVIA.ucsctracks import find_ucsctrack, find_trackhub, chunk_image, run_bigBedToBed
'''

'''

class Integration(object):
    """docstring for integration
    # Integration Attributes: Chrom ChrStart ChrEnd HotLevel VectorStart VectorEnd LTR  InsertOri   InsertGene  \
    GeneOri  ExonIntron TE TEFamily  TEClass NearestGene    NearestGeneDist  CellNumber  CellBarcodes
    # Example: hg38 chr7    114399221   114399221   5752    5752    LTR3    +   FOXP2   +   intron  None    None \
    Alu SINE    1   ATGATGCACAATAGCA    1   1   HostRead
    """

    def __init__(
        self,
        Genome='hg38',
        Chrom=None,
        ChrStart=None,
        ChrEnd=None,
        HotLevel=None,
        VectorStart=None,
        VectorEnd=None,
        LTR=None,
        InsertOri=None,
        InsertGene=None,
        GeneOri=None,
        ExonIntron=None,
        Enhancer=None,
        Elite=None,
        TEFamily=None,
        TEClass=None,
        UpstreamGene=None,
        UpstreamGeneDist=None,
        CellBarcodes=None,
        ReadInCell=None,
        CellNumber=None,
        ReadNumber=None,
        ChimeraSource=None,
        Quality = None
    ):
        super(Integration, self).__init__()
        self.Genome = Genome
        self.Chrom = Chrom
        self.ChrStart = ChrStart
        self.ChrEnd = ChrEnd
        self.VectorStart = VectorStart
        self.VectorEnd = VectorEnd
        self.LTR = LTR
        self.InsertOri = InsertOri
        self.InsertGene = InsertGene
        self.GeneOri = GeneOri
        self.ExonIntron = ExonIntron
        self.Enhancer = Enhancer
        self.Elite = Elite
        self.TEFamily = TEFamily
        self.TEClass = TEClass
        self.UpstreamGene = UpstreamGene
        self.UpstreamGeneDist = UpstreamGeneDist
        self.CellBarcodes = CellBarcodes
        self.ReadInCell = ReadInCell
        self.CellNumber = CellNumber
        self.ReadNumber = ReadNumber
        self.ChimeraSource = ChimeraSource
        self.Quality = Quality

    def annotate_TE(self):
        ChrEnd = self.ChrEnd + 1 if self.ChrEnd == self.ChrStart else self.ChrEnd
        results = find_ucsctrack('rmsk', self.Genome, self.Chrom, self.ChrStart, ChrEnd, [
                                 'repClass', 'repFamily'])

        try:
            repClass, repFamily = results[0]
        except:
            repClass = "None"
            repFamily = "None"

        self.TEFamily = repFamily
        self.TEClass = repClass
        return self

    def annotate_Gene(self, gbdb="http://hgdownload.soe.ucsc.edu/gbdb/"):
        ChrEnd = self.ChrEnd + 1 if self.ChrEnd == self.ChrStart else self.ChrEnd
        results = []
        # try:
        #     results = find_ucsctrack('knownGene', self.Genome, self.Chrom, self.ChrStart, ChrEnd, ['txStart', 'txEnd', 'strand', 'exonStarts', 'exonEnds', 'name'])
        #     # print(results)
        #     if len(results) > 0:
        #         GeneSt, GeneEd, GeneOri, Exon_Starts, Exon_Ends, GeneSymble = results[0]
        #         Exon_Starts = [int(x) for x in Exon_Starts.rstrip(",").split(",")]
        #         Exon_Ends = [int(x) for x in Exon_Ends.rstrip(",").split(",")]
        # except:
        results = run_bigBedToBed('knownGene', self.Genome, self.Chrom, self.ChrStart, ChrEnd, [1, 2, 5, 10, 11, 17], gbdb)
        if len(results) > 0:
            GeneSt, GeneEd, GeneOri, ExonLen, ExonSt, GeneSymble = results[0]
            GeneSt = int(GeneSt)
            # print(results[0])
            Exon_Lengths = [int(x) for x in ExonLen.rstrip(",").split(",")]
            Exon_Starts = [GeneSt + int(ExonSt.rstrip(",").split(",")[x]) for x in range(len(Exon_Lengths))]
            Exon_Ends = [Exon_Starts[x] + Exon_Lengths[x] for x in range(len(Exon_Lengths))]
            

            for x in range(len(Exon_Starts)):
                if self.ChrStart >= Exon_Starts[x] and self.ChrStart <= Exon_Ends[x]:
                    ExonIntron = "exon"
                elif x + 1 < len(Exon_Starts):
                    if self.ChrStart > Exon_Ends[x] and self.ChrStart <= Exon_Starts[x + 1]:
                        ExonIntron = "intron"
                    elif (self.ChrStart >= Exon_Starts[x] and self.ChrStart <= Exon_Ends[x]) and (self.ChrEnd >= Exon_Starts[x + 1] and self.ChrEnd <= Exon_Ends[x + 1]):
                        ExonIntron = 'junction'
            try:
                ExonIntron
            except:
                ExonIntron = "None"
        else:
            ExonIntron = "None"
            GeneSymble = "None"
            GeneOri = "None"

        self.ExonIntron = ExonIntron
        self.InsertGene = GeneSymble
        self.GeneOri = GeneOri
        return self

    def annotate_Enhancer(self, gbdb="http://hgdownload.soe.ucsc.edu/gbdb/"):
        ChrEnd = self.ChrEnd + 1 if self.ChrEnd == self.ChrStart else self.ChrEnd
        # try:
        #     results = find_ucsctrack('geneHancerRegElements', self.Genome, self.Chrom, self.ChrStart, ChrEnd, ['elementType', 'eliteness'])
        # except:
        results = run_bigBedToBed('geneHancerRegElements', self.Genome, self.Chrom, self.ChrStart, ChrEnd, [-2, -1], gbdb)
        
        try:
            elementType, eliteness = results[0]
        except:
            elementType = "None"
            eliteness = "None"

        self.Enhancer = elementType
        self.Elite = eliteness
        return self

    def annotate_UCSC_track(self, track):
        pass

    def nearest_gene(self):
        ChrEnd = self.ChrEnd + 1 if self.ChrEnd == self.ChrStart else self.ChrEnd
        integration_site_bed = pybedtools.BedTool("{chrom} {start} {end}".format(chrom=self.Chrom, start=self.ChrStart, end=ChrEnd), from_string=True)
        chrom_genes = find_ucsctrack('knownGene', self.Genome, self.Chrom)
        chrom_genes_list = []
        for gene in chrom_genes:
            '''
                    The result from find_ucsctrack can be either a list of dict or list of list, based on retrieve from bigBed or API.
                    1. Content of each gene as a dict:
                    {'alignID': 'uc033dmk.2',
                      'cdsEnd': 1670,
                      'cdsStart': 1670,
                      'chrom': 'chrM',
                      'exonCount': 1,
                      'exonEnds': '3229,',
                      'exonStarts': '1670,',
                      'name': 'ENST00000387347.2',
                      'proteinID': '',
                      'strand': '+',
                      'txEnd': 3229,
                      'txStart': 1670}
                    2. Content of each gene as a list:
                    [chrM, 1670, 3229, ENST00000387347.2, 3, +, 1670, 1670, 789624, 1, "1559,", "0,", ENST00000387347.2, none, none, "-1,", , MT-RNR2, , ]
                    '''
            if isinstance(gene, dict):
                chrom_genes_list.append("{}\t{}\t{}\t{}\t{}".format(gene['chrom'], gene['txStart'], gene['txEnd'], gene['name'], gene['strand']))
            if isinstance(gene, list):
                chrom_genes_list.append("{}\t{}\t{}\t{}\t{}".format(gene[0], gene[1], gene[2], gene[3], gene[5]))
        chrom_genes_string = "\n".join(chrom_genes_list)
        chrom_genes_bed = pybedtools.BedTool(chrom_genes_string, from_string=True)
        chrom_genes_bed = chrom_genes_bed.sort()
        upstream_closest = integration_site_bed.closest(chrom_genes_bed, d=True, D='ref', io=True, id=True)
        downstream_closest = integration_site_bed.closest(chrom_genes_bed, d=True, D='ref', io=True, iu=True)
        # The closest BedTool object will be in the format of:
        # chrM  1000    1300    chrM    576 647 ENST00000387314.1   +   -354
        for gene in upstream_closest:
            self.UpstreamGeneDist = gene[-1]
            self.UpstreamGeneOri = gene[-2]
            self.UpstreamGene = gene[-3]
        for gene in downstream_closest:
            self.DownstreamGeneDist = gene[-1]
            self.DownstreamGeneOri = gene[-2]
            self.DownstreamGene = gene[-3]
        return self

    def visualize_UCSC_track(self, imgdir):
        img_start = str(self.ChrStart - 10000)
        img_end = str(self.ChrEnd + 10000)
        position = "{Chrom}:{start}-{end}".format(
            Chrom=self.Chrom, start=img_start, end=img_end)
        params = {'db': self.Genome, 'position': position}
        chunk_image(params, imgdir=imgdir)

    def distance_to_known(self, konwnsites):
        pass

    def perturbation(self):
        pass

    def tostring(self):
        try:
            print_string = "\t".join([self.Genome, self.Chrom, str(self.ChrStart), str(
                self.ChrEnd), str(self.VectorStart), str(self.VectorEnd), self.LTR, self.InsertOri])
            header = "Genome\tChrom\tChrStart\tChrEnd\tVectorStart\tVectorEnd\tLTR\tInsertOri"
        except:
            print(self.Genome, self.Chrom, str(self.ChrStart), str(self.ChrEnd), str(self.VectorStart), str(self.VectorEnd), self.LTR, self.InsertOri)

        try:
            print_string += "\t" + self.InsertGene + \
                "\t" + self.GeneOri + "\t" + self.ExonIntron
            header += "\tInsertGene\tGeneOri\tExonIntron"
        except:
            print("Gene not annoated", file=sys.stderr)

        # try:
        #     print_string += "\t" + self.Enhancer + "\t" + self.Elite
        #     header += "\tEnhancer\tElite"
        # except:
        #     print("Enhancer not annoated", file=sys.stderr)

        try:
            print_string += "\t" + self.TEFamily + "\t" + self.TEClass
            header += "\tTEFamily\tTEClass"
        except:
            print("TE not annoated", file=sys.stderr)

        # try:
        #     print_string += "\t" + self.UpstreamGene + "\t" + \
        #         self.UpstreamGeneOri + "\t" + self.UpstreamGeneDist
        #     header += "\tUpstreamGene\tUpstreamGeneOri\tUpstreamGeneDist"
        # except:
        #     print("UpstreamGene not annoated", file=sys.stderr)

        # try:
        #     print_string += "\t" + self.DownstreamGene + "\t" + \
        #         self.DownstreamGeneOri + "\t" + self.DownstreamGeneDist
        #     header += "\tDownstreamGene\tDownstreamGeneOri\tDownstreamGeneDist"
        # except:
        #     print("DownstreamGene not annoated", file=sys.stderr)

        try:
            print_string += "\t" + str(self.CellNumber) + "\t" + ",".join(self.CellBarcodes) + \
                "\t" + ",".join([str(x) for x in self.ReadInCell]
                                ) + "\t" + str(self.ReadNumber)
            header += "\tCellNumber\tCellBarcodes\tReadInCell\tReadNumber"
        except:
            print("Cell and Read not calculated", file=sys.stderr)

        try:
            print_string += "\t" + ",".join(self.ChimeraSource)
            header += "\tChimeraSource"
        except:
            print("ChimeraSource not defined", file=sys.stderr)

        try:
            print_string += "\t{}".format(self.Quality)
            header += "\tQuality"
        except:
            pass

        print_string += "\n"
        header += "\n"
        return print_string, header


class HotSpot(object):
    """docstring for HotSpot"""

    def __init__(self, arg):
        super(HotSpot, self).__init__()
        self.arg = arg
