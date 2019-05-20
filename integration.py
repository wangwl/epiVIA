#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-05-10 14:41:41
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import pysnooper
from ucsctracks import find_ucsctrack, find_trackhub, chunk_image

# Integration Attributes: Chrom ChrStart ChrEnd HotLevel VectorStart VectorEnd LTR	InsertOri	InsertGene	GeneOri	 ExonIntron	TE TEFamily  TEClass NearestGene	NearestGeneDist	 CellNumber  CellBarcodes
class Integration(object):
	"""docstring for integration"""
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
			NearestGene=None,
			NearestGeneDist=None,
			CellBarcodes=None,
			ReadInCell=None,
			CellNumber=None,
			ReadNumber=None,
			ChimeraSource=None
	):
		super(Integration, self).__init__()
		self.Genome = Genome
		self.Chrom = Chrom
		self.ChrStart=ChrStart
		self.ChrEnd=ChrEnd
		self.VectorStart=VectorStart
		self.VectorEnd=VectorEnd
		self.LTR=LTR
		self.InsertOri=InsertOri
		self.InsertGene=InsertGene
		self.GeneOri=GeneOri
		self.ExonIntron=ExonIntron
		self.Enhancer=Enhancer
		self.Elite=Elite
		self.TEFamily=TEFamily
		self.TEClass=TEClass
		self.NearestGene=NearestGene
		self.NearestGeneDist=NearestGeneDist
		self.CellBarcodes=CellBarcodes
		self.ReadInCell=ReadInCell
		self.CellNumber=CellNumber
		self.ReadNumber=ReadNumber
		self.ChimeraSource=ChimeraSource

	def annotate_TE(self):
		ChrEnd = self.ChrEnd + 1 if self.ChrEnd == self.ChrStart else self.ChrEnd
		results = find_ucsctrack(self.Genome, self.Chrom, self.ChrStart, ChrEnd, 'rmsk', ['repClass', 'repFamily'])
		
		try:
			repClass, repFamily = results[0]
		except:
			repClass="None"
			repFamily="None"
		
		self.TEFamily = repFamily
		self.TEClass = repClass
		return self

	def annotate_Gene(self):
		ChrEnd = self.ChrEnd + 1 if self.ChrEnd == self.ChrStart else self.ChrEnd
		results = find_ucsctrack(self.Genome, self.Chrom, self.ChrStart, ChrEnd, 'knownGene', [1, 2, 5, 10, 11, 17])
		if len(results) > 0:
			GeneSt, GeneEd, GeneOri, ExonLen, ExonSt, GeneSymble = results[0]
			Exon_Lengths = [str(x) for x in ExonLen.rstrip(",").split(",")]
			Exon_Starts = [GeneSt + str(ExonSt.rstrip(",").split(",")[x]) for x in xrange(0,len(Exon_Lengths))]
			Exon_Ends = [GeneSt + Exon_Starts[x] + Exon_Lengths[x] for x in xrange(0,len(Exon_Lengths))]
			for x in xrange(0, len(Exon_Starts)):
				if self.ChrStart > Exon_Starts[x] and self.ChrEnd <= Exon_Ends[x]:
					ExonIntron = "exon"
				elif x +1 < len(Exon_Starts) and self.ChrStart > Exon_Ends[x] and self.ChrStart <= Exon_Starts[x+1]:
					ExonIntron = "intron"
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

	def annotate_Enhancer(self):
		ChrEnd = self.ChrEnd + 1 if self.ChrEnd == self.ChrStart else self.ChrEnd
		results = find_ucsctrack(self.Genome, self.Chrom, self.ChrStart, ChrEnd, 'geneHancerRegElements', ['elementType', 'eliteness'])
		try:
			elementType, eliteness = results[0]
		except:
			elementType = "None"
			eliteness="None"
		
		self.Enhancer = elementType
		self.Elite = eliteness
		return self

	def annotate_UCSC_track(self, track):
		pass

	def nearest_gene(self):
		pass

	def visualize_UCSC_track(self, imgdir):
		img_start = str(self.ChrStart - 10000)
		img_end = str(self.ChrEnd + 10000)
		position = "{Chrom}:{start}-{end}".format(Chrom=self.Chrom, start=img_start, end=img_end)
		params = {'db':self.Genome, 'position':position}
		chunk_image(params, imgdir=imgdir)

	def is_known(self):
		pass

	def tostring(self):
		try:
			print_string = "\t".join([self.Genome, self.Chrom, str(self.ChrStart), str(self.ChrEnd), str(self.VectorStart), str(self.VectorEnd), self.LTR, self.InsertOri])
			header = "Genome\tChrom\tChrStart\tChrEnd\tVectorStart\tVectorEnd\tLTR\tInsertOri"
		except:
			print self.Genome, self.Chrom, str(self.ChrStart), str(self.ChrEnd), str(self.VectorStart), str(self.VectorEnd), self.LTR, self.InsertOri
		
		try:
			print_string += "\t" + self.InsertGene + "\t" + self.GeneOri + "\t" + self.ExonIntron
			header += "\tInsertGene\tGeneOri\tExonIntron"
		except:
			print "Gene not annoated"
		
		try:
			print_string += "\t" + self.Enhancer + "\t" + self.Elite
			header += "\tEnhancer\tElite"
		except:
			print "Enhancer not annoated"
		
		try:
			print_string += "\t" + self.TEFamily + "\t" + self.TEClass
			header += "\tTEFamily\tTEClass"
		except:
			print "TE not annoated"
		
		try:
			print_string += "\t" + self.NearestGene + "\t" + str(self.NearestGeneDist)
			header += "\tNearestGene\tNearestGeneDist"
		except:
			print "NearestGene not defined"
		
		try:
			print_string += "\t" + str(self.CellNumber) + "\t" + ",".join(self.CellBarcodes) + "\t" + ",".join([str(x) for x in self.ReadInCell]) + "\t" + str(self.ReadNumber)
			header += "\tCellNumber\tCellBarcodes\tReadInCell\tReadNumber"
		except:
			print "Cell and Read not calculated"
		
		try:
			print_string += "\t" + ",".join(self.ChimeraSource)
			header += "\tChimeraSource"
		except:
			print "ChimeraSource not defined"
		
		print_string += "\n"
		header += "\n"
		return print_string, header
		
class HotSpot(object):
	"""docstring for HotSpot"""
	def __init__(self, arg):
		super(HotSpot, self).__init__()
		self.arg = arg
		
	