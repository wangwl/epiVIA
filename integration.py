#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-05-10 14:41:41
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
from ucsctracks import find_ucsctrack, find_trackhub, chunk_image

# Integration Attributes: Chrom ChrStart ChrEnd HotLevel VectorStart VectorEnd LTR	InsertOri	InsertGene	GeneOri	 ExonIntron	TE TEFamily  TEClass NearestGene	NearestGeneID	NearestGeneDist	 CellNumber  CellBarcodes
class IntegrationSite(object):
	"""docstring for integration"""
	def __init__(
			self,
			Genome=None,
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
			CellNumber=None,
			CellBarcodes=None,
			ReadInCell=None,
			ReadNumber=None
	):
		super(Integration, self).__init__()
		self.Genome = Genome,
		self.Chrom = Chrom,
		self.ChrStart=ChrStart,
		self.ChrEnd=ChrEnd,
		self.VectorStart=VectorStart,
		self.VectorEnd=VectorEnd,
		self.LTR=LTR,
		self.InsertOri=InsertOri,
		self.InsertGene=InsertGene,
		self.GeneOri=GeneOri,
		self.ExonIntron=ExonIntron,
		self.Enhancer=Enhancer,
		self.Elite=Elite,
		self.TEFamily=TEFamily,
		self.TEClass=TEClass,
		self.NearestGene=NearestGene,
		self.NearestGeneDist=NearestGeneDist,
		self.CellNumber=CellNumber,
		self.CellBarcodes=CellBarcodes

	def annotate_TE(self, tefile):
		results = find_ucsctrack(self.Genome, self.Chrom, self.ChrStart, self.ChrEnd, 'rmsk', ['repClass', 'repFamily'])
		repClass, repFamily = results[0]
		self.TEFamily = repFamily
		self.TEClass = repClass
		return self

	def annotate_Gene(self):
		results = find_ucsctrack(self.Genome, self.Chrom, self.ChrStart, self.ChrEnd, 'knownGene', [5, 17])
		GeneOri, GeneSymble = results[0][0]
		self.InsertGene = GeneSymble
		self.GeneOri = GeneOri
		return self

	def annotate_Enhancer():
		results = find_ucsctrack(self.Genome, self.Chrom, self.ChrStart, self.ChrEnd, 'geneHancerRegElements', ['elementType', 'eliteness'])
		elementType, eliteness = results[0]
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
		

	def tostring(self):


		
class HotSpot(object):
	"""docstring for HotSpot"""
	def __init__(self, arg):
		super(HotSpot, self).__init__()
		self.arg = arg
		
	