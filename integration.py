#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-05-10 14:41:41
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os

# Integration Attributes: Chrom ChrStart ChrEnd HotLevel VectorStart VectorEnd LTR	InsertOri	InsertGene	GeneOri	 ExonIntron	TE TEFamily  TEClass NearestGene	NearestGeneID	NearestGeneDist	 CellNumber  CellBarcodes
class IntegrationSite(object):
	"""docstring for integration"""
	def __init__(
			self,
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
			TE=None,
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
		self.chrom = Chrom,
		self.ChrStart=ChrStart,
		self.ChrEnd=ChrEnd,
		self.VectorStart=VectorStart,
		self.VectorEnd=VectorEnd,
		self.LTR=LTR,
		self.InsertOri=InsertOri,
		self.InsertGene=InsertGene,
		self.GeneOri=GeneOri,
		self.ExonIntron=ExonIntron,
		self.TE=TE,
		self.TEFamily=TEFamily,
		self.TEClass=TEClass,
		self.NearestGene=NearestGene,
		self.NearestGeneDist=NearestGeneDist,
		self.CellNumber=CellNumber,
		self.CellBarcodes=CellBarcodes

	def get_reads_num(self):
		pass

	def get_cell_barcodes(self):
		pass

	def get_metrics(self):
		pass

	def _annotate_TE(self, tefile):
		for line in TE_fh.readlines():
			chrname, chrst, chred, chrleft, strand, TEname, TEClass, TEfamily = line.rstrip().split()[5:13]
			chrst = int(chrst)
			chred = int(chred)
			if chrname in chr_block:
				for block in chr_block[chrname]:
					if (block[0] >= chrst and block[0] <= chred) or (block[1] >= chrst and block[1] <= chred):
						line2TE[block[2]] = ",".join([TEname, TEfamily, TEClass])

	def _annotate_Gene(self):
		pass

	def _is_known(self):
		pass

class HotSpot(object):
	"""docstring for HotSpot"""
	def __init__(self, arg):
		super(HotSpot, self).__init__()
		self.arg = arg
		
	