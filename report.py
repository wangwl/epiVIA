#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-05-10 14:36:28
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import sys
from collections import Counter, defaultdict

from integration import Integration

class Statistics(object):
	"""docstring for Statistics"""
	def __init__(self, arg):
		super(Statistics, self).__init__()
		self.chimeric_frags = None
		self.chimeric_pairs = None
		self.chimeric_reads = None

		
def ChimeraSummary(chimera_fragments, integrationfile):
	cellCounter = Counter()
	Integration_sites = defaultdict(list)
	# try:
	# 	cell_fh = open(cellstatfile, 'w')
	# except OSError:
	# 	print "Can't open", cellstatfile
	try:
		integration_fh = open(integrationfile, 'w')
	except OSError:
		print "Can't open", integrationfile

	total_chimera_fragments = len(chimera_fragments)
	readCounter = Counter()
	rmdupReadCounter = Counter()
	mark_dup = dict()
	for chimera in chimera_fragments:
		barcode = chimera.barcode
		readCounter[chimera.chimeric_type] += 1
		reference_name, next_reference_name = chimera.reference_name
		reference_start, next_reference_start = chimera.reference_start
		reference_end, next_reference_end = chimera.reference_end
		mark_dup_key = (barcode, reference_name, reference_start, reference_end, next_reference_start, next_reference_end)
		if mark_dup_key in mark_dup:
			continue
		else:
			mark_dup[mark_dup_key] = 1
		rmdupReadCounter[chimera.chimeric_type] += 1
		cellCounter[barcode] += 1
		integration = chimera.integration
		integration.CellBarcodes = barcode
		integration.ChimeraSource = chimera.chimeric_type
		integration_key = (integration.Chrom, integration.ChrStart, integration.ChrEnd)
		Integration_sites[integration_key].append(integration)
	
	for integration_key in Integration_sites:
		Chrom, ChrStart, ChrEnd = integration_key
		integration_hubs = Integration_sites[integration_key]
		reads_num = len(integration_hubs)
		barcode_counter = Counter()
		source_counter = Counter()
		for integration in integration_hubs:
			barcode_counter[integration.CellBarcodes] += 1
			source_counter[integration.ChimeraSource] += 1
		VectorStart = integration_hubs[0].VectorStart
		VectorEnd = integration_hubs[0].VectorEnd
		LTR = integration_hubs[0].LTR
		InsertOri = integration_hubs[0].InsertOri
		Genome = integration_hubs[0].Genome

		CellNumber = len(barcode_counter)
		CellBarcodes_tuples = sorted(barcode_counter.items(), key=lambda x:x[1], reverse=True)
		CellBarcodes = [x[0] for x in CellBarcodes_tuples]
		ReadInCell = [x[1] for x in CellBarcodes_tuples]
		ReadNumber = sum(ReadInCell)
		ChimeraSource_tuples = sorted(source_counter.items(), key=lambda x:x[1], reverse=True)
		ChimeraSource = [x[0] for x in ChimeraSource_tuples]
		ChimeraSource_count = [x[1] for x in ChimeraSource_tuples]

		LTR="None" if LTR==None else LTR
		grouped_integration = Integration(Genome=Genome, Chrom=Chrom, ChrStart=ChrStart, ChrEnd=ChrEnd, VectorStart=VectorStart, VectorEnd=VectorEnd, LTR=LTR, InsertOri=InsertOri, CellBarcodes=CellBarcodes, ReadInCell=ReadInCell, CellNumber=CellNumber, ReadNumber=ReadNumber, ChimeraSource=ChimeraSource)
		grouped_integration.annotate_TE()
		grouped_integration.annotate_Gene()
		grouped_integration.annotate_Enhancer()
		grouped_integration.nearest_gene()

		try:
			header
			integration_string, header = grouped_integration.tostring()
			integration_fh.write(integration_string)
		except:
			integration_string, header = grouped_integration.tostring()
			integration_fh.write(header)
			integration_fh.write(integration_string)
	integration_fh.close()

def VectorSummary(vector_fragments):
	pass
def MappingDepth():
	pass
