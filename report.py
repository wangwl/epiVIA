#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-05-10 14:36:28
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os

class Statistics(object):
	"""docstring for Statistics"""
	def __init__(self, arg):
		super(Statistics, self).__init__()
		self.arg = arg
		

def ChimeraSummary(chimera_fragments, cellstatfile, integrationfile):
	cellCounter = Counter()
	Integration_sites = defalutdict(list)
	try:
		cell_fh = open(cellstatfile, 'w')
	except OSError:
		print "Can't open", cellstatfile
	try:
		integration_fh = open(integrationfile, 'w')
	except OSError:
		print "Can't open", integrationfile
	for chimera in chimera_fragments:
		read = chimera.read
		barcode = read.qname.split(":")[0]
		mark_dup_key = (barcode, read.reference_name, read.reference_start, read.reference_end, read.next_reference_start, read.next_reference_end)
		if mark_dup_key in mark_dup:
			continue
		else:
			mark_dup[mark_dup_key] = 1
		cellCounter[barcode] += 1
		integration = chimera.integration
		integration_key = (integration.chrom, integration.ChrStart)
		Integration_sites[integration_key].append(integration)

	for integration_key in Integration_sites:
		chrom, ChrStart = integration_key
		integration_hubs = Integration_sites[integration_key]
		readsnum = len(integration_hubs)
		site_barcode_counter = Counter()
		for integration in integration_hubs:
			site_barcode_counter[integration.CellBarcodes[0]] += 1

def VectorSummary(vector_fragments):
	
def MappingDepth():
	pass
