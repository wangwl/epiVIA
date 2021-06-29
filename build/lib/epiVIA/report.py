#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-05-10 14:36:28
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import pandas
from collections import Counter, defaultdict

from epiVIA.integration import Integration

class Statistics(object):
	"""docstring for Statistics"""

	def __init__(self, arg):
		super(Statistics, self).__init__()
		self.chimeric_frags = None
		self.chimeric_pairs = None
		self.chimeric_reads = None


def process_integration_site(grouped_integration, gbdb):
	grouped_integration.annotate_TE()
	grouped_integration.annotate_Gene(gbdb)
	# grouped_integration.annotate_Enhancer(gbdb)
	# grouped_integration.nearest_gene()
	return grouped_integration

def ChimeraSummary(chimera_fragments, outdir, gbdb, quality):
	cellCounter = defaultdict(Counter)
	Integration_sites = defaultdict(list)
	Integrated_cells = defaultdict(list)

	try:
		integration_fh = open(os.path.join(outdir, "Integration_sites.list"), 'w')
	except OSError:
		print("Can't open", os.path.join(outdir, "Integration_sites.list"))

	mark_dup = dict()
	for chimera in chimera_fragments:
		barcode = chimera.barcode

		vector_start = chimera.vector_start
		vector_end = chimera.vector_end
		reference_name = chimera.reference_name
		reference_start = chimera.reference_start
		reference_end = chimera.reference_end
		integration = chimera.integration

		mark_dup_key = (barcode, vector_start, vector_end)
		if mark_dup_key in mark_dup:
			if reference_name == mark_dup[mark_dup_key][0] and (mark_dup[mark_dup_key][1] == reference_start or mark_dup[mark_dup_key][2] == reference_end):
				continue
			elif integration.Quality < 30:
				continue
		else:
			mark_dup[mark_dup_key] = (reference_name, reference_start, reference_end)

		cellCounter[barcode]['ChimeraSum'] +=1
		cellCounter[barcode][chimera.chimeric_type] += 1
		
		integration.CellBarcodes = barcode
		integration.ChimeraSource = chimera.chimeric_type
		integration_key = (integration.Chrom, integration.ChrStart, integration.ChrEnd)
		Integration_sites[integration_key].append(integration)
	
	grouped_integration_list = []
	for integration_key in Integration_sites:
		Chrom, ChrStart, ChrEnd = integration_key
		integration_hubs = Integration_sites[integration_key]
		
		barcode_counter = Counter()
		source_counter = Counter()
		Quality = 0
		for integration in integration_hubs:
			barcode_counter[integration.CellBarcodes] += 1
			source_counter[integration.ChimeraSource] += 1
			if Quality < integration.Quality:
				Quality = integration.Quality
		if Quality < quality:
			continue
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
		# ChimeraSource_count = [x[1] for x in ChimeraSource_tuples]
		for barcode in CellBarcodes:
			cellCounter[barcode]['integration'] += 1

		LTR = "None" if LTR == None else LTR
		grouped_integration = Integration(Genome=Genome, Chrom=Chrom, ChrStart=ChrStart, ChrEnd=ChrEnd, VectorStart=VectorStart, VectorEnd=VectorEnd, LTR=LTR, InsertOri=InsertOri, CellBarcodes=CellBarcodes, ReadInCell=ReadInCell, CellNumber=CellNumber, ReadNumber=ReadNumber, ChimeraSource=ChimeraSource, Quality=Quality)
		grouped_integration = process_integration_site(grouped_integration, gbdb)
		grouped_integration_list.append(grouped_integration)

	# p = Pool(n_jobs)
	# for grouped_integration in grouped_integration_list:
	# 	p.apply_async(process_integration_site, (grouped_integration, gbdb))
	# p.close()
	# p.close()

	for grouped_integration in grouped_integration_list:
		if grouped_integration.InsertGene == 'LINC00486':
			continue
		try:
			header
			integration_string, header = grouped_integration.tostring()
			integration_fh.write(integration_string)
		except:
			integration_string, header = grouped_integration.tostring()
			integration_fh.write(header)
			integration_fh.write(integration_string)

		for barcode in grouped_integration.CellBarcodes:
			Integrated_cells[barcode].append(grouped_integration)
	integration_fh.close()

	barcode_integration = defaultdict(dict)
	for barcode in Integrated_cells:
		if len(Integrated_cells[barcode]) == 1:
			integration = Integrated_cells[barcode][0]
			barcode_integration[barcode]['Chrom'] = integration.Chrom
			barcode_integration[barcode]['ChrStart'] = integration.ChrStart
			barcode_integration[barcode]['ChrEnd'] = integration.ChrEnd
			barcode_integration[barcode]['VectorStart'] = integration.VectorStart
			barcode_integration[barcode]['VectorEnd'] = integration.VectorEnd
			barcode_integration[barcode]['InsertGene'] = integration.InsertGene
			barcode_integration[barcode]['ExonIntron'] = integration.ExonIntron
			barcode_integration[barcode]['Enhancer'] = integration.Enhancer
			barcode_integration[barcode]['TEFamily'] = integration.TEFamily
			barcode_integration[barcode]['TEClass'] = integration.TEClass
			barcode_integration[barcode]['Quality'] = integration.Quality
		else:
			i = 0
			for integration in Integrated_cells[barcode]:
				if barcode not in barcode_integration:
					barcode_integration[barcode]['Chrom'] = integration.Chrom
					barcode_integration[barcode]['ChrStart'] = integration.ChrStart
					barcode_integration[barcode]['ChrEnd'] = integration.ChrEnd
					barcode_integration[barcode]['VectorStart'] = integration.VectorStart
					barcode_integration[barcode]['VectorEnd'] = integration.VectorEnd
					barcode_integration[barcode]['InsertGene'] = integration.InsertGene
					barcode_integration[barcode]['ExonIntron'] = integration.ExonIntron
					barcode_integration[barcode]['Enhancer'] = integration.Enhancer
					barcode_integration[barcode]['TEFamily'] = integration.TEFamily
					barcode_integration[barcode]['TEClass'] = integration.TEClass
					barcode_integration[barcode]['Quality'] = integration.Quality			
				elif integration.Chrom == barcode_integration[barcode]['Chrom'] and abs(integration.ChrStart - barcode_integration[barcode]['ChrStart']) <= 1000:
					ChrStart, ChrEnd = sorted([integration.ChrStart, integration.ChrEnd, barcode_integration[barcode]['ChrStart'], barcode_integration[barcode]['ChrEnd']])[1:3]
					barcode_integration[barcode]['ChrStart'] = ChrStart
					barcode_integration[barcode]['ChrEnd'] = ChrEnd
				else:
					i += 1
					barcode = "{}.{}".format(barcode, i)
					barcode_integration[barcode]['Chrom'] = integration.Chrom
					barcode_integration[barcode]['ChrStart'] = integration.ChrStart
					barcode_integration[barcode]['ChrEnd'] = integration.ChrEnd
					barcode_integration[barcode]['VectorStart'] = integration.VectorStart
					barcode_integration[barcode]['VectorEnd'] = integration.VectorEnd
					barcode_integration[barcode]['InsertGene'] = integration.InsertGene
					barcode_integration[barcode]['ExonIntron'] = integration.ExonIntron
					barcode_integration[barcode]['Enhancer'] = integration.Enhancer
					barcode_integration[barcode]['TEFamily'] = integration.TEFamily
					barcode_integration[barcode]['TEClass'] = integration.TEClass
					barcode_integration[barcode]['Quality'] = integration.Quality
	barcode_integration_df = pandas.DataFrame.from_dict(barcode_integration, orient='index', columns=['Chrom', 'ChrStart', 'ChrEnd', 'VectorStart', 'VectorEnd', 'InsertGene', 'ExonIntron', 'Enhancer', 'TEFamily', 'TEClass', 'Quality'])
	barcode_integration_df.index.name = 'CellBarcode'
	barcode_integration_df.to_csv(os.path.join(outdir, "barcode_integration.list"), sep="\t")
	return cellCounter

def VectorSummary(vector_fragments, vector_coverage, outdir):
	cellCounter = defaultdict(Counter)
	Uniq_depth = defaultdict(Counter)
	positions = dict()
	mark_dup = dict()
	for Vectorfrag in vector_fragments:
		barcode = Vectorfrag.barcode
		if Vectorfrag.is_host_alts():
			continue
		mark_dup_key = (Vectorfrag.reference_start, Vectorfrag.reference_end)
		if mark_dup_key in mark_dup:
			continue
		else:
			mark_dup[mark_dup_key] = 1

		cellCounter[barcode]['VectorSum'] += 1

		if Vectorfrag.is_LTR:
			cellCounter[barcode]['LTR'] += 1
			for pos in range(Vectorfrag.reference_start-1, Vectorfrag.reference_end):
				positions[pos] = 1
				Uniq_depth[barcode][pos] += 0.5
			for pos in range(Vectorfrag.alt_start-1, Vectorfrag.alt_end):
				positions[pos] = 1
				Uniq_depth[barcode][pos] += 0.5
		else:
			cellCounter[barcode]['noLTR'] += 1
			for pos in range(Vectorfrag.reference_start-1, Vectorfrag.reference_end):
				positions[pos] = 1
				Uniq_depth[barcode][pos] += 1

	CellBarcodes = cellCounter.keys()
	# csv_path=os.path.join(outdir, "VectorFragment.csv")
	depth_fh = open(vector_coverage, 'w')
	depth_fh.write("barcode,MaxDepth,FragNum,"+",".join([str(x) for x in sorted(positions.keys())]) + "\n")
	for barcode in CellBarcodes:
		max_depth = max(Uniq_depth[barcode].values())
		output = "{},{},{}".format(barcode, max_depth, cellCounter[barcode]['VectorSum'])
		for pos in sorted(positions.keys()):
			output += ",{}".format(Uniq_depth[barcode][pos])
		depth_fh.write(output+"\n")
	return cellCounter