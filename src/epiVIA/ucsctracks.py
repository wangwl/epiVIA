#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-05-15 22:58:29
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import sys
import re
import json
import yaml
import argparse
import requests
from subprocess import Popen, PIPE

from epiVIA.utils import bb_conf, get_abs_path

UCSC_API = "https://api.genome.ucsc.edu/"

class UCSCTrack(object):
	"""docstring for UCSCTrack"""
	def __init__(self, **kwargs):
		super(UCSCTrack, self).__init__()


# @pysnooper.snoop()
def run_bigBedToBed(track, genome, chrom=None, start=None, end=None, columns=None, gbdb="http://hgdownload.soe.ucsc.edu/gbdb/"):
	try:
		file_path = "{gbdb}/{genome}/{bb_path}".format(gbdb=gbdb, genome=genome, bb_path=bb_conf[track])
		cmd = "{bigBedToBed_path} {bb_path} stdout".format(bigBedToBed_path=get_abs_path('bigBedToBed'), bb_path=file_path)
		if chrom:
			cmd += " -chrom={chrom}".format(chrom=chrom)
		if start:
			cmd += " -start={start}".format(start=start)
		if end:
			cmd += " -end={end}".format(end=end)

		child = Popen(cmd, shell=True, stdout=PIPE, universal_newlines=True)
		child.wait()
		# if child.poll() != 0:
		# 	raise Exception
		output = []
		for line in child.stdout.readlines():
			row = str(line).rstrip().split("\t")
			record = []
			if len(columns) > 0:
				for x in columns:
					record.append(row[x])
			else:
				record = row
			output.append(record)
		return output
	except:
		print("fail to annotate {}, continue anyway.".format(track))
		return []
	## bigBedToBed subprocess command line
	
	

#@pysnooper.snoop()
def find_ucsctrack(track, genome, chrom=None, start=None, end=None, keys=[]):
	action = "getData"
	datatype = "track"
	track_name = track

	url = UCSC_API + "{action}/{datatype}?genome={genome};track={track}".format(action = action, datatype=datatype, genome=genome, track=track_name)
	if chrom:
		url += ";chrom={chrom}".format(chrom=chrom)
	if start:
		url += ";start={start}".format(start=start)
	if end:
		url += ";end={end}".format(end=end)
	request = requests.get(url)

	# if request.status_code != 200:
	# 	print url
	if request.status_code == 200 and request.headers['Content-Type'] == 'application/json':
		content = yaml.safe_load(request.text)
		items = content[track_name]
		output = []
		for record in items:
			line = []
			if len(keys) > 0:
				for key in keys:
					line.append(record[key])
				output.append(line)
			else:
				output.append(record)
		return output
	else:
		print("fail to annotate {} with UCSC API, continue anyway.".format(track))
		return []

def find_trackhub():
	pass

def chunk_image(imgdir, **kwargs):
	IMAGE_API = "http://genome.ucsc.edu/cgi-bin/hgRenderTracks?"
	request = requests.get(IMAGE_API, params=kwargs)
	if request.status_code != 200:
		raise Exception

	image_file = os.path.join(imgdir, kwargs['position']+".png")
	image_fh = open(image_file, 'w')
	image_fh.write(request.content)
	image_fh.close()

def main():
	parser = argparse.ArgumentParser(description='annotate chimeric region with TE and gene features')
	parser.add_argument("--bedfile", help="a bedfile with the chimeric region")
	parser.add_argument("track")
	parser.add_argument("--track_keys", nargs="+")
	parser.add_argument("--bedcolumns", nargs="+", type=int)
	parser.add_argument("--sequence")
	parser.add_argument("--hubUrl")
	parser.add_argument("--genome")
	parser.add_argument("--chrom")
	parser.add_argument("--start")
	parser.add_argument("--end")
	parser.add_argument("--imgdir")
	parser.add_argument("--gbdb", default="http://hgdownload.soe.ucsc.edu/gbdb/")
	args = parser.parse_args()
	if args.track_keys:
		items = find_ucsctrack(args.track, args.genome, args.chrom, args.start, args.end,  args.track_keys)
	else:
		items = find_ucsctrack(args.track, args.genome, args.chrom, args.start, args.end, args.bedcolumns)
	for x in items:
		print("\t".join(x))

if __name__ == '__main__':
	main()