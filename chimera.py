#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-04-28 11:12:27
# @Author  : Wenliang Wang (wangwl.me@gmail.com)
# @Link    : http://wangwl.me
# @Version : $0.1$

import os
import pysam
import argparse
from Bio.Seq import Seq
from Bio import SeqIO, SeqRecord
from collections import defaultdict, Counter
from multiprocessing import Manager, Pool, Process

from integration import IntegrationSite, HotSpot
#from report import Statistics, ChimeraSummary, VectorSummary

parser = argparse.ArgumentParser(description="Find chimera in alignment result with combined host and vector sequence as reference")
parser.add_argument("bamfile")
parser.add_argument("outdir", help="A directory where the result files (stat on cells and integrations) written")
parser.add_argument("--ProvirousSeq", help="Fasta file of provirous sequence that is integrated in the host genome")
parser.add_argument("--LTRLength", type=int)
parser.add_argument("--CombinedIndex", help="The combined Index prefix used in alignment")
parser.add_argument("--HostIndex", help="The Index prefix of Host genome")
parser.add_argument("--KnownSites", help="The known integration sites from previous studies, in the format of RID database")
parser.add_argument("--LTRClipLen", default=11, type = int, help="The soft clipped read length to search in LTR")
parser.add_argument("--HostClipLen", default=17, type = int, help="The soft clipped length to search in host genome")
parser.add_argument("--HostSeq", help="The Host fasta sequence")
parser.add_argument("--LTRseq", help="A fasta file containing the sequence of LTR after integration")
parser.add_argument("--tempdir", help="A directory for temporary files generated in the pipeline")

args = parser.parse_args()

ME_S7 = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
ME_S7_revc = str(Seq(ME_S7).reverse_complement())
S5_ME = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
S5_ME_revc = str(Seq(S5_ME).reverse_complement())
ltr_len = args.LTRLength

try:
	vector=SeqIO.parse(args.ProvirousSeq, 'fasta').next()
	vector_name = vector.name
	vector_seq = vector.seq
	vector_len = len(vector.seq)
except OSError:
	print "Faile to parse the vector fasta file", args.ProvirousSeq
# ChimericPair Attibutes:  Chr ChrStart ChrEnd VectorStart VectorEnd CellBarcode ReadID ClipHost ClipVector ClipedSeq  ClipMatched Read1 Read2

class ChimericFragment(object):
	"""docstring for Chimera"""
	def __init__(self):
		super(ChimericFragment, self).__init__()
		self.integration=None
		self.fragment_length = 0
		self.reference_name = []
		self.reference_start = []
		self.reference_end = []
		self.fragment_seq = []
		self.name = None
		self.barcode = None

class ChimericPair(ChimericFragment):
	"""docstring for ChimericPair"""
	def __init__(self, HostRead, VectorRead, integration):
		super(ChimericPair, self).__init__()
		# self.HostRead = HostRead
		# self.VectorRead = VectorRead
		self.integration = integration
		self.reference_name = [HostRead.reference_name, VectorRead.reference_name]
		self.reference_start = [HostRead.reference_start, VectorRead.reference_start]
		self.reference_end = [HostRead.reference_end, VectorRead.reference_end]
		self.fragment_seq = [HostRead.seq, VectorRead.seq]
		self.readname = HostRead.qname
		self.barcode = HostRead.qname.split(":")[0]
		
class ChimericRead(ChimericFragment):
	"""docstring for ChimericRead"""
	def __init__(self, read, integration):
		super(ChimericRead, self).__init__()
		# self.read = read
		self.reference_name = mate_map_to
		self.integration = integration

class VectorFragment(object):
	"""docstring for VectorFragment"""
	def __init__(self, arg):
		super(VectorFragment, self).__init__()
		self.first = first
		self.second = second

def getAltanatives(tags):
	tags = dict(tags)
	alt_aligns = []
	if 'XA' in tags:
		alts = tags['XA'].split(";")
		for alt in alts:
			alt_name, alt_pos, alt_cigar = alt.split(",")[0:3]
			alt_aligns.append([alt_name, alt_pos, alt_cigar])
	return alt_aligns

def GetClipFrag(read, clipcutoff):
	cigar = read.cigar
	clipped_frag = ''
	
	if cigar[0][0] == 4 and cigar[0][1] >= clipcutoff:
		clipped_len = cigar[0][1]
		clipped_frag = read.seq[0:cigar[0][1]]
	elif cigar[-1][0] == 4 and cigar[-1][1] >= clipcutoff:
		clipped_len = cigar[-1][1]
		right_clip = cigar[-1][1] * -1
		clipped_frag = read.seq[right_clip:]
	## If the clipped fragment is from Sequencing primer.
	if not read.is_reverse:
		if clipped_frag in ME_S7 or clipped_frag in S5_ME:
			clipped_frag = ''
	else:
		if clipped_frag in ME_S7_revc or clipped_frag in S5_ME_revc:
			clipped_frag = ''
	return clipped_frag

def ClipInTarget(read, clipcutoff, target):
	qname = read.qname
	barcode = qname.split(":")[0]
	cigar = read.cigar

	if args.LTRClipLen:
		clipcutoff=args.LTRClipLen
	## Only clipped on one end of the read is allowed
	clipped_frag = GetClipFrag(read, clipcutoff)
	clipped_len = len(clipped_frag)
	if clipped_len < clipcutoff:
		return

	target_seq_str = str(target.seq)
	target_revcseq_str = str(target.reverse_complement().seq)
	target_seq_len = len(target_seq_str)

	plus_hits = [x.start() + 1 for x in re.finditer(clipped_frag, target_seq_str)]
	minus_hits = [x.start() for x in re.finditer(clipped_frag, target_revcseq_str)]
	hits = plus_hits + minus_hits

	three_prime_pos = target_seq_len - clipped_len
	if len(plus_hits) == 1 and len(minus_hits) == 0:
		return hits, "+"
	elif len(minus_hits) == 1 and len(plus_hits) == 0:
		return hits, "-"
	else:
		raise Exception

def FindChimericPair(unpaired_reads, chimera_fragments, id2seq):
	for readname in unpaired_reads:
		barcode = readname.split(":")[0]
		read1, read2 = unpaired_reads[readname]
		
		if read1.reference_name == vector_name:
			VectorRead = read1
			HostRead = read2
		elif read2.reference_name == vector_name:
			VectorRead = read2
			HostRead = read1
		else:
			raise Exception
		host_alts = dict(HostRead.tags)
		if 'XA' in host_alts:
			continue
		InsertOri = "+" if HostRead.is_reverse != VectorRead.is_reverse else "-"
		Chrom=HostRead.reference_name
		ChrStart=HostRead.reference_start + 1
		ChrEnd=HostRead.reference_end
		VectorStart=VectorRead.reference_start + 1
		VectorEnd=VectorRead.reference_end

		ltr_seq = vector_seq[:ltr_len]
		hits, strand = ClipInTarget(HostRead, args.LTRClipLen, ltr_seq)

		if len(hits) == 1 and InsertOri == strand:
			if hits[0] <= 2:
				LTR = "LTR5"
				ChrStart = HostRead.reference_end
				VectorStart = 1
				VectorEnd = 1
			elif hits[0] >= three_prime_pos-2:
				LTR = "LTR3"
				ChrEnd = read.reference_start + 1
				VectorStart = vector_len
				VectorEnd = vector_len
		integration = Integration(Chrom=Chrom, ChrStart=ChrStart, ChrEnd=ChrEnd, VectorStart=VectorStart, VectorEnd=VectorEnd, LTR=LTR, InsertOri=InsertOri, CellNumber=1, CellBarcodes=[barcode], ReadInCell=[1], ReadNumber=1)
		chimera = ChimericPair(HostRead=HostRead, VectorRead=VectorRead, integration=integration)
		chimera_fragments.append(chimera)

def FindHostClip(read, clipcutoff, chimera_fragments):
	ltr_seq = vector_seq[:ltr_len]
	hits, strand = ClipInTarget(read, clipcutoff, ltr_seq)
	if len(hits) != 1:
		return
	'''
	1. Insertion Orientation
		Since the read seq in bam file is always in the same orientation with reference sequence, clipped frag strand alone can determine insertion orientation.
		"+" if the clipped fragment is on the positvie strand, "-" if it's on negative
	2. 5' LTR or 3' LTR
		a. 5' if the integration is at the first base of the sequence, 3' if it's the end
		b. integration site on host should be the end of the aligned position if it's 5' LTR, start position otherwise.
	3. why 2 or three_prime_pos-2
		A few bases (we set 2 at most here) from LTR at the integration site might be able to aligned to host genome together with the host fragment, in this case, we will lost a few bases at the start/end of LTR.
	'''
	InsertOri = strand
	if hits[0] <= 2:
		LTR = "LTR5"
		ChrStart = read.reference_end
		VectorStart = 0
	elif hits[0] >= three_prime_pos-2:
		LTR = "LTR3"
		ChrStart = read.reference_start +1
		VectorStart = vector_len
	else:
		return
	ChrEnd = ChrStart
	VectorEnd = VectorStart
	ChrEnd = ChrStart
	integration = Integration(Chrom=read.reference_name, ChrStart=ChrStart, ChrEnd=ChrEnd, VectorStart=VectorStart, VectorEnd=VectorEnd, LTR=LTR, InsertOri=InsertOri, CellNumber=1, CellBarcodes=[barcode], ReadInCell=[1], ReadNumber=1)

	if integration:
		chimera = ChimericRead(read=read, integration=integration)
		chimera_fragments.append(chimera)

def ReAlignment(fafile, index, chimera_fragments, vector_reads):
	import subprocess
	samfile = fastqfile + ".sam"
	command = "bwa mem -T {quality} -k {seed} -a -Y  -q {index} {fa} -o {sam}".format(index=index, fa=fafile, sam=samfile, quality=args.HostClipLen, seed=quality-2)
	child = subprocess.Popen(command, shell=True)
	child.wait()
	if child.poll() != 0:
		raise Exception

	re_align = pysam.AlignmentFile(samfile, 'r')
	for clipread in re_align:
		if clipread.mapq == 0:
			continue
		readname = clipread.qname
		barcode = readname.split(":")[0]
		read1, read2 = vector_reads.pop(readname)
		if clipread.is_reverse:
			clipread.seq = reverse(clipread.seq)
		if clipread.seq in read1.seq:
			read = read1
		elif clipread.seq in read2.seq:
			read = read2
		else:
			raise Exception
		clipped_loc = read.seq.index(clipread.seq)
		InsertOri = "+" if read.is_reverse == clipread.is_reverse else "-"
		if clipped_loc == 0:
			VectorStart = read.reference_start +1
			ChrStart = clipread.reference_end
		else:
			VectorStart = read.reference_end
			ChrStart = clipread.reference_start +1
		if (InsertOri == "+" and VectorStart == ltr5_start) or (InsertOri == "-" and VectorStart==ltr3_end):
			LTR = "LTR5"
		elif (InsertOri == "+" and VectorStart == ltr3_end) or (InsertOri == "-" and VectorStart==ltr5_start):
			LTR = "LTR3"
		else:
			raise Exception
		Chrom=clipread.reference_name
		ChrEnd = ChrStart
		integration = Integration(Chrom=Chrom, ChrStart=ChrStart, ChrEnd=ChrEnd, VectorStart=VectorStart, VectorEnd=VectorEnd, LTR=LTR, InsertOri=InsertOri, CellNumber=1, CellBarcodes=[barcode], ReadInCell=[1], ReadNumber=1)
		chimera = ChimericRead(read=read, integration=integration)
		chimera_fragments.append(chimera)


def CorrectVectorAlignment(vector_reads, fa_fh, vector_fragments, clipcutoff):
	for readname in vector_reads:
		barcode = readname.split(":")[0]
		read1, read2 = vector_reads[readname]
		'''
		Since there are two identical LTRs on each end of the sequence, reads that clipped at beginning of LTR3 and end of LTR5 can also be alternatively aligned to the other. If the read clipped aligned to the end of LTR3 or start of LTR5, and both reads can be alternatively aligned to the other LTR, then the altertive alignment will be chosen as primary.
		'''
		read1_alts = getAltanatives(read1.tags)
		read2_alts = getAltanatives(read2.tags)
		insert_size = abs(read1.template_length)

		## Here we only allow clipped alignment on one read, and full Alignment on the other. 2 bases bias allowed at the clipped site.
		if len(read1_alts) == 1 and len(read2_alts) == 1 and insert_size > read1.qlen:
			if read1.cigar[-1][0] == 4 and read2.cigar[0][1] == read2.qlen and abs(read1.reference_end - ltr5_end) <= 2:
				read1.reference_start = int(read1_alts[0][1].lstrip("[-+]")) +1
				read2.reference_start = int(read2_alts[0][1].lstrip("[-+]")) +1
			elif read2.cigar[-1][0] == 4 and read1.cigar[0][1] == read1.qlen and abs(read2.reference_end - ltr5_end) <= 2:
				read2.reference_start = int(read2_alts[0][1].lstrip("[-+]")) +1
				read1.reference_start = int(read1_alts[0][1].lstrip("[-+]")) +1
			elif read1.cigar[0][0] == 4 and read2.cigar[0][1] == read2.qlen and abs(read1.reference_start - ltr3_start) <= 2:
				read1.reference_start = ltr5_start
				read2.reference_start = int(read2_alts[0][1].lstrip("[-+]")) +1
			elif read2.cigar[0][0] == 4 and read1.cigar[0][1] == read1.qlen and abs(read2.reference_start - ltr3_start) <=2:
				read2.reference_start = ltr5_start
				read1.reference_start = int(read1_alts[0][1].lstrip("[-+]"))
			read1.next_reference_start = read2.reference_start
			read2.next_reference_start = read1.reference_start

		## Only clipped on one end of the read is allowed
		read1_clip_frag = GetClipFrag(read1, clipcutoff)
		read2_clip_frag = GetClipFrag(read2, clipcutoff)
		read1_clip_len = len(read1_clip_frag)
		read2_clip_len = len(read2_clip_frag)
		if read1_clip_len >= clipcutoff and read2_clip_len < clipcutoff:
			read_fa = ">" + read1.qname + "\n" + read1_clip_frag
			fa_fh.write(read_fa)
		elif read1_clip_len < clipcutoff and read2_clip_len >= clipcutoff:
			read_fa = ">" + read2.qname + "\n" + read2_clip_frag
			fa_fh.write(read_fa)

		if read1.reference_start > read2.reference_start:
			read1, read2 = read2, read1
		vector_pair = VectorFragment(first=read1, second=read2)
		vector_fragments.append(vector_pair)

def ParseBam(bamfile, unpaired_reads, chimera_fragments, vector_reads, id2seq):
	align = pysam.AlignmentFile(bamfile, 'rb')
	for read in align:
		readname = read.qname
		## filter properly paired and unmapped or SE mapped reads
		if read.flag & 2 and read.reference_name != vector_name:
			FindHostClip(read, chimera_fragments, id2seq)
		elif read.reference_name == vector_name and (read.next_reference_name == vector_name or read.mate_is_unmapped):
			vector_reads[readname].append(read)
		elif not read.flag & 14 and (read.reference_name == vector_name or read.next_reference_name == vector_name):
			unpaired_reads[readname].append(read)
	
def ParseFasta(fastafile, id2seq):
	faiter = SeqIO.parse(fastafile, 'fasta')
	for fa in faiter:
		id2seq[fa.name] = fa.seq

def main(bamfile, vector):
	chimera_fragments = []
	vector_fragments = []
	vector_reads = defalutdict(list)
	unpaired_reads = defalutdict(list)

	id2seq = dict()
	ParseFasta(args.reference, id2seq)
	ParseBam(bamfile, unpaired_reads, chimera_fragments, vector_reads)
	id2seq.clear()
	FindChimericPair(unpaired_reads, chimera_fragments)
	fapath = os.path.join(args.tempdir, "Clipped_fragment.fa")
	fa_fh = open(fapath, 'w')
	CorrectVectorAlignment(vector_reads, fa_fh, vector_fragments, args.HostClipLen)
	fa_fh.close()
	ReAlignment(fapath, args.HostIndex, chimera_fragments, vector_reads)
	
	cellstatfile = os.path.join(args.outdir, "Chimeric_cells.list")
	integrationfile = os.path.join(args.outdir, "Integration_sites.list")
	GenerateOutput(chimera_fragments, cellstatfile, integrationfile)

if __name__ == '__main__':
	main()
