#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-04-28 11:12:27
# @Author  : Wenliang Wang (wangwl.me@gmail.com)
# @Link    : http://wangwl.me
# @Version : $0.1$

import os
import sys
import re
import pysam
import pandas as pd
import time
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
from subprocess import Popen, PIPE

from epiVIA.integration import Integration
from epiVIA.report import ChimeraSummary, VectorSummary

parser = argparse.ArgumentParser(
    description="Find chimera in alignment result with combined host and vector sequence as reference.")
parser.add_argument("bamfile")
parser.add_argument("outdir", help="A directory where the result files (stat on cells and integrations) written")
parser.add_argument("--Provirus", help="Fasta file of provirus sequence that is integrated in the host genome")
parser.add_argument("--Vector", help="Fasta file of vector sequence")
parser.add_argument("--ltr5_start", type=int, default=1, help="Start position of LTR5, required if using vector")
parser.add_argument("--ltr5_end", type=int, required=True, help="End position of LTR5")
parser.add_argument("--ltr3_start", type=int, help="Start position of LTR3, required if using vector")
parser.add_argument("--ltr3_end", type=int, help="End position of LTR3, required if using vector")
parser.add_argument("--HostIndex", required=True, help="The Index prefix of Host genome")
parser.add_argument("--LTRClipLen", default=11, type=int, help="The soft clipped read length to search in LTR")
parser.add_argument("--HostClipLen", default=17, type=int, help="The soft clipped length to search in host genome")
parser.add_argument("--Host2bit", required=True, help="The 2bit (UCSC) file of host genome sequence")
parser.add_argument("--tempdir", help="A directory for temporary files generated in the pipeline")
parser.add_argument("--quality", type=int, default=0, help="Set the MAPQ for the host read that locate the integration site")
parser.add_argument("--candidate_bam", help="A bamfile to write the suspicious reads used for identifying chimeras, saving time if rerun the pipeline")
parser.add_argument("--chimera_bam", help="A bamfile to write the records that are classified as chimeras")
parser.add_argument("--gbdb", default="http://hgdownload.soe.ucsc.edu/gbdb/", help="gbdb directory of UCSC, default: http://hgdownload.soe.ucsc.edu/gbdb/.")
parser.add_argument("--genome", default='hg38', help="Genome name in UCSC, used for annotation of the integration sites")
args = parser.parse_args()

## Illumina sequence adapters

ME_S7 = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
ME_S7_revc = str(Seq(ME_S7).reverse_complement())
S5_ME = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
S5_ME_revc = str(Seq(S5_ME).reverse_complement())

######################################################   Preflight   ##################################################
if args.Vector:
    try:
        ltr5_start = args.ltr5_start
        ltr5_end = args.ltr5_end
        ltr3_start = args.ltr3_start
        ltr3_end = args.ltr3_end
        ltr_len = ltr5_end - ltr5_start + 1
    except:
        raise Exception("--ltr5_start, --ltr5_end, --ltr3_start, --ltr3_end are required if using Vector sequence.")

    if args.Provirus:
        raise Exception("Only provide provirus or vector sequence")

    try:
        for vector in SeqIO.parse(args.Vector, 'fasta'):
            vector_name = vector.name
            vector_seq = vector.seq
            vector_len = len(vector.seq)
            LTR_SEQ = vector_seq[ltr3_start-1 : ltr3_end]
    except OSError:
        print("Faile to parse the vector fasta file", args.Provirus)

if args.Provirus:
    try:
        args.ltr5_end
    except:
        raise Exception("At least ltr5_end should be provied if using provirus")

    try:
        for vector in SeqIO.parse(args.Provirus, 'fasta'):
            vector_name = vector.name
            vector_seq = vector.seq
            vector_len = len(vector.seq)
            
    except OSError:
        print("Faile to parse the vector fasta file", args.Provirus)
    ltr5_start = 1
    ltr5_end = args.ltr5_end
    ltr_len = ltr5_end
    ltr3_end = vector_len
    ltr3_start = ltr3_end - ltr_len +1
    LTR_SEQ = vector_seq[ltr3_start-1 : ltr3_end]

if not args.Vector and not args.Provirus:
    raise Exception("Vector or Provirus sequence is required")

if args.candidate_bam:
    bamfile_fh = pysam.AlignmentFile(args.bamfile, 'rb')
    candidate_fh = pysam.AlignmentFile(args.candidate_bam, 'wb', template=bamfile_fh)

if args.chimera_bam:
    bamfile_fh = pysam.AlignmentFile(args.bamfile, 'rb')
    chimera_fh = pysam.AlignmentFile(args.chimera_bam, 'wb', template=bamfile_fh)

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
if not os.path.exists(args.tempdir):
    os.makedirs(args.tempdir)

####################################################   Enjoy your trip   ##########################################################

# Chimeric fragments: (1) Chimeric Pair; (2) Host Chimeric read; (3) Vector Chimeric read
class ChimericFragment(object):
    """docstring for Chimera"""

    def __init__(self):
        super(ChimericFragment, self).__init__()
        self.integration = Integration()
        self.fragment_length = 0
        self.reference_name = None
        self.reference_start = None
        self.reference_end = None
        self.fragment_seq = []
        self.name = None
        self.barcode = None
        self.chimeric_type = None


class ChimericPair(ChimericFragment):
    """docstring for ChimericPair"""

    def __init__(self, HostRead, VectorRead, integration=None):
        super(ChimericPair, self).__init__()
        # self.HostRead = HostRead
        # self.VectorRead = VectorRead
        self.integration = integration
        self.reference_name = HostRead.reference_name
        self.reference_start = HostRead.reference_start + 1
        self.reference_end = HostRead.reference_end
        self.vector_start = VectorRead.reference_start + 1
        self.vector_end = VectorRead.reference_end
        self.fragment_seq = [HostRead.seq, VectorRead.seq]
        self.readname = HostRead.qname
        self.barcode = get_barcode(HostRead)
        self.chimeric_type = "Pair"

        if self.vector_start >= ltr3_start and self.vector_end <= ltr3_end:
            self.vector_start = self.vector_start - ltr3_start + ltr5_start
            self.vector_end = self.vector_end - ltr3_start + ltr5_start

        if args.chimera_bam:
            chimera_fh.write(HostRead)
            chimera_fh.write(VectorRead)

class ChimericRead(ChimericFragment):
    """docstring for ChimericRead"""

    def __init__(self, read, integration=None):
        super(ChimericRead, self).__init__()
        # self.read = read
        self.reference_name = read.reference_name
        self.reference_start = read.reference_start + 1
        self.reference_end = read.reference_end
        self.vector_start = None,
        self.vector_end = None,
        self.fragment_seq = [read.seq]
        self.integration = integration
        self.readname = read.qname
        self.barcode = get_barcode(read)
        if read.reference_name == vector_name:
            self.chimeric_type = "VectorRead"
        else:
            self.chimeric_type = "HostRead"

        if args.chimera_bam:
            chimera_fh.write(read)


class VectorFragment(object):
    """docstring for VectorFragment"""

    def __init__(self, first, second):
        super(VectorFragment, self).__init__()
        self.reference_start = first.reference_start + 1
        self.reference_end = second.reference_end
        self.alt_start = None
        self.alt_end = None
        self.readname = first.qname
        self.barcode = get_barcode(first)
        self.pair_cigars = [first.cigarstring, second.cigarstring]
        self.first_alts = getAltanatives(first.tags)
        self.second_alts = getAltanatives(second.tags)
        self.first_align_len = first.query_alignment_length
        self.second_align_len = second.query_alignment_length
        self.is_LTR = False

        if self.reference_start >= ltr5_start and self.reference_end <= ltr5_end:
            self.is_LTR = True
            self.alt_start = ltr3_start + self.reference_start - ltr5_start
            self.alt_end = ltr3_start + self.reference_end - ltr5_start
        elif self.reference_start >= ltr3_start and self.reference_end <= ltr3_end:
            self.is_LTR = True
            self.alt_start = self.reference_start
            self.alt_end = self.reference_end
            self.reference_start = self.reference_start - ltr3_start + ltr5_start
            self.reference_end = self.reference_end - ltr3_start + ltr5_start

    def is_host_alts(self):
        first_is_alt = False
        for alt in self.first_alts:
            if alt[0] != vector_name:
                first_is_alt = True
        second_is_alt = False
        for alt in self.second_alts:
            if alt[0] != vector_name:
                second_is_alt = True
        return first_is_alt and second_is_alt


def get_barcode(read):
    # compatible with 10x bam and snaptools bam

    tags_dict = dict(read.tags)
    first_block = read.qname.split(":")[0]
    if 'CB' in tags_dict:
        barcode = tags_dict['CB']
    elif 'CR' in tags_dict:
        barcode = tags_dict['CR']
    elif re.match("^[ATCGN]+$", first_block):
        barcode = first_block
    else:
        barcode = "None"
    return barcode


def getAltanatives(tags):
    tags = dict(tags)
    alt_aligns = []
    if 'XA' in tags:
        alts = tags['XA'].rstrip(";").split(";")
        for alt in alts:
            alt_name, alt_pos, alt_cigar = alt.split(",")[0:3]
            alt_aligns.append([alt_name, alt_pos, alt_cigar])
    return alt_aligns

# @pysnooper.snoop()


def GetClipFrag(read, clipcutoff):
    cigar = read.cigar
    clipped_frag = ''
    candidate_frag = []
    if len(cigar) <= 1:
        return candidate_frag

    if cigar[0][0] == 4 and cigar[0][1] >= clipcutoff:
        clipped_len = cigar[0][1]
        clipped_frag = read.seq[0:clipped_len]
        candidate_frag.append(clipped_frag)

    if cigar[-1][0] == 4 and cigar[-1][1] >= clipcutoff:
        clipped_len = cigar[-1][1]
        right_clip = cigar[-1][1] * -1
        clipped_frag = read.seq[right_clip:]
        candidate_frag.append(clipped_frag)
    # If the clipped fragment is from Sequencing primer.
    if not read.is_reverse:
        candidate_frag = [x for x in candidate_frag if x not in ME_S7 and x not in S5_ME]
    else:
        candidate_frag = [x for x in candidate_frag if x not in ME_S7_revc and x not in S5_ME_revc]
    return candidate_frag

def ClipInTarget(read, clipcutoff, target):
    # Only clipped on one end of the read is allowed
    candidate_frag = GetClipFrag(read, clipcutoff)

    if len(candidate_frag) == 0:
        return ([], "", 0, 0)

    target_seq_str = str(target)
    target_revcseq_str = str(target.reverse_complement())
    target_seq_len = len(target_seq_str)

    for clipped_frag in candidate_frag:
        clipped_len = len(clipped_frag)
        three_prime_pos = target_seq_len - clipped_len

        plus_hits = [x.start() + 1 for x in re.finditer(clipped_frag, target_seq_str)]
        minus_hits = [x.start() for x in re.finditer(clipped_frag, target_revcseq_str)]
        hits = plus_hits + minus_hits

        if len(plus_hits) == 1 and len(minus_hits) == 0:
            return (hits, "+", clipped_len, three_prime_pos)
        elif len(minus_hits) == 1 and len(plus_hits) == 0:
            return (hits, "-", clipped_len, three_prime_pos)
        else:
            return ([], "", 0, 0)

def FindChimericPair(unpaired_reads, chimera_fragments):
    for readname in unpaired_reads:
        try:
            read1, read2 = unpaired_reads[readname]
        except:
            continue

        if args.candidate_bam:
            candidate_fh.write(read1)
            candidate_fh.write(read2)
        if read1.reference_name == vector_name:
            VectorRead = read1
            HostRead = read2
        elif read2.reference_name == vector_name:
            VectorRead = read2
            HostRead = read1
        else:
            raise Exception("FindChimericPair is processing non-chimeric read pairs.")

        VectorRead_alts = getAltanatives(VectorRead.tags)
        HostRead_alts = getAltanatives(HostRead.tags)
        VectorRead_alts_chrs = [x[0] for x in VectorRead_alts]
        HostRead_alts_chrs = [x[0] for x in HostRead_alts]

        Quality = min(VectorRead.mapq, HostRead.mapq)
        # if HostRead.mapq == 0:
        #     continue
        if vector_name in HostRead_alts_chrs or HostRead.reference_name in VectorRead_alts_chrs:  # if there is host read alternatives on provirus sequence
            continue
        elif len([x for x in VectorRead_alts_chrs if x != vector_name]) > 0:  # if there is vector read alternatives on host genome
            continue
        elif len(VectorRead_alts_chrs) == 1 and vector_name in VectorRead_alts_chrs:
            Quality = HostRead.mapq

        if HostRead.is_reverse != VectorRead.is_reverse:
            InsertOri = "+"
        else:
            InsertOri = "-"
        Chrom = HostRead.reference_name
        ChrStart = HostRead.reference_start + 1
        ChrEnd = HostRead.reference_end
        VectorStart = VectorRead.reference_start + 1
        VectorEnd = VectorRead.reference_end
        '''
        If one of the Read is clipped and the clipped fragment can be found in LTR
        or local sequence within 200 bp of the host read aligned to, we can identify 
        the exact integration loci.
        '''
        ltr_seq = LTR_SEQ if LTR_SEQ else vector_seq[:ltr_len]
        ltr_hits, ltr_strand, clipped_len, ltr_three_prime_pos = ClipInTarget(HostRead, args.LTRClipLen, ltr_seq)

        # Retrieve the local sequence upper or downstream of Host Read aligned position.
        local_seq_st = False
        local_seq_ed = False
        if InsertOri == "+" and VectorRead.cigar[0][0] == 4 and VectorRead.cigar[0][1] >= args.LTRClipLen:
            local_seq_st = HostRead.reference_end
            local_seq_ed = local_seq_st + 200
            host_local_ltr = "LTR5"
        elif InsertOri == "+" and VectorRead.cigar[-1][0] == 4 and VectorRead.cigar[-1][1] >= args.LTRClipLen:
            local_seq_st = HostRead.reference_start - 200
            local_seq_ed = local_seq_st + 200
            host_local_ltr = "LTR3"
        elif InsertOri == "-" and VectorRead.cigar[0][0] == 4 and VectorRead.cigar[0][1] >= args.LTRClipLen:
            local_seq_st = HostRead.reference_start - 200
            local_seq_ed = local_seq_st + 200
            host_local_ltr = "LTR5"
        elif InsertOri == "-" and VectorRead.cigar[-1][0] == 4 and VectorRead.cigar[-1][1] >= args.LTRClipLen:
            local_seq_st = HostRead.reference_end
            local_seq_ed = local_seq_st + 200
            host_local_ltr = "LTR3"

        if local_seq_st and local_seq_ed:
            twoBitToFa_P = Popen("twoBitToFa -seq={chr} -start={start} -end={end} {TwoBit} stdout".format(
                chr=Chrom, start=local_seq_st, end=local_seq_ed, TwoBit=args.Host2bit), shell=True, stdout=PIPE)
            twoBitToFa_P.wait()
            
            try:
                host_local_seq = SeqIO.parse(twoBitToFa_P.stdout, 'fasta').next().seq
                host_local_hits, host_local_strand, clipped_len, host_three_prime_pos = ClipInTarget(VectorRead, args.LTRClipLen, host_local_seq)
            except:
                host_local_hits = []
                host_local_strand = ""
        else:
            host_local_hits = []
            host_local_strand = ""

        LTR = "None"
        if len(ltr_hits) == 1 and ltr_strand == InsertOri and len(host_local_hits) == 0:
            if ltr_hits[0] <= 2:
                LTR = "LTR5"
                ChrStart = HostRead.reference_end
                VectorStart = ltr5_start
                VectorEnd = ltr5_start
            elif ltr_hits[0] >= ltr_three_prime_pos - 2:
                LTR = "LTR3"
                ChrEnd = HostRead.reference_start + 1
                VectorStart = ltr3_end
                VectorEnd = ltr3_end
        elif len(host_local_hits) == 1 and host_local_strand == InsertOri and len(ltr_hits) == 0:
            if host_local_strand == "+":
                ChrStart = local_seq_st + host_local_hits[0] + 1
                ChrEnd = ChrStart
            else:
                ChrStart = 200 - host_local_hits[0] + local_seq_st
                ChrEnd = ChrStart
            if VectorRead.reference_start == 0:
                VectorStart = ltr5_start
                VectorEnd = ltr5_start
                LTR = "LTR5"
            elif VectorRead.reference_end == vector_len:
                VectorStart = ltr3_end
                VectorEnd = ltr3_end
                LTR = "LTR3"
            else:
                ChrStart = HostRead.reference_start + 1
                ChrEnd = HostRead.reference_end

        integration = Integration(Genome=args.genome, Chrom=Chrom, ChrStart=ChrStart, ChrEnd=ChrEnd,
                                  VectorStart=VectorStart, VectorEnd=VectorEnd, LTR=LTR, InsertOri=InsertOri, Quality=Quality)
        chimera = ChimericPair(HostRead=HostRead, VectorRead=VectorRead, integration=integration)
        chimera_fragments.append(chimera)

def FindHostClip(read, clipcutoff, chimera_fragments):
    ltr_seq = LTR_SEQ if LTR_SEQ else vector_seq[:ltr_len]
    hits, strand, clipped_len, three_prime_pos = ClipInTarget(read, clipcutoff, ltr_seq)
    if len(hits) != 1:
        return
    '''
    1. Insertion Orientation.
        Since the read seq in bam file is always in the same orientation with reference sequence, clipped frag strand 
        alone can determine insertion orientation.
        "+" if the clipped fragment is on the positvie strand, "-" if it's on negative
    2. 5' LTR or 3' LTR.
        a. 5' if the integration is at the first base of the sequence, 3' if it's the end
        b. integration site on host should be the end of the aligned position if it's 5' LTR, start position otherwise.
    3. why 3 or three_prime_pos-3
        A few bases (we set 3 at most here) from LTR at the integration site might be able to aligned to host genome 
        together with the host fragment, in this case, we will lost a few bases at the start/end of LTR.
    '''
    if args.candidate_bam:
        candidate_fh.write(read)
    InsertOri = strand

    if hits[0] <= 3:
        LTR = "LTR5"
        ChrStart = read.reference_end
        VectorStart = ltr5_start
        VectorEnd = hits[0] + clipped_len -1
        vector_integration_pos = VectorStart
    elif hits[0] >= three_prime_pos - 3:
        LTR = "LTR3"
        ChrStart = read.reference_start + 1
        VectorStart = hits[0]
        VectorEnd = ltr3_end
        vector_integration_pos = VectorEnd
    else:
        return
    Quality = read.mapq
    ChrEnd = ChrStart
    # VectorEnd = VectorStart
    ChrEnd = ChrStart
    integration = Integration(Genome=args.genome, Chrom=read.reference_name, ChrStart=ChrStart,
                              ChrEnd=ChrEnd, VectorStart=vector_integration_pos, VectorEnd=vector_integration_pos, LTR=LTR, InsertOri=InsertOri, Quality=Quality)
    chimera = ChimericRead(read=read, integration=integration)
    chimera.vector_start = VectorStart
    chimera.vector_end = VectorEnd
    chimera_fragments.append(chimera)

def ReAlignment(fafile, index, chimera_fragments, vector_reads):
    import subprocess
    samfile = fafile + ".sam"
    command = "bwa mem -T {quality} -k {seed} -a -Y  -q {index} {fa} -o {sam}".format(
        index=index, fa=fafile, sam=samfile, quality=args.HostClipLen, seed=args.HostClipLen - 2)
    child = subprocess.Popen(command, shell=True)
    child.wait()
    if child.poll() != 0:
        raise Exception

    re_align = pysam.AlignmentFile(samfile, 'r')
    for clipread in re_align:
        if clipread.mapq == 0:
            continue
        clipread_name = clipread.qname.split(":")
        barcode = clipread_name[0]
        readname = ":".join(clipread_name[1:])
        try:
            read1, read2 = vector_reads.pop(readname)
        except:
            print("{} has more than one non-zero hits in clipped fragment ReAlignment.".format(readname))
            continue
        
        if clipread.is_reverse:
            clipread.seq = str(Seq(clipread.seq).reverse_complement())

        if clipread.seq in read1.seq:
            read = read1
        elif clipread.seq in read2.seq:
            read = read2
        else:
            print(read2.tostring())
            print(read1.tostring())
            sys.exit()
        clipped_loc = read.seq.index(clipread.seq)
        InsertOri = "+" if read.is_reverse == clipread.is_reverse else "-"
        if clipped_loc <= 3:
            VectorStart = read.reference_start + 1
            ChrStart = clipread.reference_end
        else:
            VectorStart = read.reference_end
            ChrStart = clipread.reference_start + 1

        LTR = "None"
        if (InsertOri == "+" and abs(VectorStart - ltr5_start) <= 3) or (InsertOri == "-" and abs(VectorStart - ltr3_end) <= 3):
            LTR = "LTR5"
        elif (InsertOri == "+" and abs(VectorStart - ltr3_end) <= 3) or (InsertOri == "-" and abs(VectorStart - ltr5_start) <= 3):
            LTR = "LTR3"
        else:
            continue
        Chrom = clipread.reference_name
        ChrEnd = ChrStart
        VectorEnd = VectorStart
        integration = Integration(Genome=args.genome, Chrom=Chrom, ChrStart=ChrStart, ChrEnd=ChrEnd, VectorStart=VectorStart,
                                  VectorEnd=VectorEnd, LTR=LTR, InsertOri=InsertOri, CellNumber=1, CellBarcodes=[barcode],
                                  ReadInCell=[1], ReadNumber=1, Quality=clipread.mapq)
        chimera = ChimericRead(read=read, integration=integration)
        chimera.reference_name = clipread.reference_name
        chimera.reference_start = clipread.reference_start +1
        chimera.reference_end = clipread.reference_end
        chimera.vector_start = read1.reference_start +1
        chimera.vector_end = read2.reference_end
        chimera_fragments.append(chimera)

def CorrectLTR(read1, read2):
    read1_alts = getAltanatives(read1.tags)
    read2_alts = getAltanatives(read2.tags)
    insert_size = abs(read1.template_length)

    read1_vector_alts = [x for x in read1_alts if x[0] == vector_name]
    read2_vector_alts = [x for x in read2_alts if x[0] == vector_name]

    # Here we only allow clipped alignment on one read, and full Alignment on the other. 3 bases bias allowed at the clipped site.
    if len(read1_vector_alts) == 1 and len(read2_vector_alts) == 1 and insert_size > read1.qlen:
        if read2.cigar[-1][0] == 4 and "S" not in read1.cigarstring and abs(read2.reference_end - ltr5_end) <= 3:
            read2.reference_start = int(read2_vector_alts[0][1].lstrip("[-+]")) - 1
            read1.reference_start = int(read1_vector_alts[0][1].lstrip("[-+]")) - 1
        elif read1.cigar[0][0] == 4 and "S" not in read2.cigarstring and abs(read1.reference_start - ltr3_start) <= 3:
            read1.reference_start = ltr5_start
            read2.reference_start = int(read2_vector_alts[0][1].lstrip("[-+]")) - 1

        read1.next_reference_start = read2.reference_start
        read2.next_reference_start = read1.reference_start

    return read1, read2


def ProcessVectorAlignments(vector_reads, fa_fh, vector_fragments, clipcutoff):
    for readname in vector_reads:
        try:
            read1, read2 = vector_reads[readname]
        except:
            continue
        '''
        Since there are two identical LTRs on each end of the sequence, reads that clipped at beginning of LTR5 
        and end of LTR3 can also be alternatively aligned to the other. If the read clipped aligned to the end 
        of LTR5 or start of LTR3, and both reads can be alternatively aligned to the other LTR, then the altertive 
        alignment will be chosen as primary.
        '''
        if read1.is_unmapped or read2.is_unmapped:
            continue
        if read1.reference_start > read2.reference_start:
            read1, read2 = read2, read1

        if args.candidate_bam:
            candidate_fh.write(read1)
            candidate_fh.write(read2)
        read1, read2 = CorrectLTR(read1, read2)
        # Only clipped on one end of the read is allowed
        read1_candidate_frag = GetClipFrag(read1, clipcutoff)
        read2_candidate_frag = GetClipFrag(read2, clipcutoff)

        barcode = get_barcode(read1)
        if len(read1_candidate_frag) == 1 and len(read2_candidate_frag) == 0:
            read_fa = ">" + barcode + ":" + read1.qname + "\n" + read1_candidate_frag[0] + "\n"
            fa_fh.write(read_fa)
        elif len(read1_candidate_frag) == 0 and len(read2_candidate_frag) == 1:
            read_fa = ">" + barcode + ":" + read2.qname + "\n" + read2_candidate_frag[0] + "\n"
            fa_fh.write(read_fa)

        vector_pair = VectorFragment(first=read1, second=read2)
        vector_fragments.append(vector_pair)
    return


def ParseBam(bamfile, unpaired_reads, chimera_fragments, vector_reads):
    align = pysam.AlignmentFile(bamfile, 'rb', threads=20)
    for read in align:
        readname = read.qname
        if read.flag & 2 and read.reference_name != vector_name and "S" in read.cigarstring:
            FindHostClip(read, args.LTRClipLen, chimera_fragments)
        elif read.reference_name == vector_name and read.next_reference_name == vector_name:
            vector_reads[readname].append(read)
        elif (not read.flag & 14) and (read.reference_name == vector_name or read.next_reference_name == vector_name):
            unpaired_reads[readname].append(read)
    return


def main():
    chimera_fragments = []
    vector_fragments = []
    vector_reads = defaultdict(list)
    unpaired_reads = defaultdict(list)

    t1 = time.time()
    print("Start Parsing Bam", file=sys.stderr)
    ParseBam(args.bamfile, unpaired_reads, chimera_fragments, vector_reads)
    print("Finish Parsing Bam {}".format(time.time() - t1), file=sys.stderr)
    FindChimericPair(unpaired_reads, chimera_fragments)
    print("Finish ChimericPair {}".format(time.time() - t1), file=sys.stderr)
    fapath = os.path.join(args.tempdir, "Clipped_fragment.fa")
    fa_fh = open(fapath, 'w')
    ProcessVectorAlignments(vector_reads, fa_fh, vector_fragments, args.HostClipLen)
    fa_fh.close()
    if os.path.getsize(fapath) != 0:
        ReAlignment(fapath, args.HostIndex, chimera_fragments, vector_reads)
    print("Finish VectorFragment {}".format(time.time() - t1), file=sys.stderr)

    vector_coverage = os.path.join(args.outdir, "Vector_coverage.csv")
    report_file = os.path.join(args.outdir, "Cell_summary.txt")
    
    chimera_counter = ChimeraSummary(chimera_fragments, args.outdir, args.gbdb, args.quality)
    vector_counter = VectorSummary(vector_fragments, vector_coverage, args.outdir)
    summary_counter = defaultdict(dict)
    for barcode in chimera_counter:
        summary_counter[barcode] = chimera_counter[barcode]

    for barcode in vector_counter:
        if barcode in summary_counter:
            summary_counter[barcode].update(vector_counter[barcode])
        else:
            summary_counter[barcode] = vector_counter[barcode]

    summary_counter_df = pd.DataFrame.from_dict(summary_counter, orient = 'index')
    summary_counter_df.index.name = 'CellBarcode'
    summary_counter_df.to_csv(report_file, sep="\t", na_rep=0)
    print("Finish pipeline {}".format(time.time() - t1), file=sys.stderr)

if __name__ == '__main__':
    main()
