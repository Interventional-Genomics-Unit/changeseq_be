import argparse
import HTSeq
import os
import pyfaidx
import regex
import sys
import numpy as np
from utility import reverseComplement
import pysam
import pandas as pd
import traceback


def bam_to_dict(bam, MAPQ=0):
    output = {}
    for read in bam.fetch(region=None):  # region='chr11:46408000-46409000'):
        qname = read.query_name
        if read.is_paired:
            if read.is_unmapped:
                continue
            if read.mapping_quality < MAPQ:
                continue
            if not qname in output:
                output[qname] = {0: [], 1: []}
            if read.is_read1:
                output[qname][0].append(read)
            else:
                output[qname][1].append(read)
    return output


reads_dict = bam_to_dict(pysam.AlignmentFile(bam, "rb"), mapq_threshold)






reference_genome = "/groups/clinical/projects/clinical_shared_data/hg38/hg38.fa"
bam = "/groups/clinical/projects/Assay_Dev/changeseq_be/CD7_test//coverage/CD7_BE_rep1_identified_matched_sorted.bam"
label = "test"
control = "/groups/clinical/projects/Assay_Dev/changeseq_be/CD7_test/coverage/Control_CD7_BE_rep1_matched_sorted.bam"
targetsite = "CCCTACCTGTCACCAGGACCNGN"
search_radius = 30
window_size = 30
mapq_threshold = 50
mismatch_threshold = 6
output_dir = "/groups/clinical/projects/Assay_Dev/changeseq_be/CD7_test//aligned/"
read_count_cutoff = 6
edited_read_cutoff = 1
