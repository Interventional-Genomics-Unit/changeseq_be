import subprocess
import gzip
import os
import logging
from utility import write_default_yaml

logger = logging.getLogger('root')
logger.propagate = False

def process_refseq(tmp_output, output):
    '''
    ftp : https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz

    Makes
    A) a bed file from refseqs = 282,614 lines
    B) a bed file of the most recent genes = 66688 lines
    '''

    chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
              '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'Y', 'X']

    #in
    ref_out = open(output, 'w')

    labels = ['bin', 'id', 'chrom', 'strand', 'txStart', 'txEnd',
              'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds',
              'score', 'name', 'cdsStartStat', 'cdsEndStat','exonFrames']

    cnt = 0
    for line in gzip.open(tmp_output, 'rt'):
        tokens = line.split('\t')
        tid, chrom, strand,tstart,tend = tokens[1:6]
        cds_start, cds_end = tokens[6], tokens[7]
        exons_start, exon_end = tokens[9], tokens[10]
        gname, frames = tokens[12], tokens[-1].split('\n')[0]

        if tokens[2].replace('chr', "") in chroms:
            tid = tokens[1]
            cnt+=1
            if tid.startswith('X') == False:
                eid = '-'
                line_out = [chrom,tstart,tend,'|'.join([strand,tid,eid,gname,cds_start,cds_end, exons_start,exon_end,frames])]
                ref_out.write('\t'.join(line_out) + '\n')
    ref_out.close()

def get_refseq(ftp_path,tmp_output):
    cmd = "wget " + ftp_path + " -O " +  tmp_output
    subprocess.check_call(cmd, shell=True)


def makefiles(ftp_path,p_dir,reset_output):
    outdir = p_dir + "/data/"
    if not os.path.exists(outdir):
        logger.info("Creating data directory" + outdir)
        os.makedirs(outdir)


    if reset_output:
        logger.info("Reseting changeseqs annotation file to "+ reset_output)
        write_default_yaml('annotate_path', reset_output)
    else:
        logger.info("downloading Refseq from " + ftp_path)
        tmp_output = outdir + "tmp_ncbiRefSeq.txt.gz"
        output = outdir + "ncbiRefSeq.bed"
        get_refseq(ftp_path, tmp_output)

        logger.info("Cleaning Refseq and converting to bed file")
        process_refseq(tmp_output, output)

        logger.info("Storing file path to "+ output)
        write_default_yaml('annotate_path', output)

        logger.info("Cleaning up tmp files")
        os.remove(tmp_output)


