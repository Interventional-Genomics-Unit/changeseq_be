from subprocess import Popen, PIPE
import subprocess
import os
from Bio.Seq import Seq
import pandas as pd


class Transcript:


    tx_lib = {} #tid : []
    coord2tid = {}
    labels = ['chrom', 'txStart', 'txEnd', 'strand', 'tid', 'eid', 'name',
              'cdsStart', 'cdsEnd', 'exonStarts', 'exonEnds', 'exonFrames']


    def __init__(self,entry):
        '''
        :param entry: NCBI Transcript Entry
        '''

        self.entry = entry
        for k in ['txStart', 'txEnd','cdsStart', 'cdsEnd']:
            if k in self.entry.keys():
                self.entry[k] = int(self.entry[k])

        self.overlapping_transcripts = []

        ##output
        #mapping transcript features (relative to transcript start)
        if 'exonStarts' in self.entry.keys():
            self.exons= self.get_exons()
            self.utrs = self.get_utrs()
        self.tx_len = self.entry['txEnd'] - self.entry['txStart']  # 1608

        self.feature = None
        self.rf = None
        self.cdsseq = None
        self.txseq = None

    @classmethod
    def transcript(cls, snvcoord):
        if snvcoord in cls.coord2tid.keys():
            tid = cls.coord2tid[snvcoord]
            start, end = snvcoord.split(':')[1].split('-')
            pos = int((int(start) + int(end)) /2)
            obj = cls.tx_lib[tid]
            obj.find_feature(pos)
            return obj
        else:
            return 'Not found'


    @classmethod
    def load_transcripts(cls, annote_path,snvcoords):
        # print('DANIEL WE MIGHT NEED TO MAKE A TEMP FOLDER')
        # print('SEE load_transcripts in annotate.py')
        #temp_bedfname = "temp_in.bed"
        #temp_bedfname_sorted = "temp_in_sorted.bed"

        bed_data = ""

        #reset dict
        cls.coord2tid = {}
        cls.tx_lib = {}

        #with open(temp_bedfname, 'w') as tempbed:
        for coord in snvcoords:
            coord_field = coord.split(':')
            chrom = coord_field[0]
            start = coord_field[1].split('-')[0]
            end = coord_field[1].split('-')[1]
            line = "\t".join([chrom, start, end])
            bed_data += line + "\n"

        #cmd = 'bedtools window -w 50 -a ' + temp_bedfname + ' -b ' + annote_path

        #Sort bedfile
        sort_cmd = ['bedtools', 'sort', '-i', "stdin"]

        sorted_bed = subprocess.Popen(
            sort_cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True
        )

        sorted_output, sort_error = sorted_bed.communicate(input=bed_data)
        sorted_bed.wait()

        #find closest ORFS
        # Check for errors in the sorting step
        if sort_error:
            print("Error during sorting:", sort_error)
        else:
            # Pipe the sorted BED data directly to `bedtools window`
            window_cmd = ["bedtools", "closest", "-a",  "stdin", "-b", annote_path]
            window_output = subprocess.Popen(
                window_cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            # Run the command and capture the output
            bedtools_out, stderr = window_output.communicate(input=sorted_output)
            window_output.wait()

            if window_output.returncode == 0:
                pass
            else:
                print("Error:", stderr)

        #cmd = 'bedtools closest -a ' + temp_bedfname_sorted + ' -b ' + annote_path
        #bed_popen = Popen(cmd, shell=True, stdout=PIPE)
        #bedtools_out = bed_popen.communicate()[0]
        #et = bed_popen.wait()

        #if ret != 0:
         #   print("Bedtools process was interrupted!")

        if bedtools_out != '':
            bed_entries = bedtools_out.split('\n')[:-1]

            for bed_entry in bed_entries:
                tokens = bed_entry.split('\t')
                tokens = tokens[:-1] + tokens[-1].split('|')
                snvcoord = tokens[0] + ':' + tokens[1] + '-' + tokens[2]
                entry = dict(zip(cls.labels, tokens[-12:]))

                if entry['chrom'] != '.':
                    if snvcoord not in cls.coord2tid.keys():
                        cls.coord2tid[snvcoord] = entry['tid']
                        cls.tx_lib[entry['tid']] = cls(entry)

                    # adjusting if multiple ORFS found to retrieve the most up to date transcripts
                    else:
                        new_tid = entry['tid']
                        oldtid = cls.coord2tid[snvcoord]

                        if new_tid.startswith('NR') and oldtid.startswith('NM'):
                            pass

                        elif new_tid[0:2] == oldtid[0:2]:
                            if entry['eid'] != '-' or (float(oldtid.split('_')[-1]) < float(new_tid.split('_')[-1])):
                                cls.coord2tid[snvcoord] = new_tid
                                obj = cls(entry)
                                obj.overlapping_transcripts.append(oldtid)
                                cls.tx_lib[entry['tid']] = obj
                        else:
                            cls.coord2tid[snvcoord] = new_tid
                            obj = cls(entry)
                            obj.overlapping_transcripts.append(oldtid)
                            cls.tx_lib[entry['tid']] = obj

                else:
                    cls.coord2tid[snvcoord] = 'Error Not Found'


    def __dict__(self):
        return self.entry

    def get_utrs(self):
        exon_starts = self.entry['exonStarts'][:-1].split(',')
        exon_ends = self.entry['exonEnds'][:-1].split(',')
        exon_frames = self.entry['exonFrames'].replace("\n", "")[:-1].split(',')
        utrs = None

        if self.exons:
            ogexons = [(int(exon_starts[i]) - self.entry['txStart'], int(exon_ends[i]) - self.entry['txStart']) for i in range(len(exon_ends))]

            utrs = [(ogexons[0][0], self.exons[0][0] - 1),(self.exons[-1][1], ogexons[-1][1] - 1)]

            if exon_frames[-1] == '-1':  # -1 means entire exon is UTR
                utrs[1] = ogexons[-1]
            if exon_frames[0] == '-1':
                utrs[0] = ogexons[0]

        return utrs

    def get_exons(self):
        '''
        uses entry info to find relative exons positions
        '''
        exon_starts = self.entry['exonStarts'][:-1].split(',')
        exon_ends = self.entry['exonEnds'][:-1].split(',')
        exon_frames = self.entry['exonFrames'].replace("\n", "")[:-1].split(',')
        tx_start = self.entry['txStart']

        exons = [(int(exon_starts[i]) - tx_start, int(exon_ends[i]) - tx_start) for i in range(len(exon_ends))]
        for i in range(len(exon_frames)):
            if exon_frames[i] == '-1':  # -1 means entire exon is UTR
                exons = exons[1:]
                exon_starts = exon_starts[1:]
            else:
                break
        for i in range(1, len(exon_frames)):
            if exon_frames[-i] == '-1':
                exons = exons[0:len(exons) - 1]
                exon_ends = exon_ends[0:-1]
            else:
                break

        # Determine the stop and start of UTR
        if len(exons) > 0:
            exons[0] = (int(self.entry['cdsStart']) - int(exon_starts[0]) + exons[0][0], exons[0][1])
            exons[-1] = (exons[-1][0], exons[-1][1] - (int(exon_ends[-1]) - int(self.entry['cdsEnd'])))
        else:
            exons = None
        return exons


    def get_cdsseq(self):
        '''
        uses entry info to adjust exons positions relative to transcription start
        removes utrs from exons
        translates into cds
        '''
        # Determine the stop and start of UTR
        if self.tx.seq:
            cds = Seq(''.join([str(self.tx_seq)[a:b] for a, b in self.exons]))
            if self.entry['strand'] == '-':
                cds = cds.reverse_complement()
        else:
            cds = None
        return cds


    def tx_info(self):
        return self.entry['eid'], self.entry['tid'], self.entry['name'], self.entry['strand'], self.entry['txStart']


    def find_reading_frame(self, dist_from_cds_start):
        '''
        Finds reading frame of SNV in extracted sequence
        '''
        rf = 1 if dist_from_cds_start % 3 == 2 else 2 if dist_from_cds_start % 3 == 0 else 0
        return rf

    def find_feature(self, pos):
        feature, rf = None, None
        t_snvpos = int(pos) - self.entry['txStart']
        #cdstart, cdsend = self.cds_start, self.cds_end
        if pos not in range(int(self.entry['txStart']),int(self.entry['txEnd'])):
            dist = abs(pos - int(self.entry['txEnd'])) if abs(pos - int(self.entry['txStart'])) > abs(pos - int(self.entry['txEnd'])) else abs(pos - int(self.entry['txStart']))
            feature = 'intergenic | ' + str(dist) + 'bp from ' + self.entry['name']


        elif self.entry['tid'].startswith('NR'):
            feature = 'ncRNA'

        elif 'cdsStart' not in self.entry.keys() :
            feature = 'unknown'

        elif self.exons == None or t_snvpos < -50 or t_snvpos > self.tx_len + 50:
            # not in transcript - shouldn't happen or else no entry would be found
            feature = 'non-coding'


        elif t_snvpos in range(self.utrs[0][0],self.utrs[0][1]+1) or t_snvpos in range(self.utrs[1][0],self.utrs[1][1]+1):

            if t_snvpos in range(self.utrs[0][0],self.utrs[0][1]+1):
                feature = '5utr' if self.entry['strand'] == '+' else '3utr'
            else:
                feature = '3utr' if self.entry['strand'] == '+' else '5utr'

        elif t_snvpos in range(self.exons[0][0],self.exons[0][0]+4) or t_snvpos in range(self.exons[-1][-1]-3, self.exons[-1][-1]+1):

            if t_snvpos in range(self.exons[0][0],self.exons[0][0]+4):
                feature = "start_codon" if self.entry['strand'] == '+' else 'stop_codon'
            else:

                feature = 'stop_codon' if self.entry['strand'] == '+' else 'start_codon'

        else:# find if exon or intron
            exon_n = 0
            feature = 'intron'
            for x in self.exons:
                # if in exon find reading frame
                if abs(t_snvpos - x[0]) <=3 or abs(t_snvpos - x[1]) <=3:
                    feature = 'splice site'
                    break
                elif t_snvpos in range(x[0], x[1] + 1):
                    feature = 'exon'
                    # dist = sum([e[1] - e[0] for e in self.exons[0:exon_n]])
                    #rf = self.find_reading_frame(dist_from_cds_start)
                    break
                else:
                    pass

                exon_n += 1


        self.feature, self.rf = feature, rf





##annotation files
def annotate(identify_file, annotate_path):

    res_data = pd.read_csv(identify_file, sep='\t')
    coords = list(res_data['Genomic Coordinate'])
    Transcript.load_transcripts(annotate_path,coords)
    res_data = res_data.rename(columns={'Chromosome': 'Chr'})
    Feature,Gene_Name,Distance  = list(), list(), list()

    for coord in coords:
        tx = Transcript.transcript(coord)

        if tx == 'intergenic':
            Feature.append('intergenic')
            Gene_Name.append('-')

        else:
            Feature.append(tx.feature)
            Gene_Name.append(tx.tx_info()[2])

    res_data['Gene_Name'] = Gene_Name
    #res_data['Gene_Type'] = Gene_Type
    #res_data['Gene_ID'] = Gene_ID
    res_data['Feature'] = Feature
    res_data = res_data.sort_values('Nuclease_Read_Count',ascending = False)
    outfile_name = identify_file.strip('.txt') +'_annotated.csv'
    print(outfile_name)
    res_data.to_csv(outfile_name, index = False)



##set input

# data_dir =os.path.dirname(os.path.realpath(__file__)) +"/changeseq/test/data/MergedOutput/identified/"
# annot_dir = '/groups/clinical/projects/clinical_shared_data/hg38/annotations/'
# data_dir = "/groups/clinical/projects/Assay_Dev/CHANGEseq/CS_06/identified/"
# annotate_path = '/groups/clinical/projects/editability/tables/processed_tables/ncbiRefSeq.bed.gz'
#  identify_files = [data_dir + x for x in os.listdir(data_dir) if x.endswith('_matched.txt')]
# identify_file = "/groups/clinical/projects/Assay_Dev/CHANGEseq/CS_06/identified/spCas9_78_ALT7_R2_06_identified_matched.txt"
