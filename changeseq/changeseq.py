#!/usr/bin/env python
#-*- coding: utf-8 -*-

"""
circleseq.py as the wrapper for CIRCLE-seq analysis
"""

import argparse
import os
import subprocess
import traceback
import log
from copy import deepcopy as dp
from pprint import pformat
## Project Modules
from alignReads import alignReads
from visualization import visualizeOfftargets
from utility import get_parameters
import findCleavageSites
from annotate import annotate
from report_qc import write_qc
from multi_sample_combiner import process_results
from fastq_qc import fastqQC

logger = log.createCustomLogger('root')


try:
	global p_dir
	module_file = findCleavageSites.__file__
	p_dir = os.path.dirname(os.path.realpath(__file__))
	logger.info(f"module file path: {module_file}")
except AttributeError:
	logger.info("findCleavageSites: (built-in module or not found)")


class CircleSeq:
	def __init__(self):
		self.parameters = {}
		self.output_dir = {}
		self.findCleavageSites_input_bam = {}
		self.vis_input_tsv = {}
		self.annotation_file = {}
		self.aligned = {}

	def parseManifest(self,analysis_folder,fq_dir,manifest,settings, sample='all'):
		logger.info('Loading manifest...')

		try:
			parameters = get_parameters(analysis_folder,fq_dir,manifest,settings)
			self.parameters = dp(parameters)

			if sample != 'all':
				self.parameters['samples'] = {}
				self.parameters['samples'][sample] = parameters['samples'][sample]

			# Make folders for output
			for folder in ['preprocessed','aligned', 'raw_results', 'fastq','qc',
                           'post-process_results','post-process_results/visualization','post-process_results/tables',
                           'raw_results/visualizations','raw_results/tables']:
				self.output_dir[folder] = os.path.join(self.parameters["analysis_folder"], folder)
				if not os.path.exists(self.output_dir[folder]):
					os.makedirs(self.output_dir[folder])

			# Just to initialize some default input file names for the identify and visualization steps
			for sample in self.parameters['samples']:
				self.findCleavageSites_input_bam[sample] = [f"{self.output_dir['aligned']}/{sample}.st.bam",
                                                            f"{self.output_dir['aligned']}/Control_{sample}.st.bam"]

			# initialize some default input file names for annotation files
			for sample in self.parameters['samples']:
				self.annotation_file[sample] = f"{self.output_dir[ 'raw_results/tables']}/{sample}_annotated_results.csv"

		except Exception as e:
			logger.error(
				'Incorrect or malformed manifest file. Please ensure your manifest contains all required fields.')
			logger.error(traceback.format_exc())

	def alignReads(self):
		"""BWA mapping
		"""
		logger.info('Aligning reads...')
		for sample in self.parameters['samples']:
			try:
				self.aligned[sample]=alignReads(trim=True, output_dir=self.output_dir['aligned'],
							   trimmed_fastq_output=self.output_dir['fastq'],
							   R1=self.parameters['samples'][sample]["read1"], R2=self.parameters['samples'][sample]["read2"],
							   label=sample, **self.parameters)
				control = alignReads(trim=True, output_dir=self.output_dir['aligned'],
										 trimmed_fastq_output=self.output_dir['fastq'],
						   R1=self.parameters['samples'][sample]["controlread1"], R2=self.parameters['samples'][sample]["controlread2"],
									 label="Control_"+sample,**self.parameters)

				logger.info('Finished aligning reads to genome.')

			except Exception as e:
				logger.error('Error aligning for sample %s.'%(sample))
				logger.error(traceback.format_exc())

	def findCleavageSites(self):
		logger.info('Identifying off-target cleavage sites...')
		for sample in self.parameters['samples']:
			try:
				# logger.info(sample)
				bam,control_bam = self.findCleavageSites_input_bam[sample]

				findCleavageSites.compare(bam=bam, control=control_bam, label=sample, output_dir=self.output_dir['raw_results/tables'],
							targetsite=self.parameters['samples'][sample]["target"],**self.parameters)

				logger.info('findCleavage Sites self.parameters: %s', pformat(self.parameters))
			except Exception as e:
				logger.error('Error findCleavageSites for sample %s.'%(sample))
				logger.error(traceback.format_exc())
				quit()

	def addAnnotations(self):
		for sample in self.parameters['samples']:
			try:
				matched_file = f"{self.output_dir[ 'raw_results/tables']}/{sample}_identified_matched.txt"
				print(f"Annotating {matched_file}")
				annotate(matched_file, self.parameters['annotate_path'])
			except Exception as e:
				logger.error('Error Annotating for sample %s.' % (sample))

	def visualize(self):
		logger.info('Visualizing off-target sites')

		for sample in self.parameters['samples']:
			try:
				infile = os.path.join(self.parameters["analysis_folder"],  'raw_results/tables',sample + '_identified_matched_annotated.csv')
				outfile = os.path.join(self.parameters["analysis_folder"],  'raw_results/visualizations',sample + '_offtargets.svg')
				visualizeOfftargets(infile, outfile, title=sample,PAM=self.parameters["PAM"])
			except Exception as e:
				logger.error('Error visualizing off-target sites: %s'%(sample))
				logger.error(traceback.format_exc())
		logger.info('Finished visualizing off-target sites')

	def analyze(self):
		logger.info('Normalizing and Joining reads')
		for rep_group_name,replicates in self.parameters['replicates'].items():
			infiles = []
			qcfiles = []
			for sample in self.parameters['replicates'][rep_group_name]['sample_name']:
				infiles.append(os.path.join(self.parameters["analysis_folder"], 'raw_results/tables',
											sample + '_identified_matched_annotated.csv'))
				qcfiles.append(os.path.join(self.parameters["analysis_folder"], 'qc', sample + '_qc_report.txt'))
			logger.info('Normalizing {rep_group_name}')

			outfolder = os.path.join(self.parameters["analysis_folder"],'post-process_results/')
			process_results(rep_group_name, replicates, infiles, qcfiles,
							outfolder=outfolder,
							normalization_method=self.parameters['normalize'],
							read_threshold = int(self.parameters['read_threshold']), PAM=self.parameters["PAM"])


	def QC(self):

		try:
			for sample in self.parameters['samples']:
				logger.info('Running fastq quality and adapter analysis for {0}'.format(sample))
				fastqc_logfile = os.path.join(self.parameters["analysis_folder"], 'preprocessed', sample + '.html')
				fastqQC(self.parameters['samples'][sample]['read1'],
                        self.parameters['samples'][sample]['read2'],
                        fastqc_logfile)
		except Exception as e:
			logger.error('Error with Fastqc')

		try:
			for sample in self.parameters['samples']:
				script_path = p_dir + "/QC_matched_alignment.sh"
				logger.info('Running alignment coverage QC for {0}'.format(sample))
				coverage_command = 'sh {0} {1} {2}'.format(script_path,sample, self.parameters["analysis_folder"])
				logger.info(coverage_command)
				subprocess.check_call(coverage_command, shell=True)
				logger.info('OT Coverage for {0} completed.'.format(sample))
		except Exception as e:
			logger.error('Error with alignment coverage QC ')
			logger.error('skipping....')

		try:
			for sample in self.parameters['samples']:
				logger.info('Running  report  for {0}'.format(sample))
				preprocessed_logfile =  os.path.join(self.parameters["analysis_folder"], 'fastq', sample + '_trim_log.txt')
				coverage_stat_file  = os.path.join(self.parameters["analysis_folder"], 'qc', sample + '_aligned_stats.txt')
				qc_file = os.path.join(self.parameters["analysis_folder"], 'qc', sample + '_qc_report.txt')

				write_qc(qc_file, preprocessed_logfile, coverage_stat_file)
		except Exception as e:
			logger.error('Error with QC stats report')
			logger.error('skipping..')

	def parallel(self, manifest_path, lsf, run='all'):
		logger.info('Submitting parallel jobs')
		current_script = __file__

		try:
			for sample in self.parameters['samples']:
				cmd = ' python {0} {1} --manifest {2} --sample {3}'.format(current_script, run, manifest_path, sample)
				logger.info(cmd)
				# subprocess.call(lsf.split() + [cmd])
				# print (lsf+cmd)
				subprocess.call(lsf + cmd,shell=True)
			logger.info('Finished job submission')

		except Exception as e:
			logger.error('Error submitting jobs.')
			logger.error(traceback.format_exc())

	def skip_align_parallel(self, manifest_path, lsf, run='skip_align'):
		logger.info('Submitting parallel jobs')
		current_script = __file__

		try:
			for sample in self.parameters['samples']:
				cmd = 'python {0} {1} --manifest {2} --sample {3}'.format(current_script, run, manifest_path, sample)
				logger.info(cmd)
				subprocess.call(lsf.split() + [cmd])
			logger.info('Finished job submission')

		except Exception as e:
			logger.error('Error submitting jobs.')
			logger.error(traceback.format_exc())

def parse_args():
	parser = argparse.ArgumentParser()

	subparsers = parser.add_subparsers(description='Individual Step Commands',
									   help='Use this to run individual steps of the pipeline',
									   dest='command')

	all_parser = subparsers.add_parser('all', help='Run all steps of the pipeline')
	all_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
	all_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
	all_parser.add_argument('--analysis_folder', '-dir', help='Directory in which all pipeline outputs will be saved (aabs path)', required=True)
	all_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
	all_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	skip_align_parser = subparsers.add_parser('skip_align', help='Run all steps of the pipeline')
	skip_align_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	skip_align_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)',required=True)

	parallel_parser = subparsers.add_parser('parallel', help='Run all steps of the pipeline in parallel')
	parallel_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
	parallel_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
	parallel_parser.add_argument('--analysis_folder', '-dir', help='Directory in which all pipeline outputs will be saved (aabs path)', required=True)
	parallel_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
	parallel_parser.add_argument('--lsf', '-l', help='Specify LSF CMD', default='bsub -R rusage[mem=32000] -P Genomics -q standard')
	parallel_parser.add_argument('--run', '-r', help='Specify which steps of pipepline to run (all, align, identify, visualize, variants)', default='all')

	skip_align_parallel_parser = subparsers.add_parser('skip_align_parallel', help='Run all steps of the pipeline in parallel')
	skip_align_parallel_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	skip_align_parallel_parser.add_argument('--lsf', '-l', help='Specify LSF CMD', default='bsub -R rusage[mem=200000] -P ABE -q priority')
	skip_align_parallel_parser.add_argument('--run', '-r', help='Specify which steps of pipepline to run (all, align, identify, visualize, variants)', default='skip_align')

	align_parser = subparsers.add_parser('align', help='Run alignment only')
	align_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
	align_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
	align_parser.add_argument('--analysis_folder', '-dir', help='Directory in which all pipeline outputs will be saved (aabs path)', required=True)
	align_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
	align_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	data_parser = subparsers.add_parser('makefiles', help='Combine and normalize replicates, produces vizualations')
	data_parser.add_argument('--ftp_path', '-f', help='RefSeq FTP Path. Must be in refseq.txt format', default="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz")
	data_parser.add_argument('--reset_output', '-r', help='Change path changeseq uses for stored annotation file',default=False)

	fa_parser = subparsers.add_parser('set_fasta', help='Sets the default genome fasta path so it is not needed in sample manifest')
	fa_parser.add_argument('-fasta', '-fa',help='Specify full fasta path. make sure this is unzipped')

	identify_parser = subparsers.add_parser("identify", help='Run identification only')
	identify_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
	identify_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
	identify_parser.add_argument('--analysis_folder', '-dir', help='Directory in which all pipeline outputs will be saved (aabs path)', required=True)
	identify_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
	identify_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	visualize_parser = subparsers.add_parser('visualize', help='Run visualization only')
	visualize_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
	visualize_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
	visualize_parser.add_argument('--analysis_folder', '-dir', help='Directory in which all pipeline outputs will be saved (abs path)', required=True)
	visualize_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
	visualize_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	variants_parser = subparsers.add_parser('variants', help='Run variants analysis only')
	variants_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
	variants_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
	variants_parser.add_argument('--analysis_folder', '-dir', help='Directory in which all pipeline outputs will be saved (abs path)', required=True)
	variants_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
	variants_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	analyzer_parser = subparsers.add_parser('analyze', help='Run analysis  only')
	analyzer_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
	analyzer_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)',required=True)
	analyzer_parser.add_argument('--analysis_folder', '-dir',
							help='Directory in which all pipeline outputs will be saved (aabs path)', required=True)
	analyzer_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)',default='default')
	analyzer_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	coverage_parser = subparsers.add_parser('qc', help='Run QC analysis of matched sites')
	coverage_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
	coverage_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
	coverage_parser.add_argument('--analysis_folder', '-dir', help='Directory in which all pipeline outputs will be saved (aabs path)', required=True)
	coverage_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
	coverage_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	reference_free_parser = subparsers.add_parser('reference-free', help='Run reference-free discovery only')
	reference_free_parser.add_argument('--manifest', '-m', help='The file name of the sample manifest(.csv)', required=True)
	reference_free_parser.add_argument('--raw_fastq_folder', '-fq', help='Directory where raw fastq file are stored (abs path)', required=True)
	reference_free_parser.add_argument('--analysis_dir', '-dir', help='Directory in which all pipeline outputs will be saved (aabs path)', required=True)
	reference_free_parser.add_argument('--settings', '-stg', help='The file name of the parameters and settings are (.csv)', default = 'default')
	reference_free_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	return parser.parse_args()

def main():
	args = parse_args()

	if args.command == 'all':
		c = CircleSeq()
		c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
		c.alignReads()
		c.findCleavageSites()
		c.addAnnotations()
		c.visualize()
		c.analyze()
		c.QC()
		# c.callVariants()
		# c.extract_deamination_position()
		# c.extract_outward_reads()
	elif args.command == 'skip_align':
		c = CircleSeq()
		c.parseManifest(args.manifest, args.sample)
		c.findCleavageSites()
		c.addAnnotations()
		c.visualize()
		# c.callVariants()
	elif args.command == 'parallel':
		c = CircleSeq()
		c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
		c.parallel(args.manifest, args.lsf, args.run)
	elif args.command == 'skip_align_parallel':
		c = CircleSeq()
		c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
		c.skip_align_parallel(args.manifest, args.lsf, args.run)
	elif args.command == 'align':
		c = CircleSeq()
		c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
		c.alignReads()
	elif args.command == 'identify':
		c = CircleSeq()
		c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
		c.findCleavageSites()
		c.addAnnotations()
		c.visualize()
	elif args.command == 'visualize':
		c = CircleSeq()
		c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
		c.addAnnotations()
		c.visualize()
	elif args.command == 'analyze':
		c = CircleSeq()
		c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
		c.analyze()
	elif args.command == 'variants':
		c = CircleSeq()
		c.parseManifest(args.manifest, args.sample)
		c.callVariants()
	elif args.command == 'qc':
		c = CircleSeq()
		c.parseManifest(args.analysis_folder,args.raw_fastq_folder,args.manifest,args.settings, args.sample)
		c.QC()
	elif args.command == 'makefiles':
		from make_annotation_table import makefiles
		makefiles(args.ftp_path, p_dir, args.reset_output)


if __name__ == '__main__':
	main()
