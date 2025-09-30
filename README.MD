
[![Version][version-shield]][version-url]
[![Python versions][python-shield]][python-url]
[![Platforms][platform-shield]][python-url]


# CHANGE-seq-BE enables simultaneously sensitive and unbiased in vitro profiling of base editor genome-wide activity

This is a repository for CHANGE-seq-BE analytical software, which takes sample-specific paired-end FASTQ files as input and produces a list of CHANGE-seq-BE detected off-target cleavage sites as output.

# IGI CHANGE-seq-BE

adaptions made by T. Hudson to CHANGE-seq-BE for Off-Target Site Nomination https://github.com/tsailabSJ/changeseq

IGI Added Features:


* simple .csv sample manifest
* replicate combination
* normalization of replicates - (median-of-ratios or total count scaling)
* QC metrics of alignments and fastq circularization read composition
* identified site bam file for use with IGV
* genomic annotation to sites
* more plotting options - venn diagram, swarmplot and more
* Fixed variable analysis bug


# Summary

This package implements a pipeline that takes in reads from the CHANGE-seq-BE assay and returns detected cleavage sites as output. The individual pipeline steps are:

1. **Tn5 adapter trimming**: `cutadapt` remove adapter CTGTCTCTTATACACATCTACGTAGATGTGTATAAGAGACAG
2. **Read Alignment**: `bwa mem` paired end mapping
3. **Read Filtering**: (1) `mapq_threshold` (2) only outward, self-overlapped PE reads `BEmodel_min_overlap` and `BEmodel_max_overlap`. 
4. **Cleavage Site Identification**: Mapped sites are analyzed to determine which represent high-quality cleavage sites. (Same algorithm as before, peak is identified in a `window_size`)
5. **Off-target Identification**: Off-targets are identified by extending the cleavage site +/- `search_radius` and then fuzzy alignment to the on-target sequence `mismatch_threshold`. Gap penalty is 3 and mismatch penalty is 1. Output files are two. `identified_matched.txt` contains off-target sequences similar to on-target. `identified_unmatched.txt` contains off-target coordinates, sequences are not similar to on-target.
6. **Visualization of Results**: Identified on-target and off-target cleavage sites are rendered as a color-coded alignment map for easy analysis of results. (Same algorithm as before)

# Installation

The easiest way to install is git clone into a designated folder

```

git clone https://github.com/Interventional-Genomics-Unit/changeseq_be

python changeseq/changeseq/changeseq.py -h

## All dependencies are installed, e.g., BWA 0.7.17, samtools 1.9

```
 Create a conda environment and import all required packages in one line below:
 
 ```
 conda env create -f cseq.yaml
 conda activate cseq
 ```

In order to set up the annotation table for annotating OT sites run makefiles to download RefSeq annotations. Additional commands can change the location and FTP link (--outdir --ftp_path)

```
cd /your-path/changeseq_be/
python changeseq/changeseq.py makefiles
```
 

# Usage
The Tsai lab change-seq pipeline requires a manifest yaml file specifying input files, output directory, and pipeline parameters. 
In the IGU version this file is split into t two csv files 1) a sample manifest file <code>--manifest</code> 2) an optional settings 
file <code>--settings</code>
Once the csv files are created, users can simply run the changeseq.py with the parameters indicating the 
raw fastq directory <code>--raw_fastq_folder</code> and main analysis directory <code>--analysis_folder</code>


The bash wrapper is initiated with the parameters.csv, manifest.csv, command, and sample_name

```
cd /your-path/IGU_CHANGEseq/
python changeseq/changeseq.py all --analysis_folder /your-path/IGU_CHANGEseq/changeseq/test/Standard_Output --raw_fastq_folder /your-path/IGU_CHANGEseq/changeseq/test/input --settings unmerged_parameters.csv --manifest manifest.csv --sample all
```



# Writing A Manifest File

The original manifest.yaml is no longer used in this version. The prior manifest.yaml is now split into two .csv input files. 
The manifest file contains the sample information and an example is shown below:

 sequencing_sample_name | sample_name | target  | description | replicate_group_name |
|------------------------|-------------|---------|-------------|-------------|
| TET2_rep1              | TET2_rep1   | NTC_sample | CCAGCATAAGAGCATGAATGNGN | TET2 |
| TET2_rep2              | TET2_rep2   | NTC_sample | CCAGCATAAGAGCATGAATGNGN | TET2 |

- `sequencing_sample_name`: sample name as indicated in the raw fastq (Leave out the ("_R1_001")
- `sample_name`: The sample name as you wish to see on the output file names  
- `control_sequencing_sample_name`: No RNP control sample name as indicated in the raw fastq (Leave out the ("_R1_001")
- `target`: Target sequence for that sample. Accepts degenerate bases.
- `description`: A brief description of the sample. Usually the cell or gDNA Source name
- `replicate_group_name`: (optional) an identifier that designates replicate samples. This must be present to run normalization

 # Parameters and Setting
 This file is a .csv with the first column the parameter name and the second column the parameter name. 
 The parameters are now stored internally and only need to be given if different than default
 See IGU_CHANGEseq/changeseq/test/parameters.csv

- `bwa`: The absolute path to the `bwa` executable
- `samtools`: The absolute path to the `samtools` executable
- `cutadapt`: The absolute path to the `cutadapt` executable
- `BEmodel_min_overlap`: 2
- 'BEmodel_max_overlap': 15
- `read_threshold`: The minimum number of reads at a location for that location to be called as a site. We recommend leaving it to the default value of 4.
- `window_size`: Size of the sliding window, we recommend leaving it to the default value of 3.
- `mapq_threshold`: Minimum read mapping quality score. We recommend leaving it to the default value of 40.
- `start_threshold`: Tolerance for breakpoint location. We recommend leaving it to the default value of 1.
- `gap_threshold`: Distance between breakpoints. We recommend leaving it to the default value of 3 for Cas9.
- `mismatch_threshold`: Number of tolerated gaps in the fuzzy target search setp. We recommend leaving it to the default value of 6.
- `read_length`: Fastq file read length, default is 151.
- `PAM`: PAM sequence, default is NGG.
- `normalize`: normalization method to use for replicates "none", "median" or "rpm", default=median

# Commands

### Pipeline commands (must be run alongside samples command)

- `all`: runs all of the pipelines
- `align`: merges (if indicated in parameters.csv), trims Tn5 adapters and aligns reads to ref genome
- `identify`: identifies nominated off-target sites and annotates sites
- `visualize`:annotates identified sites and creates alignment .svg image
- `qc`: create .bam file for identified sites for IGV, alignment histogram and alignment stats file and additional QC about fastq 
- `analyze`: joines replicates if given, normalizes reads and creates additional plots
- `variants`: *currently* updating

### Parameters

- `analysis_folder`: The absolute path to the folder in which all pipeline outputs will be saved.
- `raw_fastq_folder`:The absolute path to the folder in which all raw fastq files are stored.
- `all`: runs all samples in manifest
- `sample_name`: runs the specific sample indicated

# Pipeline Output
When running the full pipeline, the results of each step are outputted to the `output_folder` in a separate folder for each step. The output folders and their respective contents are as follows:

- `output_folder/aligned`: Contains an alignment `.sam`, alignment `.bam`, sorted `bam`, and `.bai` index file for each sample.
- `output_folder/fastq`: Merged `.fastq.gz` files for each sample.
- `output_folder/identified`: Contains tab-delimited `.txt` files for each sample containing the identified DSBs, control DSBs, filtered DSBs, and read quantification.
- `output_folder/visualization`: Contains a `.svg` vector image representing an alignment of all detected off-targets to the targetsite for each sample.

# FAQ

None yet, we will keep this updated as needed.

[version-shield]: https://img.shields.io/conda/v/tsailabsj/changeseq.svg
[version-url]: https://anaconda.org/tsailabSJ/changeseq
[python-shield]: https://img.shields.io/pypi/pyversions/changeseq.svg
[python-url]: https://pypi.python.org/pypi/changeseq
[platform-shield]: https://anaconda.org/tsailabsj/changeseq/badges/platforms.svg
