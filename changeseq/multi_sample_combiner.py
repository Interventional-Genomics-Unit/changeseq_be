import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from matplotlib_venn import venn2
from visualization import *

global colors
colors = ['#406061', '#96BA66', '#B45C20','#718CA0',"#5C842F"]

def check_file(file):
    if os.path.isfile(file):
        pass
    else:
        logger.info("{file} is not your directory")

def get_read_depth(qcfiles):
    depths = []
    for file in qcfiles:
        check_file(file)
        line = open(file).readline()
        if 'total_reads:' in line:
            depths.append(int(line.split(':')[-1].strip().replace(",","")))
        else:
            logger.info('Total read depth for {file} could not be found for tpm normalization')
            logger.info('Re-run QC or choose another normalization method')
            exit()
    return depths

def find_mm_and_bulge(df):
    subs,insertions,deletions = [],[],[]
    for i,data in df.iterrows():
        db = 0
        rb = 0
        mm=0
        try:
            mm = int(data['Site_Substitution_Number'])

        except ValueError:
            for i,j in zip(data['Site_Sequence_Gaps_Allowed'][:-3],data['Aligned_Target_Sequence'][:-3]):
                if i==j:
                    pass
                elif i == "-":
                    rb = 1
                elif j =='-':
                    db = 1
                else:
                    mm+=1
        subs.append(mm)
        insertions.append(rb)
        deletions.append(db)
    return subs,insertions,deletions


def parse_df(df,sample):

    ## Filling in Aligned and ungapped sequence
    df.loc[:,['Aligned_Site_Sequence']]= df['Site_Sequence_Gaps_Allowed'].fillna(df['Site_Sequence'])
    df.loc[:,'Site_Sequence'] = df['Site_Sequence'].fillna(df['Site_Sequence_Gaps_Allowed'].str.replace("-", ""))

    df.loc[:, ['Aligned_Target_Sequence']] = df['Realigned_Target_Sequence'].replace('none', np.nan)
    df.loc[:,['Aligned_Target_Sequence']] = df.loc[:,'Aligned_Target_Sequence'].fillna(df['Target_Sequence'])
    df.loc[:,['RNA_Bulge']] = df['RNA_Bulge'].fillna(0)
    df.loc[:,['DNA_Bulge']] = df['DNA_Bulge'].fillna(0)
    df.loc[:,'Genomic Coordinate'] = df['Genomic Coordinate'].str.split('-',expand=True)[0] + df['Strand']
    df = df.loc[:,['Genomic Coordinate','Nuclease_Read_Count', 'Control_Read_Count', 'Site_Sequence',
             'Site_Substitution_Number','Aligned_Site_Sequence',
             'DNA_Bulge', 'RNA_Bulge', 'Target_Sequence','Aligned_Target_Sequence',
             'Gene_Name', 'Feature','Cell','MappingPositionStart',
       'MappingPositionEnd', 'WindowSequence']].copy()
    df = df.rename(columns={'Nuclease_Read_Count': 'Nuclease_Read_Count.' + sample,
                                         'Control_Read_Count':'Control_Read_Count.' +sample,
                                'Description':'Description.' +sample})

    return df


def join_replicates(sample1_df, sample2_df, suffixes):
    df = sample1_df.merge(sample2_df, on=['Genomic Coordinate', 'Site_Sequence'], how='outer',
                          suffixes=suffixes)

    keep_cols = ['Genomic Coordinate', 'Nuclease_Read_Count' + suffixes[0], 'Nuclease_Read_Count'+ suffixes[1],
                 'Control_Read_Count'+ suffixes[0], 'Control_Read_Count'+ suffixes[1],'Site_Sequence',
                 'Aligned_Site_Sequence','Site_Substitution_Number',
                 'DNA_Bulge', 'RNA_Bulge', 'Target_Sequence', 'Aligned_Target_Sequence', 'Gene_Name', 'Feature', 'Cell',
                 'MappingPositionStart','MappingPositionEnd', 'WindowSequence'
                 ]

    for col in df.columns:
        if "Read_Count" in col:
            df.loc[:,col] = df[col].fillna(0)

    for col in keep_cols[keep_cols.index( 'Aligned_Site_Sequence'):]:
        df.loc[:,col] = df[f'{col}{suffixes[0]}'].fillna(df[f'{col}{suffixes[1]}'])
        df = df.drop(columns=[f'{col}{suffixes[0]}', f'{col}{suffixes[1]}'])

    nuclease_read_count_cols = [col for col in df.columns if "Nuclease_Read_Count" in col]
    df.loc[:,'Number of Replicates Sites found'] = (df[nuclease_read_count_cols] > 0).sum(1)

    return df



def swarm_plot(joined_normalized, name, figout):
    sns.set_style('white')
    columns = [col for col in joined_normalized.columns if col.startswith('Nuc')]
    count_dict = joined_normalized.loc[:,joined_normalized.columns.str.startswith('Nuc')].to_dict('list')

    # Finding the mean reading excluding 0 sites
    mean = []
    counts = np.array(list(count_dict.values()))
    for i in range(len(counts[0,:])):
        count = counts[:,i]
        x = sum(count[count>0])
        mean.append(x/sum(count>0))

    df = joined_normalized.copy()
    df.loc[:,'Log2 Mean Reads'] = np.log2(np.array(mean))
    df.loc[:,'Guide Name'] = [name] * len(df['Log2 Mean Reads'])
    df.loc[:,'Shared'] = [f"found in {x} samples" for x in df['Number of Replicates Sites found']]

    # Highlight on target if present
    ont = df.loc[:,['Site_Substitution_Number', 'RNA_Bulge', 'DNA_Bulge']].sum(1)==0
    df.loc[ont, 'Shared'] = 'on target'
    df = df.sort_values('Shared', ascending=False)
    palette = colors[0:len(columns)]
    if 'on target' in list(df.Shared):
        palette = ['#B45C20'] + palette

    plt.figure(figsize=(6, 6))
    g = sns.swarmplot(data=df,
                      x='Guide Name',
                      y='Log2 Mean Reads',
                      hue='Shared',
                      palette=palette
                      )
    plt.yticks(np.arange(np.log2(10), np.log2(100000), step=np.log2(10)), [10, 100, 1000, 10000, 100000])
    plt.title(name)
    plt.ylabel("log2 read counts")
    plt.xlabel("")
    legend = g.get_legend()
    if legend:
        legend.set_title("")
    # plt.show()
    plt.savefig(figout)
    plt.close(figout)

def calc_jaccard(Rep1_unique, Rep2_unique, shared):
    total = Rep1_unique + Rep2_unique + shared
    return round(float(shared) / float(total) * 100, 2)

def create_pseudo_sample(data_with_pseudocount):
    # Compute geometric means
    geometric_means = np.exp(
        np.mean([np.log(counts) for counts in data_with_pseudocount.values()], axis=0)
    )
    return geometric_means


def median_normalization(count_dict):
    epsilon = 1
    # https://divingintogeneticsandgenomics.com/post/details-in-centered-log-ratio-clr-normalization-for-cite-seq-protein-count-data/

    data_with_pseudocount = {sample: np.array(counts) + epsilon for sample, counts in count_dict.items()}
    pseudo_sample = create_pseudo_sample(data_with_pseudocount)  # (log_cnt1,log_cnt2)
    ratios = {sample: np.array(counts) / pseudo_sample for sample, counts in data_with_pseudocount.items()}

    # Compute scaling factors
    scaling_factors = {
        sample: np.median(ratios[sample])
        for sample in count_dict}
    # Normalize counts
    normalized_data = {
        sample: (counts / scaling_factors[sample]).round(0)
        for sample, counts in count_dict.items()
    }
    return normalized_data

def rpm(count_dict,depths):
    # uses a total count scaling method

    min_depth = np.min(depths)

    scaling_factors = {
        sample: min_depth / depths[i]
        for i,sample in enumerate(list(count_dict.keys()))
    }

    normalized_data = {
        sample: (np.array(counts) * scaling_factors[sample]).round(0)
        for sample, counts in count_dict.items()
    }
    return normalized_data


def scatter_plot(x1,x2,name,figout):
    x = np.log2(np.array(x1)+1)
    y = np.log2(np.array(x2)+1)

    pearR = np.corrcoef(x1, x2)[1, 0]
    plt.figure(figsize=(4, 4))
    plt.scatter(x[x * y>0], y[x * y>0],color = colors[4],label="rho= %s" % (round(pearR,3)))
    plt.xticks(np.arange(np.log2(10), np.log2(100000), step=np.log2(10)),[10,100,1000,10000,100000])
    plt.yticks(np.arange(np.log2(10), np.log2(100000), step=np.log2(10)),[10,100,1000,10000,100000])
    plt.legend(loc=1)
    plt.title(name)
    #plt.show()
    plt.savefig(figout, bbox_inches='tight')
    plt.close(figout)

def vennplot_replicates(x1,x2,sample1,sample2,figout):
    font_size = 8
    x1, x2 =np.array(x1),np.array(x2)
    Rep1_unique,Rep2_unique  = sum(x1>0),sum(x2 > 0)
    shared = sum(x1 * x2 > 0)
    ja= calc_jaccard(Rep1_unique,Rep2_unique,shared)
    values = (Rep1_unique, Rep2_unique,shared)

    plt.figure(figsize=(4, 4))
    v =venn2(subsets = values,
             set_labels=(sample1,sample2),
             set_colors=(colors[1],colors[4]),alpha = 0.5)
    for l in ["A","B"]:
        v.get_label_by_id(l).set_fontsize(font_size)
        v.get_label_by_id(l).set_y(0.6)
        v.get_label_by_id(l).set_x(len(sample1)/100.0-0.4)

    plt.annotate("% replicate rites overlap " + str(ja), xy=v.get_label_by_id('010').get_position() +
                                           np.array([0, -0.5]), xytext=(-60, -30), ha='center',
                 textcoords='offset points')
    plt.savefig(figout, bbox_inches='tight')
    #plt.show()
    return ja
def clean_normalized(joined_normalized,read_threshold):
    read_columns = joined_normalized.columns.str.startswith('Nuc')

    # remove sites below threshold post normalization
    joined_normalized.loc[:, read_columns] = joined_normalized.loc[:,read_columns].applymap(
        lambda x: 0 if x and x < read_threshold else x)
    keep_rows = joined_normalized.loc[:, read_columns].sum(1) != 0
    joined_normalized = joined_normalized.loc[keep_rows, :]
    joined_normalized.loc[:,'Number of Replicates Sites found'] = (joined_normalized.loc[:, read_columns] > 0).sum(1)
    joined_normalized.loc[:,'Percent Total Reads'] =joined_normalized.loc[:, read_columns].sum(1) / joined_normalized.loc[:, read_columns].sum(1).sum()
    joined_normalized.loc[:, 'Percent Total Reads'] =  joined_normalized['Percent Total Reads'].round(5) *100
    joined_normalized.loc[:, 'Distance'] = joined_normalized.loc[:,['Site_Substitution_Number', 'RNA_Bulge', 'DNA_Bulge']].sum(1)

    # make a simplified version
    read_cols = [x for x in joined_normalized.columns if 'Nuc' in x]
    keep_cols = ['Genomic Coordinate'] + read_cols + ['Aligned_Site_Sequence',
                                                      'Site_Substitution_Number', 'DNA_Bulge', 'RNA_Bulge',
                                                      'Number of Replicates Sites found', 'Percent Total Reads',
                                                      'Aligned_Target_Sequence', 'Gene_Name', 'Feature', 'Cell']

    joined_simplified_report =  joined_normalized.loc[:, keep_cols ]
    return joined_normalized, joined_simplified_report


def normalize(joined,qcfiles,normalization_method,read_threshold):

    count_dict = joined.loc[:,joined.columns.str.startswith('Nuc')].to_dict('list')
    if normalization_method == "median":
        norm_count_dict = median_normalization(count_dict)
    elif normalization_method == "rpm":
        depths =get_read_depth(qcfiles)
        norm_count_dict =rpm(count_dict,depths)
    else:
        norm_count_dict = count_dict
    for k,v in norm_count_dict.items():
        joined.loc[:, k] = v
    joined_normalized, simplified_report = clean_normalized(joined,read_threshold)

    return joined_normalized,simplified_report

def make_offtarget_dict(joined_normalized,subset):
    sample_df = joined_normalized.loc[joined_normalized[subset]>0].reset_index().copy()
    offtargets = []
    total_seq = sample_df.shape[0]
    for i,row in sample_df.iterrows():
        offtarget_reads = row[subset]
        annot = ""
        target_seq = row['Target_Sequence'].replace("-", "")
        if int(row['RNA_Bulge']) + int(row['DNA_Bulge']) == 0:

            no_bulge_offtarget_sequence = row['Site_Sequence']

            bulge_offtarget_sequence = ""
            realigned_target_seq = row['Aligned_Target_Sequence']

        else:
            no_bulge_offtarget_sequence = ""
            realigned_target_seq =  row['Aligned_Target_Sequence']
            bulge_offtarget_sequence = row['Aligned_Site_Sequence']

        coord = row['Genomic Coordinate']

        try:
            if "intergenic" not in row['Feature']:
                annot = row['Gene_Name'] + "," + row['Feature'].replace("non-coding RNA","ncRNA")  # gene name and feature

        except:
            pass

        if no_bulge_offtarget_sequence != '' or bulge_offtarget_sequence != '':
            if no_bulge_offtarget_sequence:
                total_seq += 1
            if bulge_offtarget_sequence:
                total_seq += 1
            offtargets.append({'seq': no_bulge_offtarget_sequence.strip(),
                               'bulged_seq': bulge_offtarget_sequence.strip(),
                               'reads': int(offtarget_reads),
                               'target_seq': target_seq.strip(),
                               'coord': coord,
                               'annot': str(annot),
                               'realigned_target_seq': realigned_target_seq.strip()
                               })
    offtargets = sorted(offtargets, key=lambda x: x['reads'], reverse=True)
    return offtargets, target_seq, total_seq

def process_results(rep_group_name,replicates,infiles,qcfiles,outfolder, normalization_method,read_threshold,PAM):

    processed_outfile =outfolder + "/tables/"+ rep_group_name +'_joined.csv'
    simplified_report_outfile = processed_outfile.replace('.csv','_simplified.csv')
    swarm_plot_out = outfolder + "/visualization/"+ rep_group_name + "_postprocess_swarmplot.png"

    first_file = True
    for i,infile in enumerate(infiles):
        sample = replicates['sample_name'][i]
        df = pd.read_csv(infile)
        df = parse_df(df,sample)

        if first_file:
            joined = df
            previous_suffix = '.' + sample
            first_file = False
        else:
            suffix = '.' +sample
            joined = join_replicates(joined, df, suffixes = [previous_suffix,suffix])
            previous_suffix = suffix


    joined_normalized, simplified_report = normalize(joined,qcfiles,normalization_method,read_threshold)
    joined_normalized.to_csv(processed_outfile, index = False)
    simplified_report.to_csv(simplified_report_outfile, index = False)

    ## plotting
    swarm_plot(joined_normalized, rep_group_name, swarm_plot_out)

    for sample in replicates['sample_name']:
        offtargets, target_seq, total_seq = make_offtarget_dict(joined_normalized,subset='Nuclease_Read_Count.' + sample)
        alignment_plot = outfolder +"/visualization/"+ sample.replace(" ", "_") + "_postprocess_alignment_plot.svg"
        draw_plot(target_seq, offtargets, total_seq, outfile=alignment_plot, title=sample, PAM=PAM)

    for i in range(len(replicates['sample_name'])-1):
        sample_1= replicates['sample_name'][i]  #'Nuclease_Read_Count.' + sample
        for j in range(1,len(replicates['sample_name'])):
            sample_2 = replicates['sample_name'][j]
            x1, x2 = list(joined_normalized[f'Nuclease_Read_Count.{sample_1}']), list(
                joined_normalized[f'Nuclease_Read_Count.{sample_2}'])
            scatter_out = f"{outfolder}/visualization/{sample_1}_&_{sample_2}_postprocess_scatterplot.png"
            venn_out = f"{outfolder}/visualization/{sample_1}_&_{sample_2}_postprocess_venn.png"

            scatter_plot(x1, x2,f"{sample_1} & {sample_2}", scatter_out)
            sim = vennplot_replicates(x1,x2,sample_1, sample_2, venn_out)


def parse_args():
    mainParser = argparse.ArgumentParser()
    mainParser.add_argument('--replicate_sample_names', '-r', help='replicate sample names as indicated on manifest, comma seperated ex. CCR5_rep1,CCR5_rep2')
    mainParser.add_argument('--output', '-o', help='output directory')
    mainParser.add_argument('--name', '-n', help='Sample name for labeling graphs')
    mainParser.add_argument('--infiles', '-f', help='list of i_identified_matched_annotated.csv files in raw_output folder ')
    mainParser.add_argument('--normalization_method', '-nm',
                            help='')
    mainParser.add_argument('--analysis_folder', '-p',
                            help='')
    mainParser.add_argument('--PAM', '-p',
                            help='')
    mainParser.add_argument('--read_threshold', '-rt', help='limit sample combining to a min read count threshold',
                            default=6)

    return mainParser.parse_args()


def main():
    args = parse_args()
    try:
        replicates = {'sample_name': args.replicate_sample_names.split(",")}
        process_results(args.name,replicates,args.infiles,qcfiles,args.analysis_folder, args.normalization_method,args.read_threshold,PAM)
    except Exception as e:
        logger.error('Error combined replicates')
