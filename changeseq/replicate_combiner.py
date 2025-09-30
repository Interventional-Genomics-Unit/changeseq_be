import os
import pandas as pd
import numpy as np
import matplotlib
#matplotlib.use('Agg')  # Use a non-GUI backend
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

import visualization
from visualization import *

global colors
colors = {'sample1': '#8FBC8F', 'sample2': '#029386'}#, 'T': '#D8BFD8', 'C': '#8FBC8F', 'N': '#AFEEEE', 'R': '#3CB371', '-': '#E6E6FA'}

def check_file(file):
    if os.path.isfile(file):
        pass
    else:
        print(file + " is not in the your-CHANGEseq-analysis-path/identified/ directory")
        
def find_mm_and_bulge(df):
    subs,insertions,deletions = [],[],[]
    for i,data in df.iterrows():
        db = 0
        rb = 0
        mm=0
        try:
            mm = int(data['Site_Substitution_Number'])

        except ValueError:
            for i,j in zip(data['Site_Sequence_Gaps_Allowed'][:-3],data['Realigned_Target_Sequence'][:-3]):
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



def parse_df(file, threshold = 6):
    df = pd.read_csv(file)
    df = df[df['Nuclease_Read_Count'] >= threshold]

    ## Removing Gaps and merging
    df['Aligned_Site_Sequence'] = df['Site_Sequence'].fillna(df['Site_Sequence_Gaps_Allowed'])
    df['Aligned_Site_Sequence'] = df['Site_Sequence'].fillna(df['Site_Sequence_Gaps_Allowed'])
    df = df[['Genomic Coordinate','Site_Sequence','Nuclease_Read_Count','Control_Read_Count',
        'Aligned_Site_Sequence','Target_Sequence','Site_Substitution_Number',
        'DNA_Bulge','RNA_Bulge','Gene_Name','Feature']]

    return df

def join_replicates(sample_files,threshold = 6):
    keep_cols = ['Genomic Coordinate','Nuclease_Read_Count.Rep1','Nuclease_Read_Count.Rep2',
           'Control_Read_Count.Rep1','Control_Read_Count.Rep2','Site_Sequence',
        'Aligned_Site_Sequence','Target_Sequence','Site_Substitution_Number',
        'DNA_Bulge','RNA_Bulge','Gene_Name','Feature']

    sample1_df = parse_df(sample_files[0], threshold)
    sample2_df = parse_df(sample_files[1], threshold)
    df = sample1_df.merge(sample2_df, on=['Genomic Coordinate', 'Site_Sequence'], how='outer',
                          suffixes=[".Rep1", ".Rep2"])

    # df = sample1_df.join(sample2_df.set_index('Site_Sequence_NoGaps'),on = 'Site_Sequence_NoGaps',how = 'outer',lsuffix=".Rep1",rsuffix=".Rep2")
    for col in df.columns:
        if "Read_Count" in col:
            df[col] = df[col].fillna(0)


    for col in keep_cols[keep_cols.index('Aligned_Site_Sequence'):]:
        df[col] = df[f'{col}.Rep1'].fillna(df[f'{col}.Rep2'])
        df = df.drop(columns=[f'{col}.Rep1', f'{col}.Rep2'])

    nuclease_read_count_cols = [col for col in df.columns if "Nuclease_Read_Count" in col]
    df['Number of Replicates Sites found'] = (df[nuclease_read_count_cols]>0).sum(1)


    return df

def scatter_by_overlap(x1,x2,name):
    x = np.log2(np.array(x1) + 1)
    y = np.log2(np.array(x2) + 1)
    pearR = np.corrcoef(x1, x2)[1, 0]
    colors= ['#8FBC8F' if j != 0 and k !=0 else 'gray' for j,k in zip(x,y) ]
    plt.figure(figsize=(4, 4))

    plt.scatter(x, y, c=colors, label="r= %s" % (round(pearR, 3)))

    #plt.scatter(x[x * y > 0], y[x * y > 0], color=colors.values()[1], label="r= %s" % (round(pearR, 3)))
    plt.xticks(np.arange(np.log2(1), np.log2(100000), step=np.log2(10)), [0,10, 100, 1000, 10000, 100000])
    plt.yticks(np.arange(np.log2(1), np.log2(100000), step=np.log2(10)), [0,10, 100, 1000, 10000, 100000])
    plt.legend(loc=1)
    #A = (x + y) / 2
    #M = x - y
    #plt.scatter(x=A, y=M, s=10, c = colors)  # s is point size
    plt.title(name)
    #plt.show()

def swarm_plot(df,name,figout):
    x1, x2 = list(df['Nuclease_Read_Count.Rep1']), list(df['Nuclease_Read_Count.Rep2'])
    y = [(x + y) / 2 if (x * y) > 0 else x + y for x, y in zip(x1, x2)]
    y = np.log2(np.array(y))
    df['Log2 Mean Reads'] = y
    df['Guide Name'] = [name] * len(df['Log2 Mean Reads'])
    df['Shared'] = [True if (i * j) > 0 else False for i,j in zip(x1,x2)]
    palette=['lightgrey','royalblue']

    # Highlight on target if present

    if len(df.loc[df['Site_Substitution_Number'] + df['RNA_Bulge'].fillna(0) + df['DNA_Bulge'].fillna(0) ==0,'Shared']):
        df.loc[(df['Site_Substitution_Number'] + df['RNA_Bulge'].fillna(0) + df['DNA_Bulge'].fillna(0) == 0), 'Shared'] = 'on target'
        df = df.sort_values('Shared',ascending = False)
        palette = ['#6D091F','royalblue', 'lightgrey']
    plt.figure(figsize = (5,5))
    g = sns.swarmplot(data = df,
                      x= 'Guide Name',
                      y= 'Log2 Mean Reads',
                      hue='Shared',
                      palette=palette
                      )
    plt.yticks(np.arange(np.log2(10), np.log2(100000), step=np.log2(10)), [10, 100, 1000, 10000, 100000])
    plt.title(name)
    plt.ylabel("Read Counts")
    plt.xlabel("")
    #plt.show()
    plt.savefig(figout, bbox_inches='tight')
    plt.close(figout)


def scatter_plot(x1,x2,name,figout):
    x = np.log2(np.array(x1)+1)
    y = np.log2(np.array(x2)+1)

    pearR = np.corrcoef(x1, x2)[1, 0]
    plt.figure(figsize=(4, 4))
    plt.scatter(x[x * y>0], y[x * y>0],color = list(colors.values())[1],label="rho= %s" % (round(pearR,3)))
    plt.xticks(np.arange(np.log2(10), np.log2(100000), step=np.log2(10)),[10,100,1000,10000,100000])
    plt.yticks(np.arange(np.log2(10), np.log2(100000), step=np.log2(10)),[10,100,1000,10000,100000])
    plt.legend(loc=1)
    plt.title(name)
    #plt.show()
    plt.savefig(figout, bbox_inches='tight')
    plt.close(figout)

def make_offtarget_dict(joined_normalized,subset):
    sample_df = joined_normalized.loc[joined_normalized[subset]>0].reset_index().copy()
    offtargets = []
    total_seq = sample_df.shape[0]
    for i,row in sample_df.iterrows():
        offtarget_reads = row[subset]
        annot = ""
        if int(row['RNA_Bulge']) + int(row['DNA_Bulge']) == 0:
            no_bulge_offtarget_sequence = row['Aligned_Site_Sequence']
            target_seq = row['Target_Sequence']
            bulge_offtarget_sequence = ""
            realigned_target_seq = ""
        else:
            bulge_offtarget_sequence = row['Aligned_Site_Sequence']
            realigned_target_seq = row['Target_Sequence']
            no_bulge_offtarget_sequence = ""

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
                               'overlaps': 'x',
                               'realigned_target_seq': realigned_target_seq.strip()
                               })
    offtargets = sorted(offtargets, key=lambda x: x['reads'], reverse=True)
    return offtargets, target_seq, total_seq

def calc_jaccard(Rep1_unique,Rep2_unique,shared):
    total = Rep1_unique + Rep2_unique + shared
    return round(float(shared)/float(total)*100,2)

def create_pseudo_sample(data_with_pseudocount):
    # Compute geometric means
    geometric_means = np.exp(
        np.mean([np.log(counts) for counts in data_with_pseudocount.values()], axis=0)
    )

    return geometric_means

def median_normalization(read_counts1,read_counts2):
    # Similar to DESeq2
    #https://divingintogeneticsandgenomics.com/post/details-in-centered-log-ratio-clr-normalization-for-cite-seq-protein-count-data/
    epsilon = 1
    count_dict = {'Sample1':read_counts1,'Sample2':read_counts2}
    data_with_pseudocount = {sample: np.array(counts) + epsilon for sample, counts in count_dict.items()}
    pseudo_sample = create_pseudo_sample(data_with_pseudocount)#(log_cnt1,log_cnt2)
    ratios = {sample: np.array(counts) / pseudo_sample for sample, counts in data_with_pseudocount.items()}

    # Compute scaling factors
    scaling_factors = {
        sample: np.median(ratios[sample])
        for sample in count_dict
    }

    # Normalize counts
    normalized_data = {
        sample: (counts / scaling_factors[sample]).round(0)
        for sample, counts in count_dict.items()
    }
    new_read_counts1, new_read_counts2 = normalized_data['Sample1'],normalized_data['Sample2']

    return new_read_counts1, new_read_counts2

def normalize(joined,threshold=6):
    read_counts1, read_counts2 = list(joined['Nuclease_Read_Count.Rep1']), list(joined['Nuclease_Read_Count.Rep2'])
    rep1,rep2 = median_normalization(read_counts1, read_counts2)

    joined_normalized = joined.copy()
    joined_normalized['Nuclease_Read_Count.Rep1'] = [x if x >=threshold else 0 for x in rep1]
    joined_normalized['Nuclease_Read_Count.Rep2'] = [x if x >= threshold else 0 for x in rep2]
    joined_normalized = joined_normalized.loc[(joined_normalized['Nuclease_Read_Count.Rep1'] + joined_normalized['Nuclease_Read_Count.Rep2'] !=0),:]
    #add replicate
    nuclease_read_count_cols = [col for col in joined_normalized.columns if "Nuclease_Read_Count" in col]
    joined_normalized['Number of Replicates Sites found'] = (joined_normalized[nuclease_read_count_cols]>0).sum(1)
    mean_counts=joined_normalized[nuclease_read_count_cols].apply(lambda x: (np.exp(np.mean(np.log(x+1)))-1).round(0),axis=1)
    i = list(joined_normalized.columns).index(nuclease_read_count_cols[0])
    joined_normalized.insert(i,'Mean Normalized Reads',mean_counts)
    #add median
    return joined_normalized

def vennplot_replicates(joined,sample1,sample2,figout):
    Rep1_unique = len(joined[joined['Nuclease_Read_Count.Rep2'] == 0])
    Rep2_unique =len(joined[joined['Nuclease_Read_Count.Rep1'] == 0])
    shared = len(joined[joined['Nuclease_Read_Count.Rep2'] * joined['Nuclease_Read_Count.Rep1'] > 0])
    ja= calc_jaccard(Rep1_unique,Rep2_unique,shared)
    values = (Rep1_unique, Rep2_unique,shared)
    names = (sample1,sample2)
    plt.figure(figsize=(4, 4))
    v =venn2(subsets = values,set_labels=names,set_colors=(colors['sample1'],colors['sample2']),alpha = 0.5)
    v.get_label_by_id("A").set_fontsize(8)
    v.get_label_by_id("B").set_fontsize(8)
    v.get_label_by_id("A").set_y(0.6)
    v.get_label_by_id("B").set_y(0.6)
    v.get_label_by_id("A").set_x(len(sample1)/100.0-0.4)
    v.get_label_by_id("B").set_x(len(sample1)/100.0)
    plt.annotate("% Replicate Sites Overlap " + str(ja), xy=v.get_label_by_id('010').get_position() +
                                           np.array([0, -0.5]), xytext=(-60, -30), ha='center',
                 textcoords='offset points')

    plt.savefig(figout, bbox_inches='tight')
    #plt.show()
    return ja



def repCombiner(samples,name,analysis_folder,PAM='NGG',read_threshold = 6):
    ## Input
    # samples  = ["CPS1_ABE_rep1", "CPS1_ABE_rep2"]
    # analysis_folder = '/groups/clinical/projects/Assay_Dev/changeseq_be/CSB_03/'
    # name = "NGC-ABE8e-V106W CPS1"  # for labeling
    # read_threshold = 2
    # PAM = 'NGN'
    print("Creating PostProcess Folder")
    output_folder = os.path.join(analysis_folder, "postprocess")
    output_plot_folder = os.path.join(analysis_folder, "postprocess/plots")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(output_plot_folder):
        os.makedirs(output_plot_folder)

    ##  Inputs
    sample_files = []
    for sample in samples:
        file = analysis_folder + 'identified/' + sample + '_identified_matched_annotated.csv'
        sample_files.append(analysis_folder + 'identified/' + sample + '_identified_matched_annotated.csv')
        check_file(file)


    ## raw read counts
    joined_out = os.path.join(analysis_folder + 'identified/',"raw_counts_joined_" ) +  name.replace(" ","_") +  '.csv'
    venn_out_without_normalize = analysis_folder + 'visualization/' + name.replace(" ",
                                                                                   "_") + "_raw_counts_venndiagram.png"
    #normalized
    joined_normalized_out = output_folder + "/processed_counts_joined_" + name.replace(" ", "_") + '.csv'
    scatter_out = output_plot_folder + name.replace(" ","_") + "_replicate_scatterplot.png"
    venn_out = output_plot_folder + name.replace(" ","_") + "_postprocess_venn.png"
    swarm_out = output_plot_folder + name.replace(" ","_") + "_postprocess_swarmplot.png"

    # Join without normalizing
    print("Joining...")
    for sample in samples:
        print(sample)


    joined = join_replicates(sample_files, threshold=read_threshold)
    joined.to_csv(joined_out, index=False)


    print("Writing raw_count vendiagram...")
    print(venn_out_without_normalize)
    sim = vennplot_replicates(joined, samples[0],samples[1], venn_out_without_normalize)
    print(sim)

    joined_normalized = normalize(joined, threshold=read_threshold)
    joined_normalized.to_csv(joined_normalized_out, index=False)

    # plotting
    sim = vennplot_replicates(joined_normalized, samples[0],samples[1], venn_out)
    print(sim)

    x1, x2 = list(joined_normalized['Nuclease_Read_Count.Rep1']), list(joined_normalized['Nuclease_Read_Count.Rep2'])

    scatter_plot(x1, x2, name, scatter_out)
    swarm_plot(joined_normalized, name, swarm_out)

    for i in range(len(samples)):
        offtargets, target_seq, total_seq = make_offtarget_dict(joined_normalized, subset = f'Nuclease_Read_Count.Rep{i+1}')
        alignment_plot = output_plot_folder + samples[i].replace(" ", "_") + "_postprocess_alignment_plot.svg"
        draw_plot(target_seq, offtargets, total_seq, outfile = alignment_plot, title=samples[i], PAM=PAM)

    offtargets, target_seq, total_seq = make_offtarget_dict(joined_normalized, subset='Mean Normalized Reads')
    alignment_plot = output_plot_folder + "sample_means" + "_postprocess_alignment_plot.svg"
    draw_plot(target_seq, offtargets, total_seq, outfile=alignment_plot, title=name, PAM=PAM)

'''
data_similarity = []
read_thresholds = [6,12,24,48,96,192]
normalize_data = True
for r in read_thresholds:
    joined = join_replicates(sample1_identified_file,sample2_identified_file,threshold = r)
    if normalize_data:
        joined_normalized = normalize(joined,threshold=r)
        sim =vennplot_replicates(joined_normalized, sample1_name, sample2_name)
    else:
        sim = vennplot_replicates(joined,sample1,sample2)
    data_similarity.append(sim)
    
df = pd.read_csv("/groups/clinical/projects/Assay_Dev/CHANGEseq/CASAFE/casoffinder/hg38_CasSAFE_casoffinder.txt",sep = "\t")
df = df.drop(columns = ["Guide_ID","crRNA"])
joined_normalized = joined_normalized.rename(columns = {'Genomic Coordinate':'Coordinates'})

df2 = df.join(joined_normalized.set_index('Coordinates'), on = 'Coordinates')
out = os.path.join(analysis_folder, 'identified',
                                         name) + '_JOINED_NORMALIZED_casoffinder.csv'

df2[~df2['DNA'].isna()]   
    
'''


def parse_args():
    mainParser = argparse.ArgumentParser()
    mainParser.add_argument('--sample1', '-s1', help='name of replicate 1. Must match sample in manifest')
    mainParser.add_argument('--sample2', '-s2', help='name of replicate 2. Must match sample in manifest')
    #mainParser.add_argument('--file1', '-f1', help='absolute path of sample1 identified_annotated.csv file')
    #mainParser.add_argument('--file2', '-f2', help= 'absolute path of sample2 identified_annotated.csv file')
    mainParser.add_argument('--output', '-o', help='output directory')
    mainParser.add_argument('--name', '-n', help='Sample name for labeling graphs')
    mainParser.add_argument('--read_threshold', '-rt', help='limit sample combining to a min read count threshold',default =6)

    return mainParser.parse_args()


def main():
    args = parse_args()
    repCombiner(args.sample1, args.sample2,args.name, args.output, int(args.read_threshold))

if __name__ == '__main__':
    main()