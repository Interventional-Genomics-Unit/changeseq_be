import os
import logging


logger = logging.getLogger('root')
logger.propagate = False


def write_qc(qc_file,preprocessed_logfile,coverage_stat_file):
    #for sample in c.parameters['samples']:
    # preprocessed_logfile = os.path.join(c.parameters["analysis_folder"], 'preprocessed', sample + '.txt')
    # coverage_stat_file  = os.path.join(c.parameters["analysis_folder"], 'qc', sample + '_aligned_stats.txt')
    #
    # qc_file = os.path.join(c.parameters["analysis_folder"], 'qc', sample + '_qc_report.txt')
    metrics = ['total_reads','total_reads_tn5_filtered',
               'total_reads_aligned',
                'total_reads_mapped_outie (derived from circles)',
               'total_reads_mapped_inner (derived from linear)',
               'mapped to different chromosome (concatenating)',
               'average read length']
    metrics_dict =dict()
    metrics_dict.fromkeys(metrics)
    if os.path.isfile(preprocessed_logfile):
        lines = open(preprocessed_logfile,'r').readlines()

        for i,line in enumerate(lines[0:100]):
            if 'Total read pairs processed:' in line:
                total = int(line.split(':')[-1].strip().replace(",", ""))*2
                metrics_dict['total_reads'] = "{:,}".format(total)

            if 'Read 1 with adapter:' in line:
                r1_n,r1_p = line.split(':')[-1].strip().split(" ")
                r2_n,r2_p = lines[i+1].split(':')[-1].strip().split(" ")
                n = int(r1_n.replace(",", "")) +  int(r2_n.replace(",", ""))
                metrics_dict['total_reads_tn5_filtered']= "{:,}".format(n) + " (" +  str(round((n/total)*100,1)) + "%, of total reads)"

    if os.path.isfile(coverage_stat_file):
        lines = open(coverage_stat_file,'r').readlines()

        for i,line in enumerate(lines[0:100]):
            if 'raw total sequences:' in line:
                total = int(line.split(':')[-1].strip().replace(",", ""))
                metrics_dict['total_reads']= "{:,}".format(total)
            if 'reads mapped:' in line:
                n_aligned = int(line.split(':')[-1].strip())
                metrics_dict[ 'total_reads_aligned'] = "{:,}".format(n_aligned) + " (" +str(round((n_aligned/total)*100,1)) + "%, of total reads)"
            if 'outward oriented pairs' in line:
                out_aligned = int(line.split(':')[-1].strip())
                metrics_dict['total_reads_mapped_outie (derived from circles)'] = "{:,}".format(out_aligned) + " (" + str(
                    round((out_aligned / n_aligned) * 100, 1)) + "%, of aligned reads)"
            if 'inward oriented pairs' in line:
                in_aligned = int(line.split(':')[-1].strip())
                metrics_dict['total_reads_mapped_inner (derived from linear)'] = "{:,}".format(
                    in_aligned) + " (" + str(
                    round((in_aligned / n_aligned) * 100, 1)) + "%, of aligned reads))"

    with open(qc_file, "w") as out:
        for k,v in metrics_dict.items():
            out.write(f"{k}: {v}\n")

