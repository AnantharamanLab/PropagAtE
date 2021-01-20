#! /usr/bin/env python3
# Author: Kristopher Kieft, UW-Madison
# kieft@wisc.edu


# PropagAtE: Prophage Activity Estimator
# Version: v1.0.0
# Release date: January 21 2021
# dev: 22

# Imports
try:
    import warnings
    warnings.filterwarnings("ignore")
    import argparse
    import subprocess
    import sys
    import random
    import time
    import datetime
    from datetime import date
    from scipy import stats
    import statistics
    from scipy.stats import mannwhitneyu
    from scipy.stats import combine_pvalues
    import math
    import uuid
    import logging
    import os
    from Bio.SeqIO.FastaIO import SimpleFastaParser
except Exception as e:
    sys.stderr.write("\nError: please verify dependancy imports are installed and up to date:\n\n")
    sys.stderr.write(str(e) + "\n\n")
    exit(1)

# Set up variables
start = time.time()
start_time = str(datetime.datetime.now().time()).rsplit(".",1)[0]
u = str(uuid.uuid1()).split("-")[0]

propagate = argparse.ArgumentParser(description='Example 1: python3 PropagAtE_run.py -f scaffolds.fasta -r forward.fastq reverse.fastq -o output.tsv -v prophage_coordinates.tsv -t threads -i ID  .  .  .  Example 2: python3 PropagAtE_run.py -s alignment.sam -o output.tsv -v prophage_coordinates.tsv')
propagate.add_argument('--version', action='version', version='PropagAtE v1.0.0')

# Input required
propagate.add_argument('-s', type=str, nargs=1, default = [''], help='input SAM (.sam) sequence alignment file (skip -sb/-b/-r/-f).')
propagate.add_argument('-b', type=str, nargs=1, default = [''], help='input BAM (.bam) sequence alignment file (skip -sb/-s/-r/-f).')
propagate.add_argument('-sb', type=str, nargs=1, default = [''], help='input sorted BAM (.bam) sequence alignment file (skip -b/-s/-r/-f).')
propagate.add_argument('-r', type=str, nargs=2, default = ['',''], help='input forward and reverse read files (.fastq) separated by a space (skip -sb/-b/-s, requires -f)')
propagate.add_argument('-v', type=str, nargs=1, required = True, help='VIBRANT "integrated_prophage_coordinates" file or custom prophage coordinates file.')

# Optional
propagate.add_argument('-o', type=str, nargs=1, default = ['PropagAtE_result.' + str(u) + '.tsv'], help='name of output tab-separated file for results.')
propagate.add_argument('-f', type=str, nargs=1, default = [''], help='if input is -r provide genomes/scaffolds (.fna, .fastq, .fa, .fsa).')
propagate.add_argument('-t', type=str, nargs=1, default = ['1'], help='threads used for Bowtie2 mapping (if input is -r and -f).')
propagate.add_argument('-i', type=str, nargs=1, default = [''], help='unique ID to append to Bowtie2 SAM output for reference to read sets used (not requied, useful when analyzing multiple read sets).')
propagate.add_argument('-p', type=str, nargs=1, default = ['0.05'], help='p-value cutoff for Mann-Whitney statistical test (default=0.05).')
propagate.add_argument('-g', type=str, nargs=1, default = ['1'], choices=["0", "1", "2", "3"], help='number of gap extensions allowed in read alignment (per read) [default=1, options=0-3].')
propagate.add_argument('-m', type=str, nargs=1, default = ['3'], choices=["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"], help='number of mismatches allowed in read alignment (per read) [default=3, options=0-10].')
propagate.add_argument('-n', type=str, nargs=1, default = ['5'], help='number of coverage sub-samples to take for Mann-Whitney statistical analysis [default=5, minimum=5].')
propagate.add_argument('-a', type=str, nargs=1, default = ['4'], choices=["0", "1", "2", "3", "4"], help='remove outlier coverage values "a" standard deviations from the average [default=4].')
propagate.add_argument('-e', type=str, nargs=1, default = ['0.75'], help="minimum effect size for significance by Cohen's d test [default=0.75, minimum=0.6].")
propagate.add_argument('-c', type=str, nargs=1, default = ['1.75'], help="minimum prophage:host coverage ratio for significance [default=1.75, minimum=1.5].")
propagate.add_argument('-y', type=str, nargs=1, default = [''], help='term to distinguish host from prophage name in genome/scaffold definition lines [default="_fragment_"].')
propagate.add_argument('-x', type=str, nargs=1, default = [''], help="path to Samtools executable if not in PATH.")
propagate.add_argument('-clean', action='store_true', help='use this setting to remove (if applicable) any generated SAM, unsorted BAM and Bowtie2 index files. Will retain any user input files, raw reads and sorted BAM files (default=off).')
propagate.add_argument('-all', action='store_true', help='use this setting to keep all aligned reads (ignore -g and -m, default=off).')
propagate.add_argument('-z', type=str, nargs=1, default = [''], help='character(s) used to replace spaces in genome/scaffold names. Use this setting if spaces were replaced in the alignment file (SAM/BAM) but not in the prophage coordinates file.')

# Parse arguments
args = propagate.parse_args()
samfile = str(args.s[0])
bamfile = str(args.b[0])
sortbamfile = str(args.sb[0])
fasta = str(args.f[0])
forward = str(args.r[0])
reverse = str(args.r[1])
threads = str(args.t[0])
vibe = str(args.v[0])
id = str(args.i[0])
pval = float(args.p[0])
pval = "%.3f" % pval
gaps = int(args.g[0])
mismatches = int(args.m[0])
subs = int(args.n[0])
outliers = int(args.a[0])
effect = float(args.e[0])
ratio_cutoff = float(args.c[0])
termsplit = str(args.y[0])
outfile = str(args.o[0])
outpath = str(os.path.dirname(os.path.abspath(outfile))) + "/"
try:
    temp_out = str(outfile).rsplit("/",1)[1]
except Exception:
    temp_out = str(outfile)
outfile = outpath + str(temp_out)
sampath = str(args.x[0])
spaces = str(args.z[0])
replace = False
if termsplit == '':
    termsplit = "_fragment_"
if spaces == '':
    spaces == '~~'
else:
    replace = True
if sampath != '':
    if "samtools" == sampath[-8:]:
        sampath = sampath[:-8]
    if sampath[-1] != '/':
        sampath += '/'

# log file
logfilename = str(outfile).rsplit(".",1)[0] + ".log"
subprocess.run("rm " + logfilename + " 2> /dev/null", shell=True)
logging.basicConfig(filename=logfilename, level=logging.INFO, format='%(message)s')

# verify inputs and versions
if termsplit[0] == " ":
    sys.stderr.write("\nThe -y input term cannot start with a space. If you tried '-y -gene' then try '-y-gene' instead. Exiting.\n")
    logging.info("\nThe -y input term cannot start with a space. If you tried '-y -gene' then try '-y-gene' instead. Exiting.\n")
    exit(1)
if outfile.rsplit(".",1)[1] == "log":
    sys.stderr.write("\nThe output file cannot end in .log, suggested is .tsv. Exiting.\n")
    logging.info("\nThe output file cannot end in .log, suggested is .tsv. Exiting.\n")
    exit(1)
if (samfile != '' and bamfile != '') or (samfile != '' and forward != '') or (bamfile != '' and forward != '') or (samfile != '' and sortbamfile != '') or (sortbamfile != '' and bamfile != '') or (forward != '' and sortbamfile != ''):
    sys.stderr.write("\nOnly one input file (-s, -sb, -b, -r) is allowed. Exiting.\n")
    logging.info("\nOnly one input file (-s, -sb, -b, -r) is allowed. Exiting.\n")
    exit(1)
try:
    if (forward.rsplit('.',1)[1] != 'fastq' or reverse.rsplit('.',1)[1] != 'fastq') and ('.'.join(forward.rsplit('.',2)[1:]) != 'fastq.gz' or '.'.join(reverse.rsplit('.',2)[1:]) != 'fastq.gz'):
        sys.stderr.write("\nError: Provided reads files must both have the extension .fastq or .fastq.gz. Exiting.\n")
        logging.info("\nError: Provided reads files must both have the extension .fastq or .fastq.gz. Exiting.\n")
        exit(1)
except Exception:
    pass
try:
    if forward != '' and reverse != '' and fasta.rsplit('.',1)[1] != 'fa' and fasta.rsplit('.',1)[1] != 'fna' and fasta.rsplit('.',1)[1] != 'fasta' and fasta.rsplit('.',1)[1] != 'fsa':
        sys.stderr.write("\nError: Provided fasta file must have the extension .fasta, .fna, .fa or .fsa. Exiting.\n")
        logging.info("\nError: Provided fasta file must have the extension .fasta, .fna, .fa or .fsa. Exiting.\n")
        exit(1)
except Exception:
    pass
if forward != '' and reverse != '' and fasta == '':
    sys.stderr.write("\nError: Reads (-r) were provided but a sequence file (-f) was not. Exiting.\n")
    logging.info("\nError: Reads (-r) were provided but a sequence file (-f) was not. Exiting.\n")
    exit(1)
try:
    if samfile.rsplit('.',1)[1] != 'sam':
        sys.stderr.write("\nError: Provided sam file must have the extension .sam. Exiting.\n")
        logging.info("\nError: Provided sam file must have the extension .sam. Exiting.\n")
        exit(1)
except Exception:
    pass
try:
    if bamfile.rsplit('.',1)[1] != 'bam':
        sys.stderr.write("\nError: Provided bam file must have the extension .bam. Exiting.\n")
        logging.info("\nError: Provided bam file must have the extension .bam. Exiting.\n")
        exit(1)
except Exception:
    pass
try:
    if sortbamfile.rsplit('.',1)[1] != 'bam':
        sys.stderr.write("\nError: Provided bam (sorted) file must have the extension .bam. Exiting.\n")
        logging.info("\nError: Provided bam (sorted) file must have the extension .bam. Exiting.\n")
        exit(1)
except Exception:
    pass

if not os.path.exists(vibe):
    sys.stderr.write("\nError: could not identify prophage coordinates file (-v). Exiting.\n")
    logging.info("\nError: could not identify prophage coordinates file (-v). Exiting.\n")
    exit(1)
if effect < 0.6:
    sys.stderr.write("\nError: Cohen's d effect size should not be set below 0.6. Exiting.\n")
    logging.info("\nError: Cohen's d effect size should not be set below 0.6. Exiting.\n")
    exit(1)
if ratio_cutoff < 1.5:
    sys.stderr.write("\nError: ratio cutoff should not be set below 1.5. Exiting.\n")
    logging.info("\nError: ratio cutoff should not be set below 1.5. Exiting.\n")
    exit(1)
if subs < 5:
    sys.stderr.write("\nError: subsample number should not be set below 5. Exiting.\n")
    logging.info("\nError: subsample number should not be set below 5. Exiting.\n")
    exit(1)
try:
    subprocess.check_output("which " + sampath + "samtools", shell=True)
except Exception:
    sys.stderr.write("\nError: samtools does not appear to be installed, is not in the system's PATH or is not in the indicated path (-x). Exiting.\n")
    logging.info("\nError: samtools does not appear to be installed, is not in the system's PATH or is not in the indicated path (-x). Exiting.\n")
    exit(1)
try:
    subprocess.check_output(sampath + "samtools depth -a placeholder.bam 2> " + outpath + u + "_samtools_checkfile.txt", shell=True)
except Exception as e:
    with open(outpath + u + "_samtools_checkfile.txt",'r') as checkfile:
        for line in checkfile:
            if "depth: invalid option -- 'a'" in line:
                sys.stderr.write("\nError: samtools version is incorrect. Please update and verify samtools '-a' flag is available.")
                sys.stderr.write("\n       Alternatively use the '-x' flag to provide the location of the correct samtools executable.")
                sys.stderr.write("\n       Example update using conda: 'conda install -c bioconda samtools'. Exiting.\n")
                logging.info("\nError: samtools version is incorrect. Please update and verify samtools '-a' flag is available.")
                logging.info("\n       Alternatively use the '-x' flag to provide the location of the correct samtools executable.")
                subprocess.run("rm " + outpath + "samtools_checkfile.txt 2> /dev/null", shell=True)
                exit(1)
subprocess.run("rm " + outpath + u + "_samtools_checkfile.txt 2> /dev/null", shell=True)

if forward == '' and fasta != '':
    logging.info("Note: Reads not provided, ignoring input fasta file.\n")
if forward == '' and threads != '1':
    logging.info("Note: Reads not provided, using 1 thread.\n")
if forward == '' and id != '':
    logging.info("Note: Reads not provided, ignoring identifier ID input.\n")

# Statistics
def MWdownsampling(phage, background):
    n = 1
    pvalue_list = []
    while n <= subs:
        cov_p_100 = random.sample(phage, 100)
        cov_h_100 = random.sample(background, 100)
        if sum(cov_p_100) != 0 and sum(cov_h_100) != 0:
            stat, pvalue = mannwhitneyu(cov_p_100, cov_h_100)
            pvalue_list.append(pvalue)
        n += 1
    result = combine_pvalues(pvalue_list)
    return result

def cohenD(phage_mean, phage_sd, host_mean, host_sd):
    pool = math.sqrt((phage_sd**2+host_sd**2)/2)
    d = (host_mean-phage_mean)/pool
    return abs(d)


##### ----------------------------------------------------------------------------------------------------------------------- #####
logging.info("Command:    %s" % ' '.join(sys.argv))
logging.info("Outfolder:  %s" % outpath)
logging.info("")
logging.info("Date:       %s" % str(date.today()))
logging.info("Time:       %s" % str(datetime.datetime.now().time()).rsplit(".",1)[0])
logging.info("Program:    PropagAtE v1.0.0")
logging.info("Run ID:     %s\n" % u)

logging.info("Time (min) |  Log                                                   ")
logging.info("--------------------------------------------------------------------")

reads_input = False
sam_input = False
bam_input = False
sortbam_input = False

# If input is reads/fasta run Bowtie2
if forward != '' and fasta != '':
    reads_input = True
    logging.info("%s         Reads input identified, using %s threads to run Bowtie2." % (str(round((time.time() - float(start))/60,1)),threads))

    try:
        subprocess.check_output("which bowtie2", shell=True)
    except Exception:
        sys.stderr.write("\nError: Bowtie2 does not appear to be installed or is not in the system's PATH. Exiting.\n")
        logging.info("\nError: Bowtie2 does not appear to be installed or is not in the system's PATH. Exiting.\n")
        exit(1)
    try:
        temp = fasta.rsplit("/",1)[1]
        base = temp.rsplit(".",1)[0]
    except Exception:
        base = fasta.rsplit(".",1)[0]

    logging.info("%s         Checking format of FASTA file, replacing spaces if necessary." % str(round((time.time() - float(start))/60,1)))
    with open(fasta, 'r') as fasta_in, open(outpath + fasta.rsplit(".",1)[0] + ".no-spaces." + fasta.rsplit(".",1)[1], 'w') as fasta_out:
        for name,seq in SimpleFastaParser(fasta_in):
            if " " in name:
                replace = True
                name = name.replace(" ", str(spaces))
            fasta_out.write(">" + str(name) + "\n" + str(seq) + "\n")
    if replace == True:
        logging.info('%s         Spaces identified in sequence names. Replaced with "~~".' % str(round((time.time() - float(start))/60,1)))
        fasta = outpath + fasta.rsplit(".",1)[0] + ".no-spaces." + fasta.rsplit(".",1)[1]

    if not os.path.exists(outpath + str(base)+".bowtie2.index.1.bt2"):
        logging.info("%s         Building index ..." % str(round((time.time() - float(start))/60,1)))
        build = True
        subprocess.run("bowtie2-build " + str(fasta) + " " + outpath + str(base)+".bowtie2.index > /dev/null 2> /dev/null", shell=True)
    else:
        build = False
        logging.info("%s         Existing index identified, skipping build ..." % str(round((time.time() - float(start))/60,1)))
    logging.info("%s         Mapping reads ..." % str(round((time.time() - float(start))/60,1)))
    if id != '':
        subprocess.run("bowtie2 -x " + outpath + str(base)+".bowtie2.index -1 " + str(forward) + " -2 " + str(reverse) + " -S " + outpath + str(base)+"." + str(id) + ".sam -q -p " + str(threads) + " --no-unal > /dev/null 2> /dev/null", shell=True)
    else:
        subprocess.run("bowtie2 -x " + outpath + str(base)+".bowtie2.index -1 " + str(forward) + " -2 " + str(reverse) + " -S " + outpath + str(base)+".sam -q -p " + str(threads) + " --no-unal > /dev/null 2> /dev/null", shell=True)
    if args.clean == True and build == True:
        subprocess.run("rm " + outpath + str(base)+".bowtie2.index.* 2> /dev/null", shell=True)
    logging.info("%s         Bowtie2 done." % str(round((time.time() - float(start))/60,1)))
    if str(id) != "":
        samfile = outpath + str(base)+"." + str(id) + ".sam"
    else:
        samfile = outpath + str(base)+".sam"
    subprocess.run("rm " + outpath + fasta.rsplit(".",1)[0] + ".no-spaces." + fasta.rsplit(".",1)[1] + " 2> /dev/null", shell=True)

# Read in prophage coordinate data
prophage_dict = {}
prophage_dict_frags = {}
vibrant_genomes_full = []
logging.info("%s         Generating a list of all prophage regions ..." % str(round((time.time() - float(start))/60,1)))
logging.info("%s         Distinguishing host from prophage names by the term '%s' ..." % (str(round((time.time() - float(start))/60,1)),str(termsplit)))
with open(vibe, 'r') as vibrant:
    phage_list = vibrant.read().replace("\n","\t").split("\t")
    if phage_list[0] == "scaffold" and phage_list[1] == "fragment" and phage_list[5] == "nucleotide start" and phage_list[6] == "nucleotide stop":
        n = 9
        while n < len(phage_list):
            name = str(phage_list[n]).split(termsplit)[0]
            if replace == True:
                name = name.replace(" ", str(spaces))
            if name[0] == "'" or name[0] == '"':
                name = name[1:]
            if replace == True:
                frag = str(phage_list[n]).replace(" ", str(spaces))
            else:
                frag = str(phage_list[n])
            if frag[0] == "'" or frag[0] == '"':
                frag = frag[1:-1]

            if n == 9: # in case of error
                save_name = name
                save_frag = frag
            if name != '' and frag != '' and name == frag: # names were unable to be split
                sys.stderr.write("\nError: prophage names could not be distinguished from host names. Check accuracy of -y flag input. Below are the names that caused the issue. Exiting.\n")
                sys.stderr.write("Host name detected:       %s\n" % name)
                sys.stderr.write("Prophage name detected:   %s\n" % frag)
                sys.stderr.write("-y flag input detected:   %s\n" % termsplit)
                logging.info("\nError: prophage names could not be distinguished from host names. Check accuracy of -y flag input. Below are the names that caused the issue. Exiting.\n")
                logging.info("Host name detected:       %s" % name)
                logging.info("Prophage name detected:   %s" % frag)
                logging.info("-y flag input detected:   %s" % termsplit)
                exit(1)

            prophage_dict.update({str(name) + "~" + str(phage_list[n+4]):str(frag)})
            prophage_dict.update({str(name) + "~" + str(phage_list[n+5]):str(frag)})
            prophage_dict_frags.update({str(name):str(frag)})
            vibrant_genomes_full.append(str(name))

            n += 8

    elif phage_list[0] == "scaffold" and phage_list[1] == "fragment" and phage_list[2] == "start" and phage_list[3] == "stop":
        n = 5
        while n < len(phage_list):
            name = str(phage_list[n]).split(termsplit)[0]
            if replace == True:
                name = name.replace(" ", str(spaces))
            if name[0] == "'" or name[0] == '"':
                name = name[1:]
            if replace == True:
                frag = str(phage_list[n]).replace(" ", str(spaces))
            else:
                frag = str(phage_list[n])
            if frag[0] == "'" or frag[0] == '"':
                frag = frag[1:-1]

            if n == 5: # in case of error
                save_name = name
                save_frag = frag
            if name != '' and frag != '' and name == frag: # names were unable to be split
                sys.stderr.write("\nError: prophage names could not be distinguished from host names. Check accuracy of -y flag input. Below are the names that caused the issue. Exiting.\n")
                sys.stderr.write("Host name detected:       %s\n" % name)
                sys.stderr.write("Prophage name detected:   %s\n" % frag)
                sys.stderr.write("-y flag input detected:   %s\n" % termsplit)
                logging.info("\nError: prophage names could not be distinguished from host names. Check accuracy of -y flag input. Below are the names that caused the issue. Exiting.\n")
                logging.info("Host name detected:       %s" % name)
                logging.info("Prophage name detected:   %s" % frag)
                logging.info("-y flag input detected:   %s" % termsplit)
                exit(1)

            prophage_dict.update({str(name) + "~" + str(phage_list[n+1]):str(frag)})
            prophage_dict.update({str(name) + "~" + str(phage_list[n+2]):str(frag)})
            prophage_dict_frags.update({str(name):str(frag)})
            vibrant_genomes_full.append(str(name))

            n += 4

    else:
        sys.stderr.write("\nError: -v coordinates file appears to be incorrect. Either input a VIBRANT 'integrated_prophage_coordinates' file or a custom tab separated coordinates file with the header 'scaffold/fragment/start/stop'. Exiting.\n")
        logging.info("\nError: -v coordinates file appears to be incorrect. Either input a VIBRANT 'integrated_prophage_coordinates' file or a custom tab separated coordinates file with the header 'scaffold/fragment/start/stop'. Exiting.\n")
        exit(1)

    number_hosts = list(set(vibrant_genomes_full))
    logging.info("%s         Number of prophage regions identified: %s" % (str(round((time.time() - float(start))/60,1)),len(vibrant_genomes_full)))
    logging.info("%s         Number of unique host regions identified: %s" % (str(round((time.time() - float(start))/60,1)),len(number_hosts)))
    vibrant_genomes_full = list(set(vibrant_genomes_full))

# process sam/bam files
logging.info("%s         Processing alignment (.sam/.bam) file." % str(round((time.time() - float(start))/60,1)))

# Gap/mismatch processing for sam file, convert sam to bam
processed = False
if samfile != '':
    sam_input = True
    try:
        temp = samfile.rsplit("/",1)[1]
        base = temp.rsplit(".",1)[0]
    except Exception:
        base = samfile.rsplit(".",1)[0]

    if args.all == False:
        processed = True
        logging.info("%s         Alignment .sam input identified, processing gaps/mismatches ..." % str(round((time.time() - float(start))/60,1)))
        tossed = 0
        kept = 0
        with open(samfile, 'r') as sf, open(outpath + str(base)+'.gap-' + str(gaps) + '_mm-' + str(mismatches) + '.sam', 'w') as psf:
            for line in sf:
                if "XG:i:" in line and "XM:i:" in line:
                    xg = line.split("\t")[15]
                    xg = xg.split(":")[2]
                    xm = line.split("\t")[13]
                    xm = xm.split(":")[2]
                    if int(xg) <= gaps and int(xm) <= mismatches:
                        psf.write(line)
                        kept += 1
                    else:
                        tossed += 1
                else:
                    psf.write(line)
        if reads_input == True and args.clean == True:
            subprocess.run("rm " + str(samfile), shell=True) # remove non-processed sam
        samfile = outpath + str(base)+'.gap-' + str(gaps) + '_mm-' + str(mismatches) + '.sam'
        try:
            temp = samfile.rsplit("/",1)[1]
            base = temp.rsplit(".",1)[0]
        except Exception:
            base = samfile.rsplit(".",1)[0]
        logging.info("%s         %s of %s total aligned reads were tossed for having too many gaps/mismatches ..." % (str(round((time.time() - float(start))/60,1)),str(tossed),str(tossed+kept)))
    else:
        logging.info("%s         Alignment .sam input identified, skipping gaps/mismatches processing ..." % str(round((time.time() - float(start))/60,1)))

    logging.info("%s         Converting .sam to .bam ..." % str(round((time.time() - float(start))/60,1)))
    subprocess.run(sampath + "samtools view -S -b " + str(samfile) + " > " + outpath + str(base) + ".bam 2> /dev/null", shell=True)
    logging.info("%s         Conversion done." % str(round((time.time() - float(start))/60,1)))
    bamfile = outpath + str(base) + ".bam"

samfile_cleaned = False
if reads_input == True and args.clean == True:
    samfile_cleaned = True
    subprocess.run("rm " + str(samfile), shell=True) # remove processed or non-processed sam

# Gap/mismatch processing for bam file
if bamfile != '':
    bam_input = True
    if sam_input == False and args.all == False:
        try:
            temp = bamfile.rsplit("/",1)[1]
            base = temp.rsplit(".",1)[0]
        except Exception:
            base = bamfile.rsplit(".",1)[0]
        logging.info("%s         Alignment .bam input identified, coverting to .sam for gap/mismatch filtering ..." % str(round((time.time() - float(start))/60,1)))
        subprocess.run(sampath + "samtools view -h -o " + outpath + str(base) + ".sam " + str(bamfile) + " 2> /dev/null", shell=True)

        samfile = outpath + str(base) + ".sam"
        processed = True
        logging.info("%s         Alignment .sam input identified, processing gaps/mismatches ..." % str(round((time.time() - float(start))/60,1)))
        tossed = 0
        kept = 0
        with open(samfile, 'r') as sf, open(outpath + str(base)+'.gap-' + str(gaps) + '_mm-' + str(mismatches) + '.sam', 'w') as psf:
            for line in sf:
                if "XG:i:" in line and "XM:i:" in line:
                    xg = line.split("\t")[15]
                    xg = xg.split(":")[2]
                    xm = line.split("\t")[13]
                    xm = xm.split(":")[2]
                    if int(xg) <= gaps and int(xm) <= mismatches:
                        psf.write(line)
                        kept += 1
                    else:
                        tossed += 1
                else:
                    psf.write(line)
        if args.clean == True:
            subprocess.run("rm " + str(samfile), shell=True) # remove non-processed sam
        samfile = outpath + str(base)+'.gap-' + str(gaps) + '_mm-' + str(mismatches) + '.sam'
        try:
            temp = samfile.rsplit("/",1)[1]
            base = temp.rsplit(".",1)[0]
        except Exception:
            base = samfile.rsplit(".",1)[0]
        logging.info("%s         %s of %s total aligned reads were tossed for having too many gaps/mismatches ..." % (str(round((time.time() - float(start))/60,1)),str(tossed),str(tossed+kept)))
        logging.info("%s         Converting processed .sam to .bam ..." % str(round((time.time() - float(start))/60,1)))
        subprocess.run(sampath + "samtools view -S -b " + str(samfile) + " > " + outpath + str(base) + ".bam 2> /dev/null", shell=True)
        logging.info("%s         Conversion done." % str(round((time.time() - float(start))/60,1)))
        bamfile = outpath + str(base) + ".bam"

        if args.clean == True:
            samfile_cleaned = True
            subprocess.run("rm " + str(samfile), shell=True) # remove processed or non-processed sam

    if sam_input == False and args.all == True:
        logging.info("%s         Alignment .bam input identified, sorting .bam file ..." % str(round((time.time() - float(start))/60,1)))
    try:
        temp = bamfile.rsplit("/",1)[1]
        base = temp.rsplit(".",1)[0]
    except Exception:
        base = bamfile.rsplit(".",1)[0]
    subprocess.run(sampath + "samtools sort -o " + outpath + str(base) + ".sorted.bam " + str(bamfile) + " 2> /dev/null", shell=True)
    logging.info("%s         Sorting done." % str(round((time.time() - float(start))/60,1)))
    sortbam = outpath + str(base) + ".sorted.bam"

bamfile_cleaned = False
if sam_input == True and args.clean == True:
    bamfile_cleaned = True
    subprocess.run("rm " + str(bamfile), shell=True)

# Gap/mismatch processing for sorted bam file
if sortbamfile != '':
    if args.all == False and processed == False:
        try:
            temp = sortbamfile.rsplit("/",1)[1]
            base = temp.rsplit(".",1)[0]
        except Exception:
            base = sortbamfile.rsplit(".",1)[0]

        logging.info("%s         Sorted .bam input identified, coverting to .sam for gap/mismatch filtering ..." % str(round((time.time() - float(start))/60,1)))
        subprocess.run(sampath + "samtools view -h -o " + outpath + str(base) + ".sam " + str(sortbamfile) + " 2> /dev/null", shell=True)
        samfile = outpath + str(base) + ".sam"
        processed = True
        logging.info("%s         Processing gaps/mismatches ..." % str(round((time.time() - float(start))/60,1)))
        tossed = 0
        kept = 0
        with open(samfile, 'r') as sf, open(outpath + str(base)+'.gap-' + str(gaps) + '_mm-' + str(mismatches) + '.sam', 'w') as psf:
            for line in sf:
                if "XG:i:" in line and "XM:i:" in line:
                    xg = line.split("\t")[15]
                    xg = xg.split(":")[2]
                    xm = line.split("\t")[13]
                    xm = xm.split(":")[2]
                    if int(xg) <= gaps and int(xm) <= mismatches:
                        psf.write(line)
                        kept += 1
                    else:
                        tossed += 1
                else:
                    psf.write(line)

        if args.clean == True:
            subprocess.run("rm " + str(samfile), shell=True) # remove non-processed sam
        samfile = outpath + str(base)+'.gap-' + str(gaps) + '_mm-' + str(mismatches) + '.sam'
        try:
            temp = samfile.rsplit("/",1)[1]
            base = temp.rsplit(".",1)[0]
        except Exception:
            base = samfile.rsplit(".",1)[0]
        logging.info("%s         %s of %s total aligned reads were tossed for having too many gaps/mismatches ..." % (str(round((time.time() - float(start))/60,1)),str(tossed),str(tossed+kept)))

        logging.info("%s         Converting processed .sam to .bam ..." % str(round((time.time() - float(start))/60,1)))
        subprocess.run(sampath + "samtools view -S -b " + str(samfile) + " > " + outpath + str(base) + ".bam 2> /dev/null", shell=True)
        logging.info("%s         Conversion done." % str(round((time.time() - float(start))/60,1)))
        bamfile = outpath + str(base) + ".bam"

        if args.clean == True:
            samfile_cleaned = True
            subprocess.run("rm " + str(samfile), shell=True) # remove processed or non-processed sam

        logging.info("%s         Sorting .bam file ..." % str(round((time.time() - float(start))/60,1)))
        try:
            temp = bamfile.rsplit("/",1)[1]
            base = temp.rsplit(".",1)[0]
        except Exception:
            base = bamfile.rsplit(".",1)[0]
        subprocess.run(sampath + "samtools sort -o " + outpath + str(base) + ".sorted.bam " + str(bamfile) + " 2> /dev/null", shell=True)
        logging.info("%s         Sorting done." % str(round((time.time() - float(start))/60,1)))
        sortbam = outpath + str(base) + ".sorted.bam"
        if args.clean == True:
            samfile_cleaned = True
            subprocess.run("rm " + str(bamfile), shell=True) # remove processed or non-processed sam

    sortbam = sortbamfile
    if args.all == True:
        try:
            temp = sortbam.rsplit("/",1)[1]
            base = temp.rsplit(".",1)[0]
        except Exception:
            base = sortbam.rsplit(".",1)[0]
        logging.info("%s         Sorted alignment .bam input identified, skipping sorting and gap/mismatch processing." % str(round((time.time() - float(start))/60,1)))

# Extract coverage for each nucleotide
temp = outpath + str(base) + ".temp-list.txt"
subprocess.run("echo '" + str(sortbam) + "' > " + str(temp), shell=True)
cov_file = str(sortbam).rsplit(".",1)[0] + ".coverage-by-nt.tsv"
logging.info("%s         Indexing sorted .bam file ..." % str(round((time.time() - float(start))/60,1)))
subprocess.run(sampath + "samtools index " + str(sortbam), shell=True) # index bam file
logging.info("%s         Extracting coverage information from sorted .bam file ..." % str(round((time.time() - float(start))/60,1)))

# take the whole scaffold
for genome in vibrant_genomes_full:
    subprocess.run(sampath + "samtools depth -r " + str(genome) + ":0-999999999 -a -f " + str(temp) + " >> " + cov_file + " 2> /dev/null", shell=True) # extract coverage info

subprocess.run("rm " + str(temp) + " 2> /dev/null", shell=True)
subprocess.run("rm " + str(sortbam) + ".bai 2> /dev/null", shell=True) # remove index
subprocess.run("echo 'next\t0\t0\n' >> " + str(cov_file), shell=True)
logging.info("%s         Coverage extraction done." % str(round((time.time() - float(start))/60,1)))

if int(os.stat(cov_file).st_size) < 50: # Exit if file is empty (< 50 bytes)
    if samfile == "":
        samfile = "N/A"
    if bamfile == "":
        bamfile = "N/A"
    subprocess.run("rm " + str(cov_file), shell=True)
    logging.info("%s         No coverage identified for scaffolds with integrated prophages. Analysis finished." % str(round((time.time() - float(start))/60,1)))
    logging.info("")
    logging.info("")
    logging.info("CAUTION:")
    logging.info("         'No coverage identified' may be due to an incorrect -y flag.")
    logging.info("          Please verify -y input correctly distinguishes host from prophage sequences.")
    logging.info("          This can be verified by checking the -v or -f inputs.")
    logging.info("          Example host name detected:       %s" % save_name)
    logging.info("          Example prophage name detected:   %s" % save_frag)
    logging.info("          -y flag input:                    %s" % termsplit)
    logging.info("")
    logging.info("")
    logging.info("Log file:               %s" % logfilename)
    logging.info("Output results file:    N/A")
    if samfile_cleaned == False:
        logging.info("sam file:               %s" % samfile)
    else:
        logging.info("sam file:               N/A")
    if bamfile_cleaned == False:
        logging.info("bam file:               %s" % bamfile)
    else:
        logging.info("bam file:               N/A")
    logging.info("sorted bam file:        %s" % sortbam)
    logging.info("")
    logging.info("")
    logging.info("Number of active prophages identified: 0")
    logging.info("")
    logging.info('                                                               ##')
    logging.info('                                                             ##  ##')
    logging.info('                                                           ##      ##')
    logging.info('######   ##  ##     ##     #######   ######    #####       ##      ##')
    logging.info('##  ##   ##  ##   ##  ##   ##        ##       ##             ##  ##')
    logging.info('######   ######   ######   ##  ###   ######    ###             ##')
    logging.info('##       ##  ##   ##  ##   ##   ##   ##           ##           ##')
    logging.info('##       ##  ##   ##  ##   #######   ######   #####            ##')
    logging.info('                                                            #  ##  #')
    logging.info('                                                           # # ## # #')
    logging.info('                                                          #   #  #   #')
    logging.info('                                                         #            #')
    logging.info("")
    exit(0)

# calculate average read length
total = subprocess.check_output(sampath + "samtools view -c " + str(sortbam), shell=True) # if there are less than 10k reads, use total reads to length check
total = int(float(str(total.strip()).split("'")[1]))
if int(total) < 10000:
    head = total
else:
    head = 10000
length = subprocess.check_output(sampath + "samtools view " + str(sortbam) + " | head -n " + str(head) + " | cut -f 10 | wc -c", shell=True)
length = int(float(str(length.strip()).split("'")[1])/int(head))
logging.info("%s         Number of valid reads used: %s." % (str(round((time.time() - float(start))/60,1)),str(total)))
logging.info("%s         Estimated read length to trim by: %s." % (str(round((time.time() - float(start))/60,1)),str(length)))

# grab only relevant genomes and count genome lengths
prophage_coverage = {}
host_coverage = {}
temp_host = []
temp_phage = []
with open(cov_file, 'r') as covfile:
    logging.info("%s         Identifying coverages for host/prophage regions ..." % str(round((time.time() - float(start))/60,1)))
    prev_genome = ''
    n = 0
    host = False
    phage = False
    host_start = False
    phage_start = False
    host_name = False
    phage_name = False
    for line in covfile:
        if str(line.split("\t")[0]) != str(prev_genome):
            # scaffold is done. Update coverage
            if host == True:
                if len(temp_host) > 0:
                    host_coverage.update({str(prev_genome):temp_host[:-length]}) # trim the last N nucleotides from the last section added
                if len(temp_phage) > 0:
                    prophage_coverage.update({prophage_dict_frags[str(prev_genome)]:temp_phage})
            elif phage == True:
                if len(temp_phage) > 0:
                    prophage_coverage.update({prophage_dict_frags[str(prev_genome)]:temp_phage[:-length]}) # trim the last N nucleotides from the last section added
                if len(temp_host) > 0:
                    host_coverage.update({str(prev_genome):temp_host})
            n = 0
            temp_phage = []
            temp_host = []
            host = False
            phage = False
            host_name = False
            phage_name = False
        prev_genome = str(line.split("\t")[0])

        n += 1 # needs to come after previous statement of writing and length counting

        if (line.split("\t")[0] in vibrant_genomes_full or line.split("\t")[0] == "next"): # trim by first length, only include phages, "next" is for including last scaffold
            if str(line.split("\t")[0]) + "~" + str(line.split("\t")[1]) not in prophage_dict.keys() and n == 1: # is the first nucleotide host?
                if host_name == False:
                    temp_host.append(str(line.split("\t")[0]))
                    host_name = True
                temp_host.append(str(line.split("\t")[2]).strip("\n"))
                host = True
                phage = False
            elif str(line.split("\t")[0]) + "~" + str(line.split("\t")[1]) in prophage_dict.keys() and n == 1:
                if phage_name == False:
                    temp_phage.append(prophage_dict[str(line.split("\t")[0]) + "~" + str(line.split("\t")[1])])
                    phage_name = True
                temp_phage.append(str(line.split("\t")[2]).strip("\n").strip("\n"))
                phage = True
                host = False

            if str(line.split("\t")[0]) + "~" + str(line.split("\t")[1]) not in prophage_dict.keys() and n > 1:
                if phage == True:
                    if phage_name == False:
                        temp_phage.append(prophage_dict[str(line.split("\t")[0]) + "~" + str(line.split("\t")[1])])
                        phage_name = True
                    temp_phage.append(str(line.split("\t")[2]).strip("\n"))
                elif host == True:
                    if host_name == False:
                        temp_host.append(str(line.split("\t")[0]))
                        host_name = True
                    temp_host.append(str(line.split("\t")[2]).strip("\n"))

            if str(line.split("\t")[0]) + "~" + str(line.split("\t")[1]) in prophage_dict.keys() and n > 1: # is the nucleotide phage?
                if host == True: # must mean that host is finished and phage starting
                    host = False
                    phage = True
                    if phage_name == False:
                        temp_phage.append(prophage_dict[str(line.split("\t")[0]) + "~" + str(line.split("\t")[1])])
                        phage_name = True
                    temp_phage.append(str(line.split("\t")[2]).strip("\n"))

                elif phage == True: # phage must be done
                    if phage_name == False:
                        temp_phage.append(prophage_dict[str(line.split("\t")[0]) + "~" + str(line.split("\t")[1])])
                        phage_name = True
                    temp_phage.append(str(line.split("\t")[2]).strip("\n"))
                    prophage_coverage.update({prophage_dict[str(line.split("\t")[0]) + "~" + str(line.split("\t")[1])]:temp_phage}) # write to phage coverage, reset phage
                    temp_phage = []
                    host = True
                    phage = False
                    phage_name = False

# Statistics and metrics
logging.info("%s         Coverage identification done." % str(round((time.time() - float(start))/60,1)))
logging.info("%s         Performing statistical tests ..." % str(round((time.time() - float(start))/60,1)))
total = []
total_adj = []
with open(outfile, 'w') as output:
    output.write("prophage\thost\tactive\tMW_pvalue\tCohenD\tprophage-host_ratio\tmean_difference\tprophage_len\tprophage_mean_cov\tprophage_median_cov\tprophage_sd_cov\thost_len\thost_mean_cov\thost_median_cov\thost_sd_cov\n")
    for key in prophage_coverage.keys():
        host = str(key.split(termsplit)[0])
        if host in host_coverage.keys():
            cov_p = [int(i) for i in prophage_coverage[key][1:]]
            cov_h = [int(i) for i in host_coverage[host][1:]]
            avg_p = statistics.mean(cov_p)
            sd_p = statistics.stdev(cov_p)
            avg_h = statistics.mean(cov_h)
            sd_h = statistics.stdev(cov_h)
            if outliers > 0:
                cov_p = [x for x in cov_p if x >= (avg_p-outliers*sd_p) and x <= (avg_p+outliers*sd_p)]
                cov_h = [x for x in cov_h if x >= (avg_h-outliers*sd_h) and x <= (avg_h+outliers*sd_h)]
                #
                avg_p = statistics.mean(cov_p)
                sd_p = statistics.stdev(cov_p)
                avg_h = statistics.mean(cov_h)
                sd_h = statistics.stdev(cov_h)
            med_h = statistics.median(cov_h)
            med_p = statistics.median(cov_p)
            mean_diff = avg_p - avg_h
            result = MWdownsampling(cov_p, cov_h)
            if mean_diff > 0:
                cd = cohenD(avg_p, sd_p, avg_h, sd_h)
            else:
                cd = 0
            host_zero_ratio = False
            if avg_h > 0:
                ratio = avg_p/avg_h
            else:
                ratio = avg_p
                if avg_p >= 0.5:
                    host_zero_ratio = True
            #
            active = "no"
            if (float(result[1]) <= float(pval) or str(result[1]) == 'nan') and cd >= effect and mean_diff > 0 and (ratio >= ratio_cutoff or host_zero_ratio == True):
                active = "yes"
                total.append("yes")
            if len(cov_h) < 1000 or len(cov_p) < 1000: # minimum host/phage length is 1000 bp
                active = 'NA'
                active_adj = 'NA'
            if replace == True:
                output.write(str(key).replace(str(spaces), " ")+"\t"+str(host).replace(str(spaces), " ")+"\t"+str(active)+"\t"+str(result[1])+"\t"+str(cd)+"\t"+str(ratio)+"\t"+str(mean_diff)+"\t"+str(len(cov_p))+"\t"+str(avg_p)+"\t"+str(med_p)+"\t"+str(sd_p)+"\t"+str(len(cov_h))+"\t"+str(avg_h)+"\t"+str(med_h)+"\t"+str(sd_h)+"\n")
            else:
                output.write(str(key)+"\t"+str(host)+"\t"+str(active)+"\t"+str(result[1])+"\t"+str(cd)+"\t"+str(ratio)+"\t"+str(mean_diff)+"\t"+str(len(cov_p))+"\t"+str(avg_p)+"\t"+str(med_p)+"\t"+str(sd_p)+"\t"+str(len(cov_h))+"\t"+str(avg_h)+"\t"+str(med_h)+"\t"+str(sd_h)+"\n")

        else:
            if replace == True:
                output.write(str(key).replace(str(spaces), " ")+"\t"+str(host).replace(str(spaces), " ")+"\tno_coverage\tna\tna\tna\tna\tna\tna\tna\tna\tna\tna\tna\tna\n")
            else:
                output.write(str(key)+"\t"+str(host)+"\tno_coverage\tna\tna\tna\tna\tna\tna\tna\tna\tna\tna\tna\tna\n")

if samfile == "":
    samfile = "N/A"
if bamfile == "":
    bamfile = "N/A"

subprocess.run("rm " + str(cov_file), shell=True)
logging.info("%s         Statistical tests done. Analysis finished." % str(round((time.time() - float(start))/60,1)))
logging.info("")
logging.info("")
logging.info("Log file:               %s" % logfilename)
logging.info("Output results file:    %s" % outfile)
if samfile_cleaned == False:
    logging.info("sam file:               %s" % samfile)
else:
    logging.info("sam file:               N/A")
if bamfile_cleaned == False:
    logging.info("bam file:               %s" % bamfile)
else:
    logging.info("bam file:               N/A")
logging.info("sorted bam file:        %s" % sortbam)
logging.info("")
logging.info("")
logging.info("Number of active prophages identified: %s" % len(total))
logging.info("")
logging.info('                                                               ##')
logging.info('                                                             ##  ##')
logging.info('                                                           ##      ##')
logging.info('######   ##  ##     ##     #######   ######    #####       ##      ##')
logging.info('##  ##   ##  ##   ##  ##   ##        ##       ##             ##  ##')
logging.info('######   ######   ######   ##  ###   ######    ###             ##')
logging.info('##       ##  ##   ##  ##   ##   ##   ##           ##           ##')
logging.info('##       ##  ##   ##  ##   #######   ######   #####            ##')
logging.info('                                                            #  ##  #')
logging.info('                                                           # # ## # #')
logging.info('                                                          #   #  #   #')
logging.info('                                                         #            #')
logging.info("")


#
#
#
