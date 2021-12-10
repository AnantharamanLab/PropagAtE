#! /usr/bin/env python3
# Author: Kristopher Kieft, UW-Madison
# kieft@wisc.edu

# PropagAtE: Prophage Activity Estimator
# Version: v1.1.0
# Release date: January 2022

import subprocess
import os
import sys
import pysam
import numpy as np
from numba import jit


def prophages_check(vibe):
    with open(vibe, 'r') as phage_list:
        header = phage_list.readline().strip('\n').split('\t')
        check_line = phage_list.readline().strip('\n')
        if not check_line:
            return 'no prophages'
        
        try:
            if header[0] == "scaffold" and header[1] == "fragment" and header[5] == "nucleotide start" and header[6] == "nucleotide stop":
                return 'vibrant header'
        except IndexError:
            pass
        try:
            if header[0] == "scaffold" and (header[1] == "fragment" or header[1] == "prophage") and header[2] == "start" and header[3] == "stop":
                return 'custom header'
        except IndexError:
            pass

        sys.stderr.write("\nError: -v coordinates file is formatted incorrectly. See README for details. Exiting.\n")
        return ''

def does_exist(file, descript):
    if os.path.exists(file):
        sys.stderr.write(f"\nError: the {descript} already exists. Exiting.\n\n")
        return True
    return False

def not_exist(file, descript):
    if not os.path.exists(file):
        sys.stderr.write(f"\nError: the {descript} file does not exist. Exiting.\n\n")
        return True
    return False

@jit(nopython=True)
def quick_stats(depth, start, stop, min_cov, all_phages, l, host_depth):
    host_depth[start:stop] = np.nan # do not count any prophages in host coverage
    all_phages += l 
    if l >= 1000: # minimum length allowed
        cov = depth[start:stop]
        a = np.nanmean(cov)
        m = np.nanmedian(cov)
        s = np.nanstd(cov)
        d = np.count_nonzero(cov >= min_cov) / l
    
    return a,m,s,d,all_phages,host_depth

@jit(nopython=True)
def host_stats(depth):
    eff_host = depth.size
    if eff_host > 0:
        a = np.nanmean(depth)
        m = np.nanmedian(depth)
        s = np.nanstd(depth)
        return a,m,s,eff_host
    else:
        return 0,0,0,0

def coverage_stats(depth, genome, prophage_dict, prophage_lengths, min_cov, length):
    '''
    Calculate average, median and standard deviation
    '''
    prophages = prophage_dict[genome]
    phage_covs = {}
    all_phages = 0 # length of all phages on this scaffold
    host_depth = depth.copy() # in case of prophage overlap
    for p in prophages:
        phage,start,stop = p
        l = prophage_lengths[phage]
        a,m,s,d,all_phages,host_depth = quick_stats(depth, start, stop, min_cov, all_phages, l, host_depth)
        phage_covs[phage] = (a,m,s,l,d)

    host_depth_filter = host_depth[~np.isnan(host_depth)]
    a,m,s,eff_host = host_stats(host_depth_filter) # host
    return phage_covs, a, m, s, eff_host

@jit(nopython=True)
def add_depth(depth, start, end, ed, rl, read_id):
    '''
    Add depth with read alignment filter
    '''
    if ed/rl <= read_id:
        for i in range(start,end):
            depth[i] += 1
    return depth

@jit(nopython=True)
def add_depth_no_ed(depth, start, end):
    '''
    Add depth, ignore read alignment identity
    '''
    for i in range(start,end):
        depth[i] += 1
    return depth

def fasta_parse(infile):
    '''
    Standard fasta parser.
    '''
    with open(infile, 'r') as fasta:
        for line in fasta:
            if line[0] == '>':
                try:
                    yield header, seq
                except NameError:
                    pass # first line
                seq = ''
                header = line[1:].strip("\n")
            else:
                seq += line.strip("\n")

        # last one
        yield header, seq

def bowtie2_build(fasta, base, outpath, spaces, u):
    '''
    Build bowtie2 index.
    '''
    build = False
    if not os.path.exists(fasta.rsplit('.',1)[0] + ".bowtie2.index.1.bt2"):
        if spaces:
            with open(f"{outpath}{base}.{u}.fasta", 'w') as fasta_out:
                for name,seq in fasta_parse(fasta):
                    name = name.replace(" ", "~!!~")
                    fasta_out.write(f">{name}\n{seq}\n")
            fasta = f"{outpath}{base}.{u}.fasta"

        subprocess.run(f"bowtie2-build {fasta} {outpath}{base}.bowtie2.index > /dev/null 2> /dev/null", shell=True)
        build = True
    
    return build

def run_bowtie2_paired(base, outpath, threads, forward, reverse):
    subprocess.run(f"bowtie2 -x {outpath}{base}.bowtie2.index -1 {forward} -2 {reverse} -S {outpath}{base}.sam -q -p {threads} --no-unal --no-discordant > /dev/null 2> /dev/null", shell=True)

def run_bowtie2_interleaved(base, outpath, threads, interleaved):
    subprocess.run(f"bowtie2 -x {outpath}{base}.bowtie2.index --interleaved {interleaved} -S {outpath}{base}.sam -q -p {threads} --no-unal --no-discordant > /dev/null 2> /dev/null", shell=True)

def run_bowtie2_unpaired(base, outpath, threads, unpaired):
    subprocess.run(f"bowtie2 -x {outpath}{base}.bowtie2.index -U {unpaired} -S {outpath}{base}.sam -q -p {threads} --no-unal > /dev/null 2> /dev/null", shell=True)

def post_bowtie2(outpath, base, clean, u, build, threads):
    '''
    Clean up bowtie2 and SAM -> BAM
    '''
    sam = f'{outpath}{base}.sam'
    bam = f'{outpath}{base}.bam'
    subprocess.run(f"samtools view -@ {threads} -h -b {sam} > {bam} 2> /dev/null", shell=True)
    
    if clean:
        subprocess.run(f"rm {outpath}{base}.bowtie2.index.* 2> /dev/null", shell=True)
        subprocess.run(f"rm {sam} 2> /dev/null", shell=True)
        if build:
            subprocess.run(f"rm {outpath}{base}.{u}.fasta 2> /dev/null", shell=True)

    return bam

def sam_bam(sam, outpath, threads):
    '''
    SAM -> BAM
    '''
    try:
        temp = sam.rsplit("/",1)[1]
        base = temp.rsplit(".",1)[0]
    except Exception:
        base = sam.rsplit(".",1)[0]

    bam = f'{outpath}{base}.bam'
    subprocess.run(f"samtools view -@ {threads} -h -b {sam} > {bam} 2> /dev/null", shell=True)

    return bam


def sort_bam(bam, outpath, threads, clean, clean_bam):
    '''
    Sort BAM file. 
    '''
    if int(os.stat(bam).st_size) == 0:
        return False, None
    # check if anything aligned
    check_aligned = subprocess.check_output(f"samtools view -@ {threads} {bam} | head -n 1", shell=True)
    if len(check_aligned) == 0:
        return False, None
    
    try:
        temp = bam.rsplit("/",1)[1]
        base = temp.rsplit(".",1)[0]
    except Exception:
        base = bam.rsplit(".",1)[0]

    sort_check = False
    if clean_bam: # only need to check if bamfile was the input
        try:
            check = subprocess.check_output(f"samtools view -@ {threads} -H {bam} | grep '@HD'", shell=True)
            if "coordinate" in str(check):
                sort_check = True
        except Exception:
            # no @HD line, retain sort_check = False
            pass
    if not sort_check:
        subprocess.run(f"samtools sort -@ {threads} -o {outpath}{base}.sorted.bam {bam} 2> /dev/null", shell=True)

        if clean and not clean_bam:
            subprocess.run(f'rm {bam}', shell=True)

        bam = f'{outpath}{base}.sorted.bam'

    return True, bam

@jit(nopython=True)
def get_len(start,stop):
    '''
    Length of prophages.
    Python is zero-based and end-exclusive in range(),
    so even though start = start-1, leave this math
    as stop-start to include the start and stop base precisely
    '''
    return stop-start


@jit(nopython=True)
def reverse_and_zero_based(start, stop):
    '''
    Reverse start and stop if the prophage is in reverse coordinates.
    -1 to start for Python zero-based format.
    Leave stop for Python range exclusive format.
    '''
    if start > stop:
        start, stop = stop, start
    start -= 1
    return start,stop


def process_coordinates_file(vibe, spaces, vibe_header):
    '''
    Convert coordinates file into usable dictionary and list of names.
    '''
    prophage_dict = {}
    prophage_lengths = {}
    prophage_dict_frags = {}
    if vibe_header == 'vibrant header':
        x = 5
        y = 6
    elif vibe_header == 'custom header':
        x = 2
        y = 3

    with open(vibe, 'r') as phage_list:
        next(phage_list)
        for line in phage_list:
            line = line.strip('\n').split('\t')
            if not line: continue
            name = line[0]
            frag = line[1]

            if spaces:
                name = name.replace(" ", "~!!~")
                frag = frag.replace(" ", "~!!~")

            start = int(line[x])
            stop = int(line[y])
            start, stop = reverse_and_zero_based(start, stop)
            prophage_dict.setdefault(name, []).append((frag, start, stop))
            prophage_lengths[frag] = get_len(start,stop) # stop-start+1 # Python zero-based
            prophage_dict_frags[frag] = name
        
        genomes_full = set(list(prophage_dict.keys()))
        prophages = len(prophage_dict_frags.keys())

        return True, prophage_dict, prophage_dict_frags, genomes_full, prophages, prophage_lengths

def get_lengths(fasta, spaces, genomes_full):
    '''
    Get genome lengths for initialization of np.zeros(length)
    '''
    lengths = {}
    if spaces:
        for name,seq in fasta_parse(fasta):
            if name in genomes_full:
                name = name.replace(" ", "~!!~")
                lengths[name] = len(seq)
    else:
        for name,seq in fasta_parse(fasta):
            if name in genomes_full:
                lengths[name] = len(seq)

    return lengths


def extract_coverage(bam, read_id, lengths, mask, outfile, prophage_dict, effect, ratio_cutoff, prophage_dict_frags, prophage_lengths, min_breadth, min_cov, clean, spaces):
    '''
    Code mainly from vRhyme (same author).
    Extract coverage information from BAM file.
    '''
    with open(outfile, 'w') as output:
        output.write("prophage\thost\tactive\tCohenD\tprophage-host_ratio\tmean_difference\tprophage_len\tprophage_mean_cov\tprophage_median_cov\tprophage_sd_cov\tprophage_cov_breadth\thost_len\thost_mean_cov\thost_median_cov\thost_sd_cov\n")
    written = []
    total = 0

    bai = False
    if not os.path.exists(bam + '.bai'):
        subprocess.run(f'samtools index {bam}', shell=True)
        bai = True
    
    if mask != 0:
        mask_f = mask
        mask_r = -mask
    else:
        mask_f = False
        mask_r = False

    bamfile = pysam.AlignmentFile(bam, "rb")

    if read_id != 0: 
        for x in bamfile.fetch(until_eof=True):
            genome = x.reference_name
            length = lengths.get(genome,False)
            if length: # not in keep
                try:
                    if genome != prev:
                        if mask_f:
                            depth[:mask_f] = np.nan
                            depth[mask_r:] = np.nan
                        prophage_covs, avg, med, sd, eff_host = coverage_stats(depth, prev, prophage_dict, prophage_lengths, min_cov, length)
                        written, total = write_coverages(outfile, prophage_covs, prev, avg, med, sd, eff_host, effect, ratio_cutoff, written, total, min_breadth, min_cov, spaces)
                        depth = np.zeros(length)
                except NameError:
                    # prev not defined
                    depth = np.zeros(length)

                prev = genome

                ed = 0
                for t in x.tags:
                    if t[0] == 'NM':
                        ed = t[1]
                        break
                rl = x.query_length

                start = x.reference_start # 0-based
                end = x.reference_end
                if end:
                    depth = add_depth(depth, start, end, ed, rl, read_id)

        if length: # not in keep
            # last one
            if mask_f:
                depth[:mask_f] = np.nan
                depth[mask_r:] = np.nan
            prophage_covs, avg, med, sd, eff_host = coverage_stats(depth, prev, prophage_dict, prophage_lengths, min_cov, length)
            written, total = write_coverages(outfile, prophage_covs, prev, avg, med, sd, eff_host, effect, ratio_cutoff, written, total, min_breadth, min_cov, spaces)
            depth = None
    else: # no mismatches
        for x in bamfile.fetch(until_eof=True):
            genome = x.reference_name
            length = lengths.get(genome,False)
            if length: # not in keep
                try:
                    if genome != prev:
                        if mask_f:
                            depth[:mask_f] = np.nan
                            depth[mask_r:] = np.nan
                        prophage_covs, avg, med, sd, eff_host = coverage_stats(depth, prev, prophage_dict, prophage_lengths, min_cov, length)
                        written, total = write_coverages(outfile, prophage_covs, prev, avg, med, sd, eff_host, effect, ratio_cutoff, written, total, min_breadth, min_cov, spaces)
                        depth = np.zeros(length)
                except NameError:
                    depth = np.zeros(length)

                prev = genome
                start = x.reference_start # 0-based
                end = x.reference_end
                if end:
                    depth = add_depth_no_ed(depth, start, end)

        if length: # not in keep
            # last one
            if mask_f:
                depth[:mask_f] = np.nan
                depth[mask_r:] = np.nan
            prophage_covs, avg, med, sd, eff_host = coverage_stats(depth, prev, prophage_dict, prophage_lengths, min_cov, length)
            written, total = write_coverages(outfile, prophage_covs, prev, avg, med, sd, eff_host, effect, ratio_cutoff, written, total, min_breadth, min_cov, spaces)
            depth = None

    bamfile.close()
    if clean and bai:
        subprocess.run(f'rm {bam}.bai', shell=True)

    include_zeros(outfile, prophage_lengths, written, prophage_dict_frags, lengths, spaces)

    return total

@jit(nopython=True)
def cohenD(phage_mean, phage_sd, host_mean, host_sd):
    """
    Cohen's d equation
    """
    try:
        pool = ((phage_sd**2+host_sd**2)/2)**0.5
        d = abs((host_mean-phage_mean)/pool)
    except Exception:
        # host has 0 coverage because length is < 1000
        d = 0
    return d

@jit(nopython=True)
def activity(phage_mean,cov_depth,avg,d,effect,ratio_cutoff,min_breadth,min_cov,total):
    '''
    Determine activity based on cutoffs
    '''
    try:
        ratio = phage_mean/avg
    except Exception:
        ratio = phage_mean
    diff = phage_mean-avg
    active = 'dormant'
    if d >= effect and ratio >= ratio_cutoff:
        if cov_depth >= min_breadth and phage_mean >= min_cov:
            active = 'active'
            total += 1
        else:
            active = 'ambiguous'
    return active,total,diff,ratio

def write_coverages(outfile, prophage_covs, host, avg, med, sd, length, effect, ratio_cutoff, written, total, min_breadth, min_cov, spaces):
    '''
    Within the BAM loop.
    Perform statistical analyses and write out final results.
    '''
    with open(outfile, 'a') as output:
        for key,p in prophage_covs.items():
            phage_mean,phage_med,phage_sd,l,cov_depth = p
            d = cohenD(phage_mean, phage_sd, avg, sd)
            active,total,diff,ratio = activity(phage_mean,cov_depth,avg,d,effect,ratio_cutoff,min_breadth,min_cov,total)
            if spaces:
                key = key.replace("~!!~", " ")
            output.write(f'{key}\t{host}\t{active}\t{d}\t{ratio}\t{diff}\t{l}\t{phage_mean}\t{phage_med}\t{phage_sd}\t{cov_depth}\t{length}\t{avg}\t{med}\t{sd}\n')
            written.append(key)
    return written, total

def include_zeros(outfile, prophage_lengths, written, prophage_dict_frags, lengths, spaces):
    '''
    Write out all the prophages not identified within the BAM file
    '''
    written = set(written)
    with open(outfile, 'a') as output:
        for key,host in prophage_dict_frags.items():
            if spaces:
                key = key.replace("~!!~", " ")
            if key in written: continue
            l = prophage_lengths[key]
            length = lengths[host]
            output.write(f'{key}\t{host}\tno\tNA\tNA\tNA\t{l}\tNA\tNA\tNA\tNA\t{length}\tNA\tNA\tNA\n')
