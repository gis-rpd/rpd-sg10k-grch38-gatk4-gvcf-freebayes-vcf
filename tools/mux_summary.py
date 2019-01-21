#!/usr/bin/env python3
"""Creating QC CSV for GATK4 analysed SG10K health batch

Scanned directory is expected to contain one folder per sample and one indexcov folder
"""

import subprocess
import logging
import sys
import argparse
import os
import csv
import glob

#import yaml

__author__ = "Andreas WILM"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2018 Genome Institute of Singapore"
__license__ = "The MIT License (MIT)"


# global logger
LOGGER = logging.getLogger(__name__)
HANDLER = logging.StreamHandler()
HANDLER.setFormatter(logging.Formatter(
    '[{asctime}] {levelname:8s} {filename} {message}', style='{'))
LOGGER.addHandler(HANDLER)


# indexcov gender mapping
PED_SEX_MAP = {'0': 'unknown', '1': 'male', '2': 'female', '-9': 'undetected'}
# samtools stats fields to report
SAMTOOLS_STATS_KEYS = ["raw total sequences:",
    #"reads properly paired:",
    'reads properly paired [%]',# ours
    #"reads duplicated:",
    'reads duplicated [%]', # ours
    #"bases mapped (cigar):",
    #'% bases mapped',
    #'base coverage',
    "error rate:",
    'insert size average:',
    'insert size standard deviation:']
# verifybamid ethnicities
ETHNICITIES = ['CHS', 'INS', 'MAS']


def parse_summary_from_stats(statsfile):
    """FIXME:add-doc"""

    sn = dict()
    with open(statsfile) as fh:
        for line in fh:
            if not line.startswith("SN"):
                continue
            ls = line.strip().split("\t")[1:]
            k, v = ls[0], ls[1]
            if "." in ls[1]:
                v = float(v)
            else:
                v = int(float(v))# int(float(x)) helps with scientific notation
            sn[k] = v
    sn['bases mapped [%]'] = 100 * sn["bases mapped (cigar):"]/float(sn["raw total sequences:"] * sn["average length:"])
    sn['reads properly paired [%]'] = 100 * sn["reads properly paired:"]/float(sn["raw total sequences:"])
    sn['reads duplicated [%]'] = 100 * sn["reads duplicated:"]/float(sn["raw total sequences:"])
    #sn['base coverage'] = sn["bases mapped (cigar):"]/float(genome_size)
    return sn



def parse_selfsm(selfsm_file):
    """FIXME:add-doc"""

    with open(selfsm_file) as fh:
        header = fh.readline()[1:].split()
        #print(headers)
        values = fh.readline().split()
        #print(values)
    d = dict(zip(header, values))
    for k in ['AVG_DP', 'FREELK1', 'FREELK0', 'FREEMIX']:
        try:
            d[k] = float(d[k])
        except ValueError:
            pass
    for k in ['#SNPS', '#READS']:
        try:
            d[k] = int(d[k])
        except ValueError:
            pass
    return d


def samples_from_ped(pedfile):
    samples = dict()
    with open(pedfile) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            s = row['sample_id']
            del row['sample_id']
            samples[s] = row
    return samples


def main():
    """main
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help="Increase verbosity")
    parser.add_argument('-q', '--quiet', action='count', default=0,
                        help="Decrease verbosity")
    parser.add_argument('-d', '--pubdir', required=True)
    args = parser.parse_args()
    # Repeateable -v and -q for setting logging level.
    # See https://www.reddit.com/r/Python/comments/3nctlm/what_python_tools_should_i_be_using_on_every/
    # and https://gist.github.com/andreas-wilm/b6031a84a33e652680d4
    # script -vv -> DEBUG
    # script -v -> INFO
    # script -> WARNING
    # script -q -> ERROR
    # script -qq -> CRITICAL
    # script -qqq -> no logging at all
    LOGGER.setLevel(logging.WARN + 10*args.quiet - 10*args.verbose)

    pedfile = os.path.join(args.pubdir, "indexcov/all/all-indexcov.ped")
    assert os.path.exists(pedfile), ("Expected file %s missing" % pedfile)
    samples = samples_from_ped(pedfile)

    out = dict()
    for sample_id, v in samples.items():
        assert sample_id not in out
        out[sample_id] = dict()

        out[sample_id]['gender'] = PED_SEX_MAP[v['sex']]

        sample_dir = os.path.join(args.pubdir, sample_id)
        assert os.path.exists(sample_dir), ("Missing expected directory %s" % sample_dir)
        #print(sample_dir)

        # parse samtools stats and only keep values of interest
        ststatsfile_matches = glob.glob(os.path.join(sample_dir, "stats/*.stats"))
        if len(ststatsfile_matches) == 1:
            ststatsfile = ststatsfile_matches[0]
            ststats = parse_summary_from_stats(ststatsfile)
            for k in SAMTOOLS_STATS_KEYS:
                out[sample_id][k.strip(":")] = ststats[k]
        else:
            LOGGER.warning("No samtools stats file found in %s", sample_dir)
            for k in SAMTOOLS_STATS_KEYS:
                out[sample_id][k.strip(":")] = -1

        # parse verifybamid files for depth and contamination values
        selfsm_files = glob.glob(os.path.join(sample_dir, "verifybamid/*selfSM"))
        selfsm = dict()
        for f in selfsm_files:
            ethnicity = f.split(".")[-2]
            selfsm[ethnicity] = parse_selfsm(f)
        if sorted(list(selfsm.keys())) != sorted(ETHNICITIES):
            LOGGER.warning("Only found contamination files for %s in %s",
                    ', '.join(list(selfsm.keys())), sample_dir)
        for e in ETHNICITIES:
            try:
                cont = selfsm[e]['FREEMIX']
            except KeyError:
                cont = -1
            out[sample_id]['cont. {}'.format(e)] = cont
        avg_dp = set([v['AVG_DP'] for v in selfsm.values()])
        avg_dp = list(avg_dp)[0]
        out[sample_id]['average depth'] = avg_dp

        # coverage as per SOP
        covsop_matches = glob.glob(os.path.join(sample_dir, "*.cov-062017.txt"))
        if len(covsop_matches) == 1:
            covsop_file = covsop_matches[0]
            with open(covsop_file) as fh:
                l = fh.readline()
            covsop = int(l.split()[-1]) / float(3.1*10**9)
        else:
            LOGGER.warning("Expected one coverage file in %s but got %d",
                    sample_dir, len(covsop_matches))
            covsop = -1
        out[sample_id]['coverage (SOP 06-2017)'] = covsop

        # number of vcf records and vcf file size
        vcf_matches = glob.glob(os.path.join(sample_dir, "*fb.vcf.gz"))
        filesize = -1
        if len(vcf_matches) == 1:
            fbvcf = vcf_matches[0]
            filesize = os.path.getsize(fbvcf) >> 20

            cmd = ['bcftools', 'index', '--nrecords', fbvcf]
            try:
                #env = {"PATH": "/data/users/astar/gis/rpd/apps/bcftools-1.3/bin/bcftools"}
                res = subprocess.run(cmd, check=True, capture_output=True)
                #res = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
                num_vars = int(res.stdout.decode())
            except subprocess.CalledProcessError:
                LOGGER.warning("bcftools failed on %s (error was '%s')",
                        fbvcf, res.stderr.decode())
                num_vars = -1

        else:
            LOGGER.warning("Expected one vcf file in %s but got %d",
                    sample_dir, len(vcf_matches))
            num_vars = -1
        out[sample_id]['FB variants'] = num_vars
        out[sample_id]['fb-vcf filesize [MB]'] = filesize

        # cram size
        cram_matches = glob.glob(os.path.join(sample_dir, "*cram"))
        if len(cram_matches) == 1:
            filesize = os.path.getsize(cram_matches[0]) >> 20
        else:
            LOGGER.warning("Expected one cram file in %s but got %d",
                    sample_dir, len(cram_matches))
            filesize = -1
        out[sample_id]['cram filesize [MB]'] = filesize

        # gvcf size
        gvcf_matches = glob.glob(os.path.join(sample_dir, "*g.vcf.gz"))
        if len(gvcf_matches) == 1:
            filesize = os.path.getsize(gvcf_matches[0]) >> 20
        else:
            LOGGER.warning("Expected one gvcf file in %s but got %d",
                    sample_dir, len(gvcf_matches))
            filesize = -1
        out[sample_id]['gvcf filesize [MB]'] = filesize
 


    fieldnames = set([f for d in out.values() for f in d.keys()])
    #yaml.safe_dump_all(out, sys.stdout, canonical=True)#, default_flow_style=False)

    # reuse last sample_id
    fieldnames = ['sample'] + sorted(fieldnames)
    #LOGGER.warning("fieldnames = %s", fieldnames)
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)#, delimiter="\t")
    writer.writeheader()
    for sample_id, d in out.items():
        d['sample'] = sample_id
        writer.writerow(d)


if __name__ == "__main__":
    main()
