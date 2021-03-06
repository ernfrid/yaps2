#!/usr/bin/env python

from __future__ import print_function, division
import itertools, re, sys, os

if 'VIRTUAL_ENV' in os.environ:
    print('found a virtualenv -- activating: {}'.format(os.environ['VIRTUAL_ENV']))
    activation_script = os.path.join(os.environ['VIRTUAL_ENV'], 'bin', 'activate_this.py')
    execfile(activation_script, dict(__file__=activation_script))

import click
from cyvcf2 import VCF, Writer
from pandas import Series, DataFrame
import pandas as pd

def variant_missing_criteria(threshold, variant_pct):
    if variant_pct < threshold:
        return 'pass'
    else:
        return 'fail'

def get_alts(variant):
    return ','.join(variant.ALT)

def compute_missingness(variant):
    alleles = Series(list(itertools.chain.from_iterable([re.split('[/|]', x) for x in variant.gt_bases])))
    (missing, total) = ( alleles[alleles == '.'].size, alleles.size )
    missingness = (missing/total) * 100
    return (missingness, missing, total)

def load_snp_database(dbfile):
    vcf = VCF(dbfile)
    db  = set((v.CHROM, v.POS, v.REF, get_alts(v)) for v in vcf)
    return db

def check_snp_db_membership(variant, snp_db):
    snp = (variant.CHROM, variant.POS, variant.REF, get_alts(variant))
    if snp in snp_db:
        return True
    else:
        return False

def compute_alt_allele_frequency(variant):
    alleles = Series(list(itertools.chain.from_iterable([re.split('[/|]', x) for x in variant.gt_bases])))
    alt_freqs = []
    total_non_missing_alleles = alleles[alleles != '.'].size
    for alt in variant.ALT:
        alt_count = alleles[alleles == alt].size
        alt_freq = alt_count / total_non_missing_alleles
        alt_freqs.append(alt_freq)
    alt_freqs = [ "{:.7f}".format(i) for i in alt_freqs ]
    return ','.join(alt_freqs)

def record_stats(file, variant, verdict, missing_pct, missing, total, snp_db):
    filt = variant.FILTER if variant.FILTER else 'PASS'
    db_membership = check_snp_db_membership(variant, snp_db)

    alt_freq = variant.INFO.get('AF', 'NA')
    if alt_freq == 'NA':
        alt_freq = compute_alt_allele_frequency(variant)
        msg = "Missing AF for POS: {} (computed as: {})"
        print(msg.format(variant.POS, alt_freq), file=sys.stderr)

    data = [
        variant.CHROM,
        variant.POS,
        variant.REF,
        get_alts(variant),
        filt,
        alt_freq,
        missing,
        total,
        missing_pct,
        db_membership,
        verdict,
    ]
    data = [str(i) for i in data]
    print("\t".join(data), file=file)

def print_stats_headers(file):
    headers = (
        'CHROM',
        'POS',
        'REF',
        'ALT',
        'ORIG_FILTER',
        'AF',
        'missing_allele_count',
        'total_allele_count',
        'missing_pct',
        'in_SNP_DB',
        'missing_verdict',
    )
    line = "{}{}".format('#', "\t".join(headers))
    print(line, file=file)

def update_variant(variant, verdict):
    if verdict == "pass":
        if not variant.FILTER:
            variant.FILTER = ['PASS']
    else:
        if (not variant.FILTER) or (variant.FILTER == 'PASS'):
            variant.FILTER = ['MISSING']
        else:
            variant.FILTER = [ variant.FILTER, 'MISSING' ]

    return variant

def mark_missing_sites(vcffile, region, missing_threshold, soft_filter, statsfile, snp_db):
    vcf = VCF(vcffile)
    header_param = {
        'ID' : 'MISSING',
        'Description' : 'failed variant site missingness threshold ({} %)'.format(missing_threshold)
    }
    vcf.add_filter_to_header(header_param)
    out = Writer('-', vcf)
    (total_sites, noted_sites) = (0, 0)
    with open(statsfile, "w") as fstats:
        print_stats_headers(fstats)
        for variant in vcf(region):
            total_sites += 1
            (missing_pct, missing, total) = compute_missingness(variant)
            verdict = variant_missing_criteria(missing_threshold, missing_pct)
            record_stats(fstats, variant, verdict, missing_pct, missing, total, snp_db)
            variant = update_variant(variant, verdict)
            if verdict == "pass":
                noted_sites += 1
                out.write_record(variant)
            elif verdict == "fail" and soft_filter:
                out.write_record(variant)
    out.close()
    msg = "After filtering, kept {} out of a possible {} Sites ({})"
    msg = msg.format(noted_sites, total_sites, 'pass')
    print(msg, file=sys.stderr)

@click.command()
@click.option('--missing-threshold', type=float, default=2.0,
        help="the missingness threshold to discard [default: 2.0]")
@click.option('--soft', is_flag=True, default=False,
        help="soft filtering -- keep all variants and update FILTER field [default: False]")
@click.option('--stats', default="stats.out", type=click.Path(),
        help="an output file to write site stats to [default: 'stats.out']")
@click.option('--db', required=True, type=click.Path(),
        help="a database [e.g. 'dbSNP', or 'HapMap3'] to check membership")
@click.option('--region', default=None, type=click.STRING,
        help="a chromosome region to limit to [default: None]")
@click.argument('vcfs', nargs=-1, type=click.Path())
def main(missing_threshold, region, stats, db, soft, vcfs):
    snp_db = load_snp_database(db)
    for vcf in vcfs:
        mark_missing_sites(vcf, region, missing_threshold, soft, stats, snp_db)
    print("All Done!", file=sys.stderr)

if __name__ == "__main__":
    main()
