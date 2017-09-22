#!/usr/bin/env python

from __future__ import print_function, division
import sys, os, datetime

if 'VIRTUAL_ENV' in os.environ:
    print('found a virtualenv -- activating: {}'.format(os.environ['VIRTUAL_ENV']), file=sys.stderr)
    activation_script = os.path.join(os.environ['VIRTUAL_ENV'], 'bin', 'activate_this.py')
    execfile(activation_script, dict(__file__=activation_script))

import click
from cyvcf2 import VCF, Writer
import numpy as np

def log(msg):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %T")
    print('[{}] {}'.format(timestamp, msg), file=sys.stderr)

def fail_variant(variant, tag):
    if variant.FILTER is not None:
        variant.FILTER = variant.FILTER.split(';') + [tag]
    else:
        variant.FILTER = tag

def will_evaluate(variant, exclude_filters, exclude_fields):
    if variant.FILTER is not None:
        filters = set(variant.FILTER.split(';'))
        for filt_string in exclude_filters:
            if filt_string in filters:
                return False
    info_tags = set([ t[0] for t in variant.INFO ])
    for field in exclude_fields:
        if field in info_tags:
            return False
    return True

def ab_filter_rescue(variant, ab, min_ab, max_ab, dp, min_dp, filter_tag):
    if (ab >= min_ab and
            ab <= max_ab and
            dp >= min_dp):
        if variant.FILTER is not None:
            variant.INFO['AB_RESCUED'] = variant.FILTER
            variant.FILTER = 'PASS'
            return 1
        else:
            return 2
    else:
        self.fail_variant(filter_tag)
        return 0

def is_biallelic(variant):
    return True if len(variant.ALT) == 1 else False

def compute_allelic_balances(variant):
    # initial values
    (het_ab, het_hom_alt_ab) = (0.0 , 0.0)

    # gt_types is array of (0,1,2,3) == (HOM_REF, HET, UNKNOWN, HOM_ALT)

    # "0/1" genotype cases
    het_mask = variant.gt_types == 1
    ref_het_counts = np.sum(variant.format('AD')[het_mask][:,0])
    alt_het_counts = np.sum(variant.format('AD')[het_mask][:,1])

    # "1/1" genotype cases
    hom_alt_mask = variant.gt_types == 3
    ref_hom_alt_counts = np.sum(variant.format('AD')[hom_alt_mask][:,0])
    alt_hom_alt_counts = np.sum(variant.format('AD')[hom_alt_mask][:,1])

    total_het_counts = alt_het_counts + ref_het_counts
    if total_het_counts != 0 :
        het_ab = alt_het_counts / total_het_counts

    total_het_hom_alt_counts = total_het_counts + alt_hom_alt_counts + ref_hom_alt_counts
    if total_het_hom_alt_counts != 0 :
        numerator = alt_het_counts + (0.5 * alt_hom_alt_counts)
        het_hom_alt_ab = numerator / total_het_hom_alt_counts

    return (het_ab, het_hom_alt_ab, total_het_counts, total_het_hom_alt_counts)

def update_variant(variant, het_ab, het_hom_alt_ab, total_het_count, total_het_hom_alt_count):
    variant.INFO['HetAB'] = '{:.4f}'.format(het_ab)
    variant.INFO['HetHomAltAB'] = '{:.4f}'.format(het_hom_alt_ab)
    variant.INFO['HetAB_DP'] = '{}'.format(total_het_count)
    variant.INFO['HetHomAltAB_DP'] = '{}'.format(total_het_hom_alt_count)
    return variant

def filter_description(ab_tag, depth_tag, min_ab, max_ab, min_dp, exclude_filters, exclude_fields):
        desc = 'Failed allele balance filter. {1} <= {0} <= {2} and {4} >= {3}.'.format(ab_tag, min_ab, max_ab, min_dp, depth_tag)
        if exclude_filters:
            desc += ' Ignored sites with FILTER containing: {0}.'.format(','.join(exclude_filters))
        if exclude_fields:
            desc += ' Ignored sites with any of the following tags in the INFO field: {0}.'.format(','.join(exclude_fields))
        return desc

def annotate_allelic_balance(vcffile, region, het_only, min_ab, max_ab, min_dp, exclude_filters, exclude_fields):
    vcf = VCF(vcffile, strict_gt=True)

    header_hetab_param_info = {
        'ID' : 'HetAB',
        'Description' : 'heterozygous genotype allele balance',
        'Type' : 'Float',
        'Number' : '1'
    }

    header_hetab_dp_param_info = {
        'ID' : 'HetAB_DP',
        'Description' : 'heterozygous genotype read depth',
        'Type' : 'Integer',
        'Number' : '1'
    }

    header_het_hom_alt_ab_param_info = {
        'ID' : 'HetHomAltAB',
        'Description' : 'heterozygous + homozygous ALT genotype allele balance',
        'Type' : 'Float',
        'Number' : '1'
    }

    header_het_hom_alt_ab_dp_param_info = {
        'ID' : 'HetHomAltAB_DP',
        'Description' : 'heterozygous + homozygous ALT genotype read depth',
        'Type' : 'Integer',
        'Number' : '1'
    }

    header_rescued_info = {
        'ID' : 'AB_RESCUED',
        'Description' : 'Filter status before rescued based on allele balance',
        'Type' : 'String',
        'Number' : '1'
    }

    ab_tag = 'HetHomAltAB'
    depth_tag = 'HetHomAltAB_DP'
    if het_only:
        ab_tag = 'HetAB'
        depth_tag = 'HetAB_DP'

    filter_tag = 'AB_FILTERED{0}to{1}'.format(min_ab, max_ab)

    header_ab_filter = {
        'ID' : filter_tag,
        'Description' : filter_description(ab_tag, depth_tag, min_ab, max_ab, min_dp, exclude_filters, exclude_fields)
    }

    vcf.add_info_to_header(header_hetab_param_info)
    vcf.add_info_to_header(header_hetab_dp_param_info)
    vcf.add_info_to_header(header_het_hom_alt_ab_param_info)
    vcf.add_info_to_header(header_het_hom_alt_ab_dp_param_info)
    vcf.add_info_to_header(header_rescued_info)
    vcf.add_filter_to_header(header_ab_filter)
    out = Writer('-', vcf)
    (total_sites, noted_sites, evaluated_sites, passed_sites, rescued_sites, filtered_sites) = (0, 0, 0, 0, 0)

    for variant in vcf(region):
        total_sites += 1
        if is_biallelic(variant):
            noted_sites += 1
            (hetab, het_hom_alt_ab, total_het_count, total_het_hom_alt_count) = compute_allelic_balances(variant)
            variant = update_variant(variant, hetab, het_hom_alt_ab, total_het_count, total_het_hom_alt_count)
            if will_evaluate(variant, exclude_filters, exclude_fields):
                evaluated_sites += 1
                if het_only:
                    rv = ab_filter_rescue(variant, hetab, min_ab, max_ab, total_het_count, min_dp, filter_tag)  
                else:
                    rv = ab_filter_rescue(variant, het_hom_alt_ab, min_ab, max_ab, total_het_hom_alt_counts, min_dp, filter_tag)
                if rv == 1:
                    rescued_sites += 1
                elif rv == 0:
                    filtered_sites += 1
                elif rv == 2:
                    passed_sites += 1

        out.write_record(variant)

    out.close()
    msg = "Annotated {} out of a possible {} sites.\nEvaluated {} for filtering or rescue.\nPassed {} previously PASS or unfiltered.\nRescued {} previously filtered.\nFailed {} due to allele balance"
    msg = msg.format(noted_sites, total_sites, evaluated_sites, passed_sites, rescued_sites, filtered_sites)
    log(msg)

@click.command()
@click.option('--region', default=None, type=click.STRING,
        help="a chromosome region to limit to [default: None]")
@click.option('--min-allele-balance', default=0.3, type=click.FLOAT,
        help='Minimum allele balance value to rescue a variant')
@click.option('--max-allele-balance', default=0.7, type=click.FLOAT,
        help='Maximum allele balance value to rescue a variant')
@click.option('--use-het-only', is_flag=True,
        help='Filter using allele balance calculated from heterozygous genotypes only')
@click.option('--min-depth', default=10, type=click.INT,
        help='Minimum depth of non-reference samples to rescue a variant')
@click.option('--exclude-filters', default=['MISSING'], multiple=True, type=click.STRING,
        help='Filters that cannot be rescued')
@click.option('--exclude-fields', multiple=True, type=click.STRING,
        help='INFO fields that preclude a line from being rescued')
@click.argument('vcfs', nargs=-1, type=click.Path())
def main(region, min_ab, max_ab, het_only, min_dp, exclude_filters, exclude_fields, vcfs):
    for vcf in vcfs:
        log('processing: {}'.format(vcf))
        annotate_allelic_balance(vcf, region, het_only, min_ab, max_ab, min_dp, exclude_filters, exclude_fields)
    log("All Done!")

if __name__ == "__main__":
    main()
