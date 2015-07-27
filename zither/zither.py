from __future__ import print_function, absolute_import, division
import argparse
from zither import __version__
from datetime import datetime
import os
import pysam
import sys


class _BamReader(object):
    def __init__(self, bam_file_name):
        self.samfile = pysam.AlignmentFile(bam_file_name, "rb" )

    def get_depth_and_alt_freq(self, chrom, position_one_based, ref, alt):
        BASE_INDEX = {'A':0, 'C':1, 'G':2, 'T':3}        
        alt = alt.upper()
        AF = '.'
        position_zero_based = position_one_based - 1
        coverage = self.samfile.count_coverage(chr=chrom,
                                          start=position_zero_based,
                                          stop=position_one_based,
                                          quality_threshold=-1,
                                          read_callback='nofilter')
        
        total_depth = coverage[0][0] + coverage[1][0] + coverage[2][0]+ coverage[3][0]
        try:
            variant_count = coverage[BASE_INDEX[alt]][0]
            if total_depth and len(ref)==1:
                AF = str(variant_count/total_depth)
        except KeyError:
            pass #AF remains null

        return (total_depth, AF)

def _create_vcf(input_vcf):
    vcf_headers = \
'''##fileformat=VCFv4.1
##FORMAT=<ID=BDP,Number=1,Type=Integer,Description="BAM depth">
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="BAM alt frequency">'''

    with open(input_vcf, 'r') as input_file:
        print(vcf_headers)
        sample_field = []
        column_header = None
        now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        cwd = os.path.dirname(os.getcwd())
        command = ' '.join(sys.argv)
        zither = '##zither=<timestamp="{}",command="{}",cwd="{}",version="{}">'.format(now, command, cwd, __version__)
        bamReader_A = None
        bamReader_B = None
        for line in input_file.readlines():
            if not line.startswith("##"):
                if line.startswith("#"):
                    column_header = line.rstrip()
                    column_fields = column_header.split("\t")
                    sample_name_A = column_fields[9]
                    bamReader_A = _BamReader(os.path.join(os.path.dirname(input_vcf), sample_name_A + '.bam'))
                    if len(column_fields) == 11:
                        sample_name_B = column_fields[10]
                        bamReader_B = _BamReader(os.path.join(os.path.dirname(input_vcf), sample_name_B + '.bam'))

                    
                    print(zither)
                    print(column_header)
                    
                else:
                    vcf_fields = line.rstrip("\n").split("\t")
                    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = vcf_fields[0:8]
                    FORMAT = "BDP:BAF"
                    vcf_fields[5] = '.'
                    vcf_fields[6] = '.'
                    vcf_fields[7] = '.'
                    vcf_fields[8] = FORMAT
                    vcf_fields[9] = ''
                    if bamReader_A:
                        [depth,FA] = bamReader_A.get_depth_and_alt_freq(CHROM, int(POS), REF, ALT)
                        sample_field = [str(depth),FA]
                        sample_field_joint = ':'.join(sample_field)
                        vcf_fields[9] += sample_field_joint
                    if bamReader_B:
                        [depth,FA] = bamReader_B.get_depth_and_alt_freq(CHROM, int(POS), REF, ALT)
                        sample_field = [str(depth),FA]
                        sample_field_joint = ':'.join(sample_field)
                        vcf_fields[9] +='\t'+ sample_field_joint
                    a = '\t'.join(vcf_fields[0:10])
                    print(a)
                        
def _parse_command_line_args(arguments):
    parser = argparse.ArgumentParser(usage="zither [-h] [-V] input_vcf input_bam",
        description='''For all positions in VCF, pull raw depths and alt freqs from BAM file, writing output as new VCF to stdout. Type 'zither -h' for help''')

    parser.add_argument("-V", "--version", action='version', version=__version__)
    parser.add_argument('input_vcf', help="Path to input VCFs; all record locations will appear in output file")
    parser.add_argument('input_bam', help="Path to indexed BAM used to calculate raw depth and frequency")
    args = parser.parse_args(arguments)
    return args
    
        
def main():
    args = _parse_command_line_args(sys.argv[1:])
    _create_vcf(args.input_vcf)

if __name__ == '__main__':
    main()
