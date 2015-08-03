from __future__ import print_function, absolute_import, division
import argparse
import csv
from zither import __version__
from datetime import datetime
import os
import os.path
import pysam
import sys


class _MappingFileStrategy(object):
    def __init__(self, mapping_file):
        self._mapping_file = mapping_file

    def _abs_path(self, path):
        mapping_dir_path = os.path.dirname(self._mapping_file)
        if path != os.path.abspath(path):
            path = os.path.abspath(os.path.join(mapping_dir_path, path))
        return path
        
    def build_sample_bam_mapping(self):
        sample_bam_mapping = {}
        mapping_dir_path = os.path.dirname(self._mapping_file)
        with open(self._mapping_file, 'rb') as tsvfile:
            for sample_name, bam_path in csv.reader(tsvfile, delimiter='\t'):
                sample_bam_mapping[sample_name] = self._abs_path(bam_path)
        return sample_bam_mapping    
            
            # sample_bam_mapping = {k:v for (k,v) in reader}
            # sample_bam_mapping_modified = {}
            # mapping_dir_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..',"examples/mapping_files/")
            # #if the bam is not absolute, then prepend mapping file dir 
            # for bam_path in sample_bam_mapping.values():
                # if bam_path != os.path.abspath(bam_path):
                    # bam_path = mapping_dir_path + bam_path
                    # for sample_name in sample_bam_mapping.keys():
                        # sample_bam_mapping_modified[sample_name] = bam_path
                    # return sample_bam_mapping_modified
                # else:
                    # return sample_bam_mapping

class _MatchingNameStrategy(object):
    def __init__(self, sample_names, input_vcf_path):
        self._sample_names = sample_names
        self._input_vcf_path = input_vcf_path

    def build_sample_bam_mapping(self):
        sample_bam_mapping = {}
        bam_dir = os.path.dirname(self._input_vcf_path)
        for sample_name in self._sample_names:
            bam_path = os.path.join(bam_dir, sample_name + ".bam")
            sample_bam_mapping[sample_name] = bam_path
        return sample_bam_mapping

        
def _get_sample_bam_strategy(args):
    if args.mapping_file:
        return _MappingFileStrategy(args.mapping_file)
    else:
        sample_names = _get_sample_names(args.input_vcf)
        return _MatchingNameStrategy(sample_names, args.input_vcf)
            
class _BamReader(object):
    def __init__(self, bam_file_name):
        self._bam_file_name = bam_file_name
        self._bam_file = pysam.AlignmentFile(bam_file_name, "rb" )

    def __eq__(self, other):
        return isinstance(other, _BamReader) and self._bam_file_name == other._bam_file_name

    def __hash__(self):
        return hash(self._bam_file_name)
        
    def get_depth_and_alt_freq(self, chrom, position_one_based, ref, alt):
        BASE_INDEX = {'A':0, 'C':1, 'G':2, 'T':3}        
        alt = alt.upper()
        AF = '.'
        position_zero_based = position_one_based - 1
        coverage = self._bam_file.count_coverage(chr=chrom,
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

def _get_sample_names(input_vcf):
    with open(input_vcf, 'r') as input_file:
        column_header = None
        sample_names = []
        for line in input_file.readlines():
            if not line.startswith("##"):
                if line.startswith("#"):
                    column_header = line.rstrip()
                    column_fields = column_header.split("\t")
                    n = len(column_fields)
                    sample_names = column_fields[9:(n+1)]
        return sample_names

def _build_reader_dict(sample_bam_mapping):
    readers_dict = {}
    for (sample, bam_file) in sample_bam_mapping.items():
        readers_dict[sample] = _BamReader(bam_file)
    return readers_dict

def _create_vcf(input_vcf, sample_reader_dict):
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

        for line in input_file.readlines():
            if not line.startswith("##"):
                if line.startswith("#"):
                    column_header = line.rstrip()
                    sample_names = column_header.split("\t")[9:]
                    print(zither)
                    print(column_header)
                    
                else:
                    vcf_fields = line.rstrip("\n").split("\t")[0:5]
                    (CHROM, POS, ID, REF, ALT) = vcf_fields
                    FORMAT = "BDP:BAF"
                    vcf_fields.append('.')
                    vcf_fields.append('.')
                    vcf_fields.append('.')
                    vcf_fields.append(FORMAT)
                    for sample_name in sample_names:
                        bam_reader = sample_reader_dict[sample_name]
                        [depth,FA] = bam_reader.get_depth_and_alt_freq(CHROM, int(POS), REF, ALT)
                        sample_field = [str(depth),FA]
                        sample_field_joint = ':'.join(sample_field)
                        vcf_fields.append(sample_field_joint)
                    a = '\t'.join(vcf_fields)
                    print(a)
                        
def _parse_command_line_args(arguments):
    parser = argparse.ArgumentParser(usage="zither [-h] [-V] input_vcf input_bam",
        description='''For all positions in VCF, pull raw depths and alt freqs from BAM file, writing output as new VCF to stdout. Type 'zither -h' for help''')

    parser.add_argument("-V", "--version", action='version', version=__version__)
    parser.add_argument('input_vcf', help="Path to input VCFs; all record locations will appear in output file")
    parser.add_argument('--mapping_file', help="Path to tab delimited list of VCF_sample_names and BAM_file_names")
    args = parser.parse_args(arguments)
    return args
    
        
def main(command_line_args=sys.argv):
    args = _parse_command_line_args(command_line_args[1:])
    strategy = _get_sample_bam_strategy(args)
    sample_bam_mapping = strategy.build_sample_bam_mapping()
    reader_dict = _build_reader_dict(sample_bam_mapping)
    _create_vcf(args.input_vcf, reader_dict)

if __name__ == '__main__':
    main()
