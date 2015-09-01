
##   Copyright 2014 Bioinformatics Core, University of Michigan
##
##   Licensed under the Apache License, Version 2.0 (the "License");
##   you may not use this file except in compliance with the License.
##   You may obtain a copy of the License at
##
##       http://www.apache.org/licenses/LICENSE-2.0
##
##   Unless required by applicable law or agreed to in writing, software
##   distributed under the License is distributed on an "AS IS" BASIS,
##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##   See the License for the specific language governing permissions and
##   limitations under the License.

#pylint: disable=invalid-name, too-few-public-methods, too-many-locals
from __future__ import print_function, absolute_import, division
import argparse
import csv
from collections import OrderedDict
from zither import __version__
from datetime import datetime
import os.path
import pysam
import sys

_VCF_FIXED_HEADERS = ["#CHROM",
                      "POS",
                      "ID",
                      "REF",
                      "ALT",
                      "QUAL",
                      "FILTER",
                      "INFO",
                      "FORMAT"]

_NULL = "."

class _ExplicitBamFileStrategy(object):
    def __init__(self, bam_file_path):
        self._bam_file_path = bam_file_path

    def build_sample_bam_mapping(self):
        sample_file = os.path.basename(self._bam_file_path)
        sample_name = os.path.splitext(sample_file)[0]
        return {sample_name: self._bam_file_path}


class _MappingFileStrategy(object):
    def __init__(self, mapping_file):
        self._mapping_file = mapping_file

    def _abs_path(self, path):
        mapping_dir_path = os.path.dirname(self._mapping_file)
        if path != os.path.abspath(path):
            path = os.path.abspath(os.path.join(mapping_dir_path, path))
        return path

    def build_sample_bam_mapping(self):
        sample_bam_mapping = OrderedDict()
        with open(self._mapping_file, 'rb') as tsvfile:
            for sample_name, bam_path in csv.reader(tsvfile, delimiter='\t'):
                sample_bam_mapping[sample_name] = self._abs_path(bam_path)
        return sample_bam_mapping

class _MatchingNameStrategy(object):
    def __init__(self, sample_names, input_vcf_path):
        self._sample_names = sample_names
        self._input_vcf_path = input_vcf_path

    def build_sample_bam_mapping(self):
        sample_bam_mapping = OrderedDict()
        bam_dir = os.path.dirname(self._input_vcf_path)
        for sample_name in self._sample_names:
            bam_path = os.path.join(bam_dir, sample_name + ".bam")
            sample_bam_mapping[sample_name] = bam_path
        return sample_bam_mapping


def _get_sample_bam_strategy(args):
    if args.mapping_file:
        return _MappingFileStrategy(args.mapping_file)
    elif args.bam:
        return _ExplicitBamFileStrategy(args.bam)
    else:
        sample_names = _get_sample_names(args.input_vcf)
        return _MatchingNameStrategy(sample_names, args.input_vcf)


class _PileupStats(object):
    _PYSAM_BASE_INDEX = {'A':0, 'C':1, 'G':2, 'T':3}
    def __init__(self, ref, alt, unfiltered_coverage, filtered_coverage):
        (self.total_depth,
         self.total_af) = self._init_depth_freq(ref, alt, unfiltered_coverage)
        (self.filtered_depth,
         self.filtered_af) = self._init_depth_freq(ref, alt, filtered_coverage)

    def _init_depth_freq(self, ref, alt, coverage):
        alt = alt.upper()
        freq = _NULL
        depth = (coverage[0][0] +
                 coverage[1][0] +
                 coverage[2][0] +
                 coverage[3][0])
        try:
            variant_count = coverage[self._PYSAM_BASE_INDEX[alt]][0]
            if depth and len(ref)==1:
                freq = str(variant_count/depth)
        except KeyError:
            freq = _NULL
        return (depth, freq)


class _BamReader(object):
    def __init__(self, bam_file_name, basecall_quality_cutoff):
        self._bam_file_name = bam_file_name
        #make cutoff inclusive (a cutoff of 20 will include bcq of 20)
        self._basecall_quality_cutoff = basecall_quality_cutoff - 1
        #pylint: disable=no-member
        self._bam_file = pysam.AlignmentFile(bam_file_name, "rb")
        self._filtered_bam_file = pysam.AlignmentFile(bam_file_name, "rb")

    def __eq__(self, other):
        return (isinstance(other,_BamReader) and
                self._bam_file_name == other._bam_file_name and
                self._basecall_quality_cutoff == other._basecall_quality_cutoff)

    def __hash__(self):
        return hash(self._bam_file_name)

    def get_pileup_stats(self, chrom, pos_one_based, ref, alt):
        pos_zero_based = pos_one_based - 1
        try:
            coverage2 = {}
            coverage2["A"]=0
            coverage2["C"]=0
            coverage2["G"]=0
            coverage2["T"]=0

            coverage3 = {}
            coverage3["A"]=0
            coverage3["C"]=0
            coverage3["G"]=0
            coverage3["T"]=0

            pileupcolumns = self._filtered_bam_file.pileup(chrom,
                                                           pos_zero_based,
                                                           pos_one_based,
                                                           stepper='nofilter',
                                                           truncate=True)
            for pileupcolumn in pileupcolumns:
#                if pileupcolumn.reference_pos < pos_zero_based:
#                    continue
#                if pileupcolumn.reference_pos > pos_zero_based:
#                    break
                for read in pileupcolumn.pileups:
                    pos = read.query_position
                    if not read.is_del:
                        base = read.alignment.query_sequence[pos]
                        coverage2[base.upper()] += 1
                        basecall_qual = read.alignment.query_qualities[pos]
                        if basecall_qual > self._basecall_quality_cutoff:
                            coverage3[base.upper()] += 1
            coverage2 = [[coverage2["A"]],
                         [coverage2["C"]],
                         [coverage2["G"]],
                         [coverage2["T"]]]
            coverage3 = [[coverage3["A"]],
                         [coverage3["C"]],
                         [coverage3["G"]],
                         [coverage3["T"]]]

#             coverage = self._bam_file.count_coverage(chr=chrom,
#                                                      start=pos_zero_based,
#                                                      stop=pos_one_based,
#                                                      quality_threshold=-1,
#                                                      read_callback='nofilter')

#             filtered_coverage = self._filtered_bam_file.count_coverage(chr=chrom,
#                                                                        start=pos_zero_based,
#                                                                        stop=pos_one_based,
#                                                                        quality_threshold=self._basecall_quality_cutoff,
#                                                                        read_callback='nofilter')

        except ValueError as samtools_error:
            if str(samtools_error).startswith("invalid reference"):
                coverage2 = [[0], [0], [0], [0]]
                coverage3 = [[0], [0], [0], [0]]
#                 coverage = [[0], [0], [0], [0]]
#                 filtered_coverage = [[0], [0], [0], [0]]
            else:
                raise samtools_error
        return _PileupStats(ref, alt, coverage2, coverage3)

class _Tag(object):
    _METAHEADER = '##FORMAT=<ID={},Number={},Type={},Description="{}">'
    def __init__(self, vcf_id, number, tag_type, description, stats_method):
        self.metaheader = self._METAHEADER.format(vcf_id,
                                                  number,
                                                  tag_type,
                                                  description)
        self.id = vcf_id
        self._get_value_method = stats_method

    def get_value(self, pileup_stats):
        return str(self._get_value_method(pileup_stats))


total_depth = _Tag("ZTDP", "1", "Integer",
                   "Zither total (unfiltered) BAM depth",
                   lambda pileup_stats: pileup_stats.total_depth)
total_af = _Tag("ZTAF", "1", "Float",
                "Zither total (unfiltered) BAM alt frequency",
                lambda pileup_stats: pileup_stats.total_af)
filtered_depth = _Tag("ZFDP", "1", "Integer",
                      "Zither filtered BAM depth",
                      lambda pileup_stats: pileup_stats.filtered_depth)
filtered_af = _Tag("ZFAF", "1", "Float",
                   "Zither filtered BAM alt frequency",
                   lambda pileup_stats: pileup_stats.filtered_af)

DEFAULT_TAGS = [total_depth, total_af, filtered_depth, filtered_af]


def _build_execution_context(argv):
    return OrderedDict([("timestamp",
                         datetime.now().strftime('%Y-%m-%d %H:%M:%S')),
                        ("command",
                         ' '.join(argv)),
                        ("cwd",
                         os.getcwd()),
                        ("version",
                         __version__)])

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

def _build_reader_dict(sample_bam_mapping, args):
    readers_dict = OrderedDict()
    for (sample, bam_file) in sample_bam_mapping.items():
        readers_dict[sample] = _BamReader(bam_file,
                                          int(args.basecall_quality_cutoff))
    return readers_dict

def _build_column_header_line(sample_names):
    column_headers = list(_VCF_FIXED_HEADERS)
    column_headers.extend(sample_names)
    return '\t'.join(column_headers)

def _create_vcf(input_vcf, sample_reader_dict, execution_context, tags=None):
    if tags is None:
        tags = DEFAULT_TAGS
    exec_tags = ['{}="{}"'.format(k,v) for (k,v) in execution_context.items()]
    zither_metaheader = '##zither=<{}>'.format(",".join(exec_tags))
    vcf_headers = ['##fileformat=VCFv4.1']
    vcf_headers.extend([tag.metaheader for tag in tags])
    vcf_headers.append(zither_metaheader)
    FORMAT = ":".join([tag.id for tag in tags])
    with open(input_vcf, 'r') as input_file:

        print("\n".join(vcf_headers))
        print(_build_column_header_line(sample_reader_dict.keys()))
        for line in input_file.readlines():
            if not line.startswith("#"):
                vcf_fields = line.rstrip("\n").split("\t")[0:5]
                (CHROM, POS, dummy, REF, ALT) = vcf_fields
                vcf_fields.append('.')
                vcf_fields.append('.')
                vcf_fields.append('.')
                vcf_fields.append(FORMAT)
                for sample_name in sample_reader_dict.keys():
                    bam_reader = sample_reader_dict[sample_name]
                    pileup_stats = bam_reader.get_pileup_stats(CHROM,
                                                               int(POS),
                                                               REF,
                                                               ALT)
                    sample_field = [tag.get_value(pileup_stats) for tag in tags]
                    sample_field_joint = ':'.join(sample_field)
                    vcf_fields.append(sample_field_joint)
                a = '\t'.join(vcf_fields)
                print(a)

def _parse_command_line_args(arguments):
    parser = argparse.ArgumentParser( \
        usage="zither [-h] [-V] input_vcf input_bam",
        description=("For all positions in VCF, pull raw depths and alt freqs "
                     "from BAM file, writing output as new VCF to stdout. "
                     "Type 'zither -h' for help."))

    parser.add_argument("-V",
                        "--version",
                        action='version',
                        version=__version__)
    parser.add_argument('input_vcf',
                        help="path to input VCFs; all record locations will "
                        "appear in output file")
    parser.add_argument('--bam',
                        help="path to indexed BAM; used to calculate raw depth "
                        "and frequency")
    parser.add_argument('--mapping_file',
                        help="path to tab delimited list of VCF_sample_names "
                        "and BAM_file_names")
    parser.add_argument('--basecall_quality_cutoff',
                        default=20,
                        help="minimum base-call quality to be included. "
                        "Defaults to 20")
    args = parser.parse_args(arguments)
    return args


def main(command_line_args=None):
    if not command_line_args:
        command_line_args = sys.argv
    args = _parse_command_line_args(command_line_args[1:])
    execution_context = _build_execution_context(command_line_args)
    strategy = _get_sample_bam_strategy(args)
    sample_bam_mapping = strategy.build_sample_bam_mapping()
    reader_dict = _build_reader_dict(sample_bam_mapping, args)
    _create_vcf(args.input_vcf, reader_dict, execution_context)

if __name__ == '__main__':
    main(sys.argv)
