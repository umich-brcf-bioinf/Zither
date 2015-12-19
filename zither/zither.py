'''
Zither pulls raw depths and alt freqs from BAM file(s) based on loci in an
existing VCF; writes a new VCF to stdout.

'''
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

#pylint: disable=invalid-name, too-few-public-methods
from __future__ import print_function, absolute_import, division
import argparse
import csv
from collections import OrderedDict, defaultdict
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

_DEFAULT_DEPTH_CUTOFF = 100000
'''For a given position, reads will only be counted up to this cutoff.'''

_DEFAULT_BASECALL_QUALITY_CUTOFF = 20
'''Bases below this quality will be ignored in filtered tag calculations.'''

_DEFAULT_MAPPING_QUALITY_CUTOFF = 20
'''Reads below this mapping quality will be ignored in filtered tag
calculations.'''

class ZitherException(Exception):
    """Base class for all run-time exceptions in this module."""
    def __init__(self, msg, *args):
        #pylint: disable=star-args
        error_msg = msg.format(*[str(i) for i in args])
        super(ZitherException, self).__init__(error_msg)


class ZitherUsageError(ZitherException):
    """Raised for malformed command or invalid arguments."""
    def __init__(self, msg, *args):
        super(ZitherUsageError, self).__init__(msg, *args)


class _ZitherArgumentParser(argparse.ArgumentParser):
    """Argument parser that raises UsageError instead of exiting."""
    #pylint: disable=too-few-public-methods
    def error(self, message):
        '''Suppress default exit behavior'''
        raise ZitherUsageError(message)


class _ExplicitBamFileStrategy(object):
    '''Returns single entry dict of sample name and bam path derived from
    explicit bam file'''
    def __init__(self, bam_file_path):
        self._bam_file_path = bam_file_path

    def build_sample_bam_mapping(self):
        sample_file = os.path.basename(self._bam_file_path)
        sample_name = os.path.splitext(sample_file)[0]
        return {sample_name: self._bam_file_path}


class _MappingFileStrategy(object):
    '''Returns dict of sample names and bam paths derived from
    explicit tab separated mapping file'''
    def __init__(self, mapping_file):
        self._mapping_file = mapping_file

    def _abs_path(self, path):
        mapping_dir_path = os.path.dirname(self._mapping_file)
        if path != os.path.abspath(path):
            path = os.path.abspath(os.path.join(mapping_dir_path, path))
        return path

    def build_sample_bam_mapping(self):
        sample_bam_mapping = OrderedDict()
        with open(self._mapping_file, 'rt') as tsvfile:
            for sample_name, bam_path in csv.reader(tsvfile, delimiter='\t'):
                sample_bam_mapping[sample_name] = self._abs_path(bam_path)
        return sample_bam_mapping

class _MatchingNameStrategy(object):
    '''Returns dict of sample names and bam path derived from
    sample names in VCF'''
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

def _round_digits(val):
    val_str = str(val)
    try:
        if len(val_str.split(".")[1]) > 6:
            return "{0:.6f}".format(val)
    except IndexError:
        return val_str
    return val_str

class _PileupStats(object):
    ''''Calculate basic stats based on dict of counts for A,C,G,T'''
    def __init__(self, ref, alt, total_acgt, filtered_acgt):
        (self.total_depth,
         self.total_depth_acgt,
         self.total_af) = self._init_depth_freq(ref, alt, total_acgt)
        (self.filtered_depth,
         self.filtered_depth_acgt,
         self.filtered_af) = self._init_depth_freq(ref, alt, filtered_acgt)

    @staticmethod
    def _is_snp(ref, alt):
        if len(ref) != 1:
            return False
        if len(alt) == 1:
            return True
        if max(len(x) for x in alt.split(",")) == 1:
            return True
        return False

    @staticmethod
    def _init_depth_freq(ref, alts, acgt):
        alts = alts.upper()
        freq = _NULL
        depth = sum(acgt.values())
        try:
            if depth and _PileupStats._is_snp(ref, alts):
                freqs = []
                for alt in alts.split(","):
                    freqs.append(_round_digits(acgt[alt]/depth))
                freq = ",".join(freqs)
        except KeyError:
            freq = _NULL

        depth_acgt = "{},{},{},{}".format(acgt.get("A", 0),
                                          acgt.get("C", 0),
                                          acgt.get("G", 0),
                                          acgt.get("T", 0))
        return (depth, depth_acgt, freq)

def _basecall_quality_filter(basecall_quality_cutoff):
    def include(read):
        align = read.alignment
        pos = read.query_position
        basecall_qual = align.query_qualities[pos]
        return basecall_qual >= basecall_quality_cutoff
    return include

def _mapping_quality_filter(mapping_quality_cutoff):
    def include(read):
        return read.alignment.mapping_quality >= mapping_quality_cutoff
    return include

def _build_filters(args):
    filters = [_basecall_quality_filter(int(args.basecall_quality_cutoff)),
               _mapping_quality_filter(int(args.mapping_quality_cutoff))]
    def include(read):
        for filter_include in filters:
            if not filter_include(read):
                return False
        return True
    return include

class _BamReader(object):
    '''Reads pileup stats from BAM file. (Assumes BAM index is present.)'''
    def __init__(self,
                 bam_file_name,
                 depth_cutoff,
                 filter_include):
        self._bam_file_name = bam_file_name
        self._depth_cutoff = depth_cutoff + 1
        self._filter_include = filter_include
        #pylint: disable=no-member
        self._bam_file = pysam.AlignmentFile(bam_file_name, "rb")

    def __eq__(self, other):
        return (isinstance(other, _BamReader) and
                self._bam_file_name == other._bam_file_name and
                self._filter_include == other._filter_include and
                self._depth_cutoff == other._depth_cutoff)

    def __hash__(self):
        return hash(self._bam_file_name)

    def get_pileup_stats(self, chrom, pos_one_based, ref, alt):
        total_acgt = defaultdict(int)
        filtered_acgt = defaultdict(int)
        try:
            pileupcolumns = self._bam_file.pileup(chrom,
                                                  pos_one_based-1,
                                                  pos_one_based,
                                                  stepper='nofilter',
                                                  truncate=True,
                                                  max_depth=self._depth_cutoff)
            for pileupcolumn in pileupcolumns:
                for read in pileupcolumn.pileups:
                    align = read.alignment
                    pos = read.query_position
                    if not read.is_del:
                        base = align.query_sequence[pos].upper()
                        total_acgt[base] += 1
#                         basecall_qual = align.query_qualities[pos]
#                         if basecall_qual >= self._basecall_quality_cutoff:
                        if self._filter_include(read):
                            filtered_acgt[base] += 1
        except ValueError as samtools_error:
            if str(samtools_error).startswith("invalid reference"):
                #querying unknown chrom returns depth 0
                pass
            else:
                raise samtools_error

        return _PileupStats(ref, alt, total_acgt, filtered_acgt)

class _Tag(object):
    '''Holds the tag metadata along with a way to extract a value from
    pileup stats'''
    _METAHEADER = '##FORMAT=<ID={},Number={},Type={},Description="{}">'
    #pylint: disable=too-many-arguments
    def __init__(self, vcf_id, number, tag_type, description, stats_method):
        self.metaheader = self._METAHEADER.format(vcf_id,
                                                  number,
                                                  tag_type,
                                                  description)
        self.id = vcf_id
        self._get_value_method = stats_method

    def get_value(self, pileup_stats):
        return str(self._get_value_method(pileup_stats))


total_dp = _Tag("ZTDP", "1", "Integer",
                "Zither total (unfiltered) BAM depth; deletions excluded",
                lambda pileup_stats: pileup_stats.total_depth)
total_af = _Tag("ZTAF", "A", "Float",
                ("Zither total (unfiltered) BAM alt frequency; "
                 "(alt count/total count) rounded to 6 decimals"),
                lambda pileup_stats: pileup_stats.total_af)
filtered_dp = _Tag("ZFDP", "1", "Integer",
                   "Zither filtered BAM depth; deletions excluded",
                   lambda pileup_stats: pileup_stats.filtered_depth)
filtered_af = _Tag("ZFAF", "A", "Float",
                   ("Zither filtered BAM alt frequency; "
                    "(alt count/total count) rounded to 6 decimals"),
                   lambda pileup_stats: pileup_stats.filtered_af)
total_dp_acgt = _Tag("ZTDP_ACGT", "4", "Integer",
                     ("Zither total (unfiltered) BAM depths for A,C,G,T;"
                      " deletions excluded"),
                     lambda pileup_stats: pileup_stats.total_depth_acgt)
filtered_dp_acgt = _Tag("ZFDP_ACGT", "4", "Integer",
                        ("Zither filtered BAM depths for A,C,G,T;"
                         " deletions excluded"),
                        lambda pileup_stats: pileup_stats.filtered_depth_acgt)


DEFAULT_TAGS = [total_dp,
                total_dp_acgt,
                total_af,
                filtered_dp,
                filtered_dp_acgt,
                filtered_af]


def _build_execution_context(argv):
    '''Execution context is included in output VCF for reproducibility'''
    return OrderedDict([("timestamp",
                         datetime.now().strftime('%Y-%m-%d %H:%M:%S')),
                        ("command",
                         ' '.join(argv)),
                        ("cwd",
                         os.getcwd()),
                        ("version",
                         __version__)])

def _get_sample_names(input_vcf):
    '''Returns list of sample names from input VCF'''
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

def _build_reader_dict(sample_bam_mapping, filter_include, args):
    '''Given a sample name to bam path mapping, return dict of sample_name
    to BamReader'''
    readers_dict = OrderedDict()
    for (sample, bam_file) in sample_bam_mapping.items():
        readers_dict[sample] = _BamReader(bam_file,
                                          depth_cutoff=int(args.depth_cutoff),
                                          filter_include=filter_include)
    return readers_dict

def _build_column_header_line(sample_names):
    column_headers = list(_VCF_FIXED_HEADERS)
    column_headers.extend(sample_names)
    return '\t'.join(column_headers)


def _build_vcf_metaheaders(execution_context, tags):
    exec_tags = ['{}="{}"'.format(k, v) for k, v in execution_context.items()]
    zither_metaheader = '##zither=<{}>'.format(",".join(exec_tags))
    vcf_headers = ['##fileformat=VCFv4.1']
    vcf_headers.extend([tag.metaheader for tag in tags])
    vcf_headers.append(zither_metaheader)
    return vcf_headers


def _build_sample_fields(sample_reader_dict, tags, CHROM, POS, REF, ALT):
    sample_fields = []
    for sample_name in sample_reader_dict.keys():
        bam_reader = sample_reader_dict[sample_name]
        pileup_stats = bam_reader.get_pileup_stats(CHROM, int(POS), REF, ALT)
        tag_values = [tag.get_value(pileup_stats) for tag in tags]
        sample_fields.append(':'.join(tag_values))
    return sample_fields


def _create_vcf(input_vcf, sample_reader_dict, execution_context, tags=None):
    '''Reads input VCF and emits output VCF with new Zither tags.'''
    if tags is None:
        tags = DEFAULT_TAGS
    vcf_headers = _build_vcf_metaheaders(execution_context, tags)
    vcf_format = ":".join([tag.id for tag in tags])
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
                vcf_fields.append(vcf_format)
                vcf_fields.extend(_build_sample_fields(sample_reader_dict,
                                                       tags,
                                                       CHROM,
                                                       POS,
                                                       REF,
                                                       ALT))
                print('\t'.join(vcf_fields))

def _parse_command_line_args(arguments):
    parser = _ZitherArgumentParser( \
        formatter_class=argparse.RawTextHelpFormatter,
        usage="zither input_vcf",
        description=(\
'''For all positions in VCF, pull raw depths and alt freqs from BAM file,
writing output as new VCF to stdout.
By default, bams must be in same dir as VCF and bam filenames must match VCF
sample names. See below for alternative approaches (for example, --bam
and --mapping_file).'''))

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
                        default=_DEFAULT_BASECALL_QUALITY_CUTOFF,
                        help="minimum base-call quality to be included. "
                        "Defaults to " + str(_DEFAULT_BASECALL_QUALITY_CUTOFF))
    parser.add_argument('--mapping_quality_cutoff',
                        default=_DEFAULT_MAPPING_QUALITY_CUTOFF,
                        help="minimum mapping quality to be included. "
                        "Defaults to " + str(_DEFAULT_MAPPING_QUALITY_CUTOFF))
    parser.add_argument('--depth_cutoff',
                        default=_DEFAULT_DEPTH_CUTOFF,
                        help="maximum pileup depth for a given position. "
                        "Defaults to " + str(_DEFAULT_DEPTH_CUTOFF))
    args = parser.parse_args(arguments)
    return args

def main(command_line_args=None):
    '''Zither entry point.'''
    try:
        if not command_line_args:
            command_line_args = sys.argv
        args = _parse_command_line_args(command_line_args[1:])
        execution_context = _build_execution_context(command_line_args)
        strategy = _get_sample_bam_strategy(args)
        sample_bam_mapping = strategy.build_sample_bam_mapping()
        filter_include = _build_filters(args)
        reader_dict = _build_reader_dict(sample_bam_mapping,
                                         filter_include,
                                         args)
        _create_vcf(args.input_vcf, reader_dict, execution_context)
    except ZitherUsageError as usage_error:
        message = "Zither usage problem: {}".format(str(usage_error))
        print(message, file=sys.stderr)
        print("See 'zither --help'.", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main(sys.argv)
