#pylint: disable=invalid-name, too-few-public-methods
from __future__ import print_function, absolute_import, division
from argparse import Namespace
from collections import OrderedDict
import datetime
import os
import pysam
import sys
from testfixtures.tempdirectory import TempDirectory
import time
import unittest
import zither.zither as zither

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

# I would rather just say pysam.view(...), but since that global is
# added dynamically, Eclipse flags this as a compilation problem. So
# instead we connect directly to the pysam.SamtoolsDispatcher.
PYSAM_VIEW = pysam.SamtoolsDispatcher("view", None).__call__
PYSAM_INDEX = pysam.SamtoolsDispatcher("index", None).__call__

def _absolute_base_dir(*path_components):
    base_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    return os.path.join(base_path, *path_components)

def _create_file(path, filename, contents):
    filename = os.path.join(path, filename)
    with open(filename, 'w') as new_file:
        new_file.write(contents)
    return filename

def _create_bam(path, filename, sam_contents):
    sam_filename = _create_file(path, filename, sam_contents)
    bam_filename = sam_filename.replace(".sam", ".bam")
    _pysam_bam_from_sam(sam_filename, bam_filename)
    return bam_filename

def _pysam_bam_from_sam(sam_filename, bam_filename):
    temp_stdout = sys.stdout
    try:
        sys.stdout = sys.__stdout__
        bam = PYSAM_VIEW("-S", "-b", sam_filename)
        bam_file = open(bam_filename, "wb")
        bam_file.writelines(bam)
        bam_file.close()
        PYSAM_INDEX(bam_filename)

    finally:
        sys.stdout = temp_stdout

def _get_zither_metaheader(lines):
    return _get_line_starts_with(lines, "##zither=<")

def _get_column_header(lines):
    return _get_line_starts_with(lines, "#CHROM")

def _get_line_starts_with(lines, starts_with):
    for line in lines:
        if line.startswith(starts_with):
            return line
    return None


class _BamFlag(object):
    #the read is paired in sequencing, no matter whether it is mapped in a pair
    PAIRED = 1
    # the read is mapped in a proper pair
    PROPER_PAIR = 2
    #the read itself is unmapped; conflictive with PROPER_PAIR
    UNMAP = 4
    #the mate is unmapped
    MUNMAP = 8
    #the read is mapped to the reverse strand
    REVERSE = 16
    #the mate is mapped to the reverse strand
    MREVERSE = 32
    #this is read1
    READ1 = 64
    #this is read2
    READ2 = 128
    #not primary alignment
    SECONDARY = 256
    #QC failure
    QCFAIL = 512
    #optical or PCR duplicate
    DUP = 1024
    #supplementary alignment
    SUPPLEMENTARY = 2048


class ZitherBaseTestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.simple_forward_read = (_BamFlag.PAIRED
                                    + _BamFlag.PROPER_PAIR
                                    + _BamFlag.MREVERSE
                                    + _BamFlag.READ1)
        self.stdout = StringIO()
        self.saved_stdout = sys.stdout
        sys.stdout = self.stdout

    def tearDown(self):
        self.stdout.close()
        sys.stdout = self.saved_stdout
        unittest.TestCase.tearDown(self)

    def assertStartsWith(self, full_text, search_text):
        self.assertTrue(full_text.startswith(search_text),
                        "[{}] does not start with [{}]".format(full_text,
                                                               search_text))

    def _compare_lines(self, expected, actual):
        expected = expected.split("\n")
        actual = actual.split("\n")
        for i in range(len(expected)):
            if expected[i].startswith("##zither=<"):
                self.assertStartsWith(actual[i], "##zither=<")
            else:
                self.assertEquals(expected[i].rstrip(),
                                  actual[i].rstrip())


class ExplicitBamFileStrategyTestCase(ZitherBaseTestCase):
    def test_build_sample_bam_mapping(self):
        strategy = zither._ExplicitBamFileStrategy("/foo/bar/baz.bam")
        actual_mapping = strategy.build_sample_bam_mapping()
        self.assertEquals(["baz"],  list(actual_mapping.keys()))
        self.assertEquals("/foo/bar/baz.bam",  actual_mapping["baz"])

    def test_build_sample_bam_mapping_longExtension(self):
        strategy = zither._ExplicitBamFileStrategy("/foo/bar/baz.hoopy.oop.bam")
        actual_mapping = strategy.build_sample_bam_mapping()
        self.assertEquals(["baz.hoopy.oop"],  list(actual_mapping.keys()))
        self.assertEquals("/foo/bar/baz.hoopy.oop.bam",
                          actual_mapping["baz.hoopy.oop"])

    def test_build_sample_bam_mapping_noExtension(self):
        strategy = zither._ExplicitBamFileStrategy("/foo/bar/baz")
        actual_mapping = strategy.build_sample_bam_mapping()
        self.assertEquals(["baz"],  list(actual_mapping.keys()))
        self.assertEquals("/foo/bar/baz",  actual_mapping["baz"])


class MatchingNameStrategyTestCase(ZitherBaseTestCase):
    def test_build_sample_bam_mapping(self):
        strategy = zither._MatchingNameStrategy(["sA", "sB"],
                                                "/foo/bar/input.vcf")
        actual_mapping = strategy.build_sample_bam_mapping()
        self.assertEquals(["sA", "sB"],  sorted(actual_mapping.keys()))
        self.assertEquals("/foo/bar/sA.bam",  actual_mapping["sA"])
        self.assertEquals("/foo/bar/sB.bam",  actual_mapping["sB"])

    def test_build_sample_bam_mapping_emptySampleList(self):
        strategy = zither._MatchingNameStrategy([],
                                                "/foo/bar/input.vcf")
        actual_mapping = strategy.build_sample_bam_mapping()
        self.assertEquals([], list(actual_mapping.keys()))


class MappingFileStrategyTestCase(ZitherBaseTestCase):
    def test_build_sample_bam_mapping(self):
        mapping_file_contents = \
'''sA	/foo/bar/sA.bam
sB	/foo/bar/sB.bam
'''
        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            mapping_file = _create_file(tmp_path,
                                        "mapping_file.txt",
                                        mapping_file_contents)

            strategy = zither._MappingFileStrategy(mapping_file)
            actual_mapping = strategy.build_sample_bam_mapping()
            self.assertEquals(["sA", "sB"],  sorted(actual_mapping.keys()))
            self.assertEquals("/foo/bar/sA.bam",  actual_mapping["sA"])
            self.assertEquals("/foo/bar/sB.bam",  actual_mapping["sB"])

    def test_build_sample_bam_mapping_preservesSampleOrderFromMappingFile(self):
        mapping_file_contents = \
'''sX	/foo
s1a	/foo
'''
        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            mapping_file = _create_file(tmp_path,
                                        "mapping_file.txt",
                                        mapping_file_contents)

            strategy = zither._MappingFileStrategy(mapping_file)
            actual_mapping = strategy.build_sample_bam_mapping()
            self.assertEquals(["sX", "s1a"],  list(actual_mapping.keys()))

    def test_build_sample_bam_mapping_invalidMappingFile(self):
        mapping_file_contents = \
'''sA
sB	/foo/bar/sB.bam
'''
        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            mapping_file = _create_file(tmp_path,
                                        "mapping_file.txt",
                                        mapping_file_contents)

            strategy = zither._MappingFileStrategy(mapping_file)

            self.assertRaises(ValueError,
                strategy.build_sample_bam_mapping)

class PileupStats(ZitherBaseTestCase):
    def test_init(self):
        coverage = [[1], [2], [4], [8]]
        stats = zither._PileupStats("A", "C", coverage)
        self.assertEquals(15, stats.unfiltered_depth)
        self.assertEquals(str(2/15), stats.unfiltered_af)
    
    def test_no_depth(self):
        coverage = [[0], [0], [0], [0]]
        stats = zither._PileupStats("A", "C", coverage)
        self.assertEquals(0, stats.unfiltered_depth)
        self.assertEquals(".", stats.unfiltered_af)

    def test_insertion(self):
        coverage = [[1], [2], [4], [8]]
        stats = zither._PileupStats("A", "ACGT", coverage)
        self.assertEquals(15, stats.unfiltered_depth)
        self.assertEquals(".", stats.unfiltered_af)

    def test_deletion(self):
        coverage = [[1], [2], [4], [8]]
        stats = zither._PileupStats("ACGT", "A", coverage)
        self.assertEquals(15, stats.unfiltered_depth)
        self.assertEquals(".", stats.unfiltered_af)

    def test_multalt(self):
        coverage = [[1], [2], [4], [8]]
        stats = zither._PileupStats("A", "C,G", coverage)
        self.assertEquals(15, stats.unfiltered_depth)
        self.assertEquals(".", stats.unfiltered_af)

        
class BamReaderTestCase(ZitherBaseTestCase):
    def test_equals(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr1	LN:10
readA	99	chr1	10	0	1M	=	105	0	A	>
'''
        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            bam_a = _create_bam(tmp_path, "sample_A.sam", sam_contents)
            bam_b = _create_bam(tmp_path, "sample_B.sam", sam_contents)

            base = zither._BamReader(bam_a)
            base_equivalent = zither._BamReader(bam_a)
            self.assertEquals(base, base_equivalent)
            different_file = zither._BamReader(bam_b)
            self.assertNotEquals(base, different_file)

    def test_hash(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr1	LN:10
readA	99	chr1	10	0	1M	=	105	0	A	>
'''
        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            bam_A = _create_bam(tmp_path,
                                "sample_A.sam",
                                sam_contents)

            base = zither._BamReader(bam_A)
            base_equivalent = zither._BamReader(bam_A)
            self.assertEquals(base.__hash__(), base_equivalent.__hash__())

            reader_set = set()
            reader_set.add(base)
            reader_set.add(base_equivalent)
            self.assertEquals(1, len(reader_set))


    def test_get_pileup_stats(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	5M	=	105	0	AAAAA	>>>>>
readNameB	99	chr10	6	0	5M	=	106	0	CCCCC	>>>>>
'''
        with TempDirectory() as tmp_dir:
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                         pos_one_based=4,
                                                         ref="C",
                                                         alt="A")
            self.assertEquals(0, actual_stats.unfiltered_depth)
            self.assertEquals('.', actual_stats.unfiltered_af)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                         pos_one_based=5,
                                                         ref="C",
                                                         alt="A")
            self.assertEquals(1, actual_stats.unfiltered_depth)
            self.assertEquals('1.0', actual_stats.unfiltered_af)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                         pos_one_based=6,
                                                         ref="C",
                                                         alt="A")
            self.assertEquals(2, actual_stats.unfiltered_depth)
            self.assertEquals('0.5', actual_stats.unfiltered_af)

            
    def test_get_pileup_stats_ignoresDeletions(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	3M	=	105	0	AAA	>>>
readNameB	99	chr10	5	0	1M1D1M	=	105	0	CG	>>
'''
        with TempDirectory() as tmp_dir:
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)
            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=5,
                                                   ref="C",
                                                   alt="A")
            self.assertEquals(2, actual_stats.unfiltered_depth)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=6,
                                                   ref="C",
                                                   alt="A")
            self.assertEquals(1, actual_stats.unfiltered_depth)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=7,
                                                   ref="C",
                                                   alt="A")
            self.assertEquals(2, actual_stats.unfiltered_depth)

    def test_get_pileup_stats_ignoresSkippedRef(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	3M	=	105	0	AAA	>>>
readNameB	99	chr10	5	0	1M1N1M	=	105	0	CG	>>
'''
        with TempDirectory() as tmp_dir:
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)
            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=5,
                                                   ref="C",
                                                   alt="A")
            self.assertEquals(2, actual_stats.unfiltered_depth)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=6,
                                                   ref="C",
                                                   alt="A")
            self.assertEquals(1, actual_stats.unfiltered_depth)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=7,
                                                   ref="C",
                                                   alt="A")
            self.assertEquals(2, actual_stats.unfiltered_depth)

    def test_get_pileup_stats_ignoresCase(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	3M	=	105	0	AAA	>>>
'''
        with TempDirectory() as tmp_dir:
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)
            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=5,
                                                   ref="C",
                                                   alt="a")
            self.assertEquals(1, actual_stats.unfiltered_depth)
            self.assertEquals('1.0', actual_stats.unfiltered_af)
            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=5,
                                                   ref="C",
                                                   alt="A")
            self.assertEquals(1, actual_stats.unfiltered_depth)
            self.assertEquals('1.0', actual_stats.unfiltered_af)

    def test_get_pileup_stats_skipsFreqForIndels(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	3M	=	105	0	AAA	>>>
'''
        with TempDirectory() as tmp_dir:
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=5,
                                                   ref="A",
                                                   alt="ATA")
            self.assertEquals(1, actual_stats.unfiltered_depth)
            self.assertEquals('.', actual_stats.unfiltered_af)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=5,
                                                   ref="ATA",
                                                   alt="A")
            self.assertEquals(1, actual_stats.unfiltered_depth)
            self.assertEquals('.', actual_stats.unfiltered_af)

    def test_get_pileup_stats_ignoresInsertions(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	8M	=	105	0	ACGTACGT	>>>>>>>>
readNameA	99	chr10	5	0	4M2I4M	=	105	0	ACGTTTACGT	>>>>>>>>>>
'''
        with TempDirectory() as tmp_dir:
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=5,
                                                   ref="C",
                                                   alt="A")
            self.assertEquals(2, actual_stats.unfiltered_depth)
            self.assertEquals('1.0', actual_stats.unfiltered_af)
            
            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=9,
                                                   ref="C",
                                                   alt="A")
            self.assertEquals(2, actual_stats.unfiltered_depth)
            self.assertEquals('1.0', actual_stats.unfiltered_af)

    def test_get_pileup_stats_overlappingReadPairsCountedSeparately(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	4M	=	105	0	AAAA	>>>>
readNameA	147	chr10	7	0	4M	=	105	0	TTTT	>>>>
'''
        with TempDirectory() as tmp_dir:
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                            pos_one_based=5,
                                                            ref="C",
                                                            alt="A")

            self.assertEquals(1, actual_stats.unfiltered_depth)
            
            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                            pos_one_based=6,
                                                            ref="C",
                                                            alt="A")
            self.assertEquals(1, actual_stats.unfiltered_depth)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                            pos_one_based=7,
                                                            ref="C",
                                                            alt="A")

            self.assertEquals(2, actual_stats.unfiltered_depth)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                            pos_one_based=8,
                                                            ref="C",
                                                            alt="A")

            self.assertEquals(2, actual_stats.unfiltered_depth)
            
            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                            pos_one_based=9,
                                                            ref="C",
                                                            alt="A")
            self.assertEquals(1, actual_stats.unfiltered_depth)


            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                            pos_one_based=10,
                                                            ref="C",
                                                            alt="A")

            self.assertEquals(1, actual_stats.unfiltered_depth)


    def test_get_pileup_stats_ignoresBaseCallQuality(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	5M	=	105	0	AAAAA	50+&!
'''
        #Q  ASCII
        #0  !
        #5  &
        #10 +
        #15 0
        #20 5
        with TempDirectory() as tmp_dir:
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=5,
                                                   ref="C",
                                                   alt="A")

            self.assertEquals(1, actual_stats.unfiltered_depth)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=6,
                                                   ref="C", alt="A")

            self.assertEquals(1, actual_stats.unfiltered_depth)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=7,
                                                   ref="C",
                                                   alt="A")

            self.assertEquals(1, actual_stats.unfiltered_depth)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=8,
                                                   ref="C",
                                                   alt="A")

            self.assertEquals(1, actual_stats.unfiltered_depth)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=9,
                                                   ref="C",
                                                   alt="A")

            self.assertEquals(1, actual_stats.unfiltered_depth)

    def test_get_pileup_stats_duplicatesCounted(self):

        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	{}	chr10	5	0	5M	=	105	0	AAAAA	>>>>>
'''.format(self.simple_forward_read + _BamFlag.DUP)
        with TempDirectory() as tmp_dir:
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=5,
                                                   ref="C",
                                                   alt="A")
            self.assertEquals(1, actual_stats.unfiltered_depth)

    def test_get_pileup_stats_failedQualityCounted(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	{}	chr10	5	0	5M	=	105	0	AAAAA	>>>>>
'''.format(self.simple_forward_read + _BamFlag.QCFAIL)
        with TempDirectory() as tmp_dir:
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=5,
                                                   ref="C",
                                                   alt="A")
            self.assertEquals(1, actual_stats.unfiltered_depth)

    def test_get_pileup_stats_secondaryReadsCounted(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	{}	chr10	5	0	5M	=	105	0	AAAAA	>>>>>
'''.format(self.simple_forward_read + _BamFlag.SECONDARY)
        with TempDirectory() as tmp_dir:
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)

            actual_stats = reader.get_pileup_stats(chrom="chr10",
                                                   pos_one_based=5,
                                                   ref="C",
                                                   alt="A")
            self.assertEquals(1, actual_stats.unfiltered_depth)

class ZitherTestCase(ZitherBaseTestCase):
    def test_build_execution_context(self):
        argv = ["zither", "foo", "bar", "baz"]
        actual_text = zither._build_execution_context(argv)

        self.assertEquals([ "timestamp", "command", "cwd", "version"],
                          list(actual_text.keys()))

        actual_timestamp = datetime.datetime.strptime(actual_text["timestamp"],
                                                      '%Y-%m-%d %H:%M:%S')
        now_in_seconds = time.mktime(datetime.datetime.now().timetuple())
        timestamp_in_seconds = time.mktime(actual_timestamp.timetuple())
        self.assertTrue((now_in_seconds - timestamp_in_seconds) < 2)
        self.assertEquals("zither foo bar baz", actual_text["command"])
        self.assertEquals(os.getcwd(), actual_text["cwd"])
        self.assertEquals(zither.__version__, actual_text["version"])

    def test_create_vcf_singleSample(self):
        input_vcf_contents = \
'''##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_A
chr1	10	.	A	C	QAULStuff	FILTERStuff	INFOStuff	GT	0/1
chr10	10	.	C	G	.	.	.	GT	0/1
chr15	42	.	G	T	.	.	.	GT	0/1
'''

        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr1	LN:10
@SQ	SN:chr10	LN:10
@SQ	SN:chr15	LN:10
readA	99	chr1	10	0	1M	=	105	0	C	>
readB	99	chr1	10	0	1M	=	105	0	C	>
readC	99	chr10	10	0	1M	=	105	0	C	>
readD	99	chr10	10	0	1M	=	105	0	G	>
readE	99	chr15	42	0	1M	=	105	0	T	>
readF	99	chr15	42	0	1M	=	105	0	T	>
readG	99	chr15	42	0	1M	=	105	0	T	>
readH	99	chr15	42	0	1M	=	105	0	G	>
'''

        expected_vcf_contents = \
'''##fileformat=VCFv4.1
##FORMAT=<ID=BDP,Number=1,Type=Integer,Description="BAM depth">
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="BAM alt frequency">
##zither=<timestamp=...
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_A
chr1	10	.	A	C	.	.	.	BDP:BAF	2:1.0
chr10	10	.	C	G	.	.	.	BDP:BAF	2:0.5
chr15	42	.	G	T	.	.	.	BDP:BAF	4:0.75
'''

        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            input_vcf = _create_file(tmp_path, "input.vcf", input_vcf_contents)
            bam_A = _create_bam(tmp_path, "sample_A.sam", sam_contents)

            sample_reader_dict = {"sample_A": zither._BamReader(bam_A)}
            zither._create_vcf(input_vcf, sample_reader_dict, {})

            actual_output_lines = self.stdout.getvalue()

            self._compare_lines(expected_vcf_contents, actual_output_lines)


    def test_create_vcf_multipleSamples(self):
        input_vcf_contents = \
'''##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_A	sample_B
chr1	10	.	A	C	QAULStuff	FILTERStuff	INFOStuff	GT	0/1	0/1
chr10	10	.	C	G	.	.	.	GT	0/1	0/1
'''

        sam_contents_A = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr1	LN:10
@SQ	SN:chr10	LN:10
readA	99	chr1	10	0	1M	=	105	0	A	>
readB	99	chr1	10	0	1M	=	105	0	A	>
readC	99	chr1	10	0	1M	=	105	0	A	>
readD	99	chr10	10	0	1M	=	105	0	C	>
readE	99	chr10	10	0	1M	=	105	0	C	>
readF	99	chr10	10	0	1M	=	105	0	C	>
'''

        sam_contents_B = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr1	LN:10
@SQ	SN:chr10	LN:10
readA	99	chr1	10	0	1M	=	105	0	C	>
readB	99	chr1	10	0	1M	=	105	0	C	>
readC	99	chr10	10	0	1M	=	105	0	C	>
readD	99	chr10	10	0	1M	=	105	0	G	>
'''


        expected_vcf_contents = \
'''##fileformat=VCFv4.1
##FORMAT=<ID=BDP,Number=1,Type=Integer,Description="BAM depth">
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="BAM alt frequency">
##zither=<timestamp=...
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_A	sample_B
chr1	10	.	A	C	.	.	.	BDP:BAF	3:0.0	2:1.0
chr10	10	.	C	G	.	.	.	BDP:BAF	3:0.0	2:0.5
'''

        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            input_vcf = _create_file(tmp_path, "input.vcf", input_vcf_contents)
            bam_A = _create_bam(tmp_path, "sample_A.sam", sam_contents_A)
            bam_B = _create_bam(tmp_path, "sample_B.sam", sam_contents_B)

            sample_reader_dict = OrderedDict([
                                    ("sample_A", zither._BamReader(bam_A)),
                                    ("sample_B", zither._BamReader(bam_B))
                                    ])
            zither._create_vcf(input_vcf, sample_reader_dict, {})

            actual_output_lines = self.stdout.getvalue()

            self._compare_lines(expected_vcf_contents, actual_output_lines)



    def test_create_vcf_missingBAMLocationsOk(self):
        input_vcf_contents = \
'''##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_A
chr1	10	.	A	C	.	.	.	GT	0/1
chr10	10	.	C	G	.	.	.	GT	0/1
'''

        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr1	LN:10
@SQ	SN:chr10	LN:10
@SQ	SN:chr15	LN:10
readA	99	chr1	10	0	1M	=	105	0	C	>
readB	99	chr1	10	0	1M	=	105	0	C	>
'''

        expected_vcf_contents = \
'''##fileformat=VCFv4.1
##FORMAT=<ID=BDP,Number=1,Type=Integer,Description="BAM depth">
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="BAM alt frequency">
##zither=<timestamp=...
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_A
chr1	10	.	A	C	.	.	.	BDP:BAF	2:1.0
chr10	10	.	C	G	.	.	.	BDP:BAF	0:.
'''

        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            input_vcf = _create_file(tmp_path, "input.vcf", input_vcf_contents)
            bam_A = _create_bam(tmp_path, "sample_A.sam", sam_contents)

            sample_reader_dict = {"sample_A": zither._BamReader(bam_A)}
            zither._create_vcf(input_vcf, sample_reader_dict, {})

            actual_output_lines = self.stdout.getvalue()
            self._compare_lines(expected_vcf_contents, actual_output_lines)

    def test_create_vcf_multAltReturnsNullAltFreq(self):
        input_vcf_contents = \
'''##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_A
chr1	10	.	A	C,G	.	.	.	GT	0/1
'''

        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr1	LN:10
@SQ	SN:chr10	LN:10
@SQ	SN:chr15	LN:10
readA	99	chr1	10	0	1M	=	105	0	A	>
readB	99	chr1	10	0	1M	=	105	0	C	>
readC	99	chr1	10	0	1M	=	105	0	G	>
'''

        expected_vcf_contents = \
'''##fileformat=VCFv4.1
##FORMAT=<ID=BDP,Number=1,Type=Integer,Description="BAM depth">
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="BAM alt frequency">
##zither=<timestamp=...
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_A
chr1	10	.	A	C,G	.	.	.	BDP:BAF	3:.
'''

        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            input_vcf = _create_file(tmp_path, "input.vcf", input_vcf_contents)
            bam_A = _create_bam(tmp_path, "sample_A.sam", sam_contents)

            sample_reader_dict = {"sample_A": zither._BamReader(bam_A)}
            zither._create_vcf(input_vcf, sample_reader_dict, {})


            actual_output_lines = self.stdout.getvalue()

            self._compare_lines(expected_vcf_contents, actual_output_lines)


    def test_create_vcf_pullsSampleNamesFromSampleReaderDict(self):
        input_vcf_contents = \
'''##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_A
chr1	10	.	A	C	.	.	.	GT	0/1
'''

        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr1	LN:10
readA	99	chr1	10	0	1M	=	105	0	A	>
'''

        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            input_vcf = _create_file(tmp_path, "input.vcf", input_vcf_contents)
            bam_A = _create_bam(tmp_path, "sample_A.sam", sam_contents)

            sample_reader_dict = {"mySampleName": zither._BamReader(bam_A)}
            zither._create_vcf(input_vcf, sample_reader_dict, {})

        actual_output_lines = self.stdout.getvalue().split("\n")
        column_fields = _get_column_header(actual_output_lines).split("\t")

        self.assertEquals(10, len(column_fields))
        self.assertEquals("mySampleName", column_fields[9])

    def test_create_vcf_addsCommandLineMetaHeader(self):
        input_vcf_contents = \
'''##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_A
chr1	10	.	A	C	.	.	.	GT	0/1
'''

        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr1	LN:10
readA	99	chr1	10	0	1M	=	105	0	A	>
'''

        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            input_vcf = _create_file(tmp_path, "input.vcf", input_vcf_contents)
            bam_A = _create_bam(tmp_path, "sample_A.sam", sam_contents)
            sample_reader_dict = {"sample_A": zither._BamReader(bam_A)}
            execution_context=OrderedDict([("foo", "A"),
                                           ("bar", "B"),
                                           ("baz", "C")])

            zither._create_vcf(input_vcf, sample_reader_dict, execution_context)

            actual_output_lines = self.stdout.getvalue().split("\n")

            self.assertEquals('##zither=<foo="A",bar="B",baz="C">',
                              _get_zither_metaheader(actual_output_lines))


    def test_build_reader_dict(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr1	LN:10
readA	99	chr1	10	0	1M	=	105	0	A	>
'''
        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            bam_A = _create_bam(tmp_path, "sample_A.sam", sam_contents)
            bam_B = _create_bam(tmp_path, "sample_B.sam", sam_contents)

            sample_bam_mapping = {'sA':bam_A, 'sB':bam_B}
            actual_mapping = zither._build_reader_dict(sample_bam_mapping)
            self.assertEquals(["sA", "sB"],
                              sorted(actual_mapping.keys()))
            self.assertEquals(zither._BamReader(bam_A),
                              actual_mapping["sA"])
            self.assertEquals(zither._BamReader(bam_B),
                              actual_mapping["sB"])

    def test_build_reader_dict_preservesSampleOrder(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr1	LN:10
readA	99	chr1	10	0	1M	=	105	0	A	>
'''
        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            bam_A = _create_bam(tmp_path, "sample_X.sam", sam_contents)
            bam_B = _create_bam(tmp_path, "sample_1.sam", sam_contents)

            sample_bam_mapping = OrderedDict([('sX', bam_A),('s1', bam_B)])
            actual_mapping = zither._build_reader_dict(sample_bam_mapping)
            self.assertEquals(["sX", "s1"],
                              list(actual_mapping.keys()))


    def test_get_sample_bam_strategy(self):
        input_vcf_contents = \
'''##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_A
chr1	10	.	A	C	.	.	.	GT	0/1
'''
        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            vcf_filename = _create_file(tmp_path, "foo.vcf", input_vcf_contents)

            args = Namespace(input_vcf=vcf_filename,
                             mapping_file=None,
                             bam=None)
            actual_strategy = zither._get_sample_bam_strategy(args)
            self.assertIsInstance(actual_strategy, zither._MatchingNameStrategy)

            args = Namespace(input_vcf=vcf_filename,
                             mapping_file="foo.txt",
                             bam=None)
            actual_strategy = zither._get_sample_bam_strategy(args)
            self.assertIsInstance(actual_strategy, zither._MappingFileStrategy)

            args = Namespace(input_vcf=vcf_filename,
                             mapping_file=None,
                             bam="/foo/bar/baz.bam")
            actual_strategy = zither._get_sample_bam_strategy(args)
            self.assertIsInstance(actual_strategy,
                                  zither._ExplicitBamFileStrategy)


class ZitherFunctionalTestCase(ZitherBaseTestCase):
    def test_main_explicit_bam(self):
        args = ["zither",
                _absolute_base_dir("examples",
                                   "explicit_bam",
                                   "input.vcf"),
                "--bam",
                _absolute_base_dir("examples", "explicit_bam", "sample_X.bam")]
        expected_vcf_filename = _absolute_base_dir("examples",
                                                   "explicit_bam",
                                                   "expected_output.vcf")
        with open(expected_vcf_filename, 'r') as expected_vcf:
            expected_vcf_contents = expected_vcf.read()
            zither.main(args)
            actual_output_lines = self.stdout.getvalue()
            self._compare_lines(expected_vcf_contents, actual_output_lines)

    def test_main_matching_names(self):
        args = ["zither",
                _absolute_base_dir("examples",
                                   "matching_names",
                                   "input.vcf")]
        expected_vcf_filename = _absolute_base_dir("examples",
                                                   "matching_names",
                                                   "expected_output.vcf")
        with open(expected_vcf_filename, 'r') as expected_vcf:
            expected_vcf_contents = expected_vcf.read()
            zither.main(args)
            actual_output_lines = self.stdout.getvalue()
            self._compare_lines(expected_vcf_contents, actual_output_lines)


    def test_main_mapping_file(self):
        args = ["zither",
                _absolute_base_dir("examples", "mapping_files", "input.vcf") ,
                "--mapping_file",
                _absolute_base_dir("examples",
                                   "mapping_files",
                                   "mapping_file.txt")]
        expected_vcf_filename = _absolute_base_dir("examples",
                                                   "mapping_files",
                                                   "expected_output.vcf")
        with open(expected_vcf_filename, 'r') as expected_vcf:
            expected_vcf_contents = expected_vcf.read()
            zither.main(args)
            actual_output_lines = self.stdout.getvalue()
            self._compare_lines(expected_vcf_contents, actual_output_lines)
