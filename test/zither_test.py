from __future__ import print_function, absolute_import, division

from zither import __version__
import zither.zither as zither
import datetime
import filecmp
import os
from testfixtures import TempDirectory
import unittest
import re
from subprocess import call
import sys
import time

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

class _SamFlag(object): 
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

        
def _create_file(path, filename, contents):
    filename = os.path.join(path, filename)
    with open(filename, 'w') as new_file:
        new_file.write(contents)
    return filename

def _create_bam(dir, sam_contents):
    test_input_sam = os.path.join(dir, "test.sam")
    sam_file = open(test_input_sam, "w")
    sam_file.write(sam_contents)
    sam_file.close()
    test_input_bam = os.path.join(dir, "test.bam")
    call(["samtools view -S -b "+ test_input_sam + " -o " + test_input_bam + " 2>/dev/null "], shell=True)
    call(["samtools index " + test_input_bam + " 2>/dev/null "], shell=True)
    return test_input_bam
                    
def _get_zither_metaheader(lines):
    for line in lines:
        if line.startswith("##zither=<"):
            return line
    return None
    
class ReaderTestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.simple_forward_read = _SamFlag.PAIRED + _SamFlag.PROPER_PAIR + _SamFlag.MREVERSE + _SamFlag.READ1
        self.stdout = StringIO()
        self.saved_stdout = sys.stdout
        sys.stdout = self.stdout

    def tearDown(self):
        self.stdout.close()
        sys.stdout = self.saved_stdout
        unittest.TestCase.tearDown(self)

    def assertStartsWith(self, full_text, search_text):
        self.assertTrue(full_text.startswith(search_text), 
                        "[{}] does not start with [{}]".format(full_text, search_text))
        
    def _compare_lines(self, expected, actual):
        expected = expected.split("\n")
        actual = actual.split("\n")
        for i in range(len(expected)):
            if expected[i].startswith("##zither=<timestamp="):
                self.assertStartsWith(actual[i], "##zither=<timestamp=")
            else:
                self.assertEquals(expected[i].rstrip(),
                                  actual[i].rstrip())

    
    def test_get_depth_and_alt_freq(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	5M	=	105	0	AAAAA	>>>>>
readNameB	99	chr10	6	0	5M	=	106	0	CCCCC	>>>>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = zither.Reader(input_bam)
            self.assertEquals((0, '.'), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=4, ref="C", alt="A"))
            self.assertEquals((1, '1.0'), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A"))
            self.assertEquals((2, '0.5'), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=6, ref="C", alt="A"))

    def test_get_depth_and_alt_freq_ignoresDeletions(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	3M	=	105	0	AAA	>>>
readNameB	99	chr10	5	0	1M1D1M	=	105	0	CG	>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = zither.Reader(input_bam)
            self.assertEquals(2, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A")[0])
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=6, ref="C", alt="A")[0])
            self.assertEquals(2, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=7, ref="C", alt="A")[0])

    def test_get_depth_and_alt_freq_ignoresSkippedRef(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	3M	=	105	0	AAA	>>>
readNameB	99	chr10	5	0	1M1N1M	=	105	0	CG	>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = zither.Reader(input_bam)
            self.assertEquals(2, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A")[0])
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=6, ref="C", alt="A")[0])
            self.assertEquals(2, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=7, ref="C", alt="A")[0])

    def test_get_depth_and_alt_freq_ignoresCase(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	3M	=	105	0	AAA	>>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = zither.Reader(input_bam)
            self.assertEquals((1,'1.0'), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="a"))
            self.assertEquals((1,'1.0'), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A"))

    def test_get_depth_and_alt_freq_skipsFreqForIndels(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	3M	=	105	0	AAA	>>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = zither.Reader(input_bam)
            self.assertEquals((1,'.'), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="A", alt="ATA"))
            self.assertEquals((1,'.'), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="ATA", alt="A"))
 
    def test_get_depth_and_alt_freq_ignoresInsertions(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	8M	=	105	0	ACGTACGT	>>>>>>>>
readNameA	99	chr10	5	0	4M2I4M	=	105	0	ACGTTTACGT	>>>>>>>>>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = zither.Reader(input_bam)
            self.assertEquals((2,'1.0'), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A"))
            self.assertEquals((2,'1.0'), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=9, ref="C", alt="A"))

    def test_get_depth_and_alt_freq_overlappingReadPairsCountedSeparately(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	4M	=	105	0	AAAA	>>>>
readNameA	147	chr10	7	0	4M	=	105	0	TTTT	>>>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = zither.Reader(input_bam)
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A")[0])
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=6, ref="C", alt="A")[0])
            self.assertEquals(2, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=7, ref="C", alt="A")[0])
            self.assertEquals(2, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=8, ref="C", alt="A")[0])
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=9, ref="C", alt="A")[0])
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=10, ref="C", alt="A")[0])


    def test_get_depth_and_alt_freq_ignoresBaseCallQuality(self):
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
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = zither.Reader(input_bam)
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A")[0])
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=6, ref="C", alt="A")[0])
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=7, ref="C", alt="A")[0])
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=8, ref="C", alt="A")[0])
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=9, ref="C", alt="A")[0])

    def test_get_depth_and_alt_freq_duplicatesCounted(self):
        
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	{}	chr10	5	0	5M	=	105	0	AAAAA	>>>>>
'''.format(self.simple_forward_read + _SamFlag.DUP)
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = zither.Reader(input_bam)
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A")[0])
            
    def test_get_depth_and_alt_freq_failedQualityCounted(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	{}	chr10	5	0	5M	=	105	0	AAAAA	>>>>>
'''.format(self.simple_forward_read + _SamFlag.QCFAIL)
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = zither.Reader(input_bam)
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A")[0])
            
    def test_get_depth_and_alt_freq_secondaryReadsCounted(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	{}	chr10	5	0	5M	=	105	0	AAAAA	>>>>>
'''.format(self.simple_forward_read + _SamFlag.SECONDARY)
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = zither.Reader(input_bam)
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A")[0])

    def test_create_vcf(self):
        input_vcf_contents = \
'''##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_A
chr1	10	.	A	C	.	.	.	GT	0/1
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
            input_bam = _create_bam(tmp_path, sam_contents)
            reader = zither.Reader(input_bam)
            
            reader.create_vcf(input_vcf)
        
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
            input_bam = _create_bam(tmp_path, sam_contents)
            reader = zither.Reader(input_bam)
            
            reader.create_vcf(input_vcf)
        
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
            input_bam = _create_bam(tmp_path, sam_contents)
            reader = zither.Reader(input_bam)
            
            reader.create_vcf(input_vcf)
        
            actual_output_lines = self.stdout.getvalue()
            
            self._compare_lines(expected_vcf_contents, actual_output_lines)

            
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
            input_bam = _create_bam(tmp_path, sam_contents)
            reader = zither.Reader(input_bam)

            reader.create_vcf(input_vcf)
        
            actual_output_lines = self.stdout.getvalue().split("\n")
            
            tag_contents = re.search("##zither=<(.*)>", _get_zither_metaheader(actual_output_lines)).group(1)
            tags = tag_contents.split(",")
            self.assertStartsWith(tags[0], 'timestamp="')
            actual_timestamp = datetime.datetime.strptime(tags[0].split("=")[1], '"%Y-%m-%d %H:%M:%S"')
            self.assertAlmostEqual(time.mktime(datetime.datetime.now().timetuple()),
                                   time.mktime(actual_timestamp.timetuple()))
            self.assertEquals('command="/usr/bin/nosetests"', tags[1])
            self.assertStartsWith(tags[2], 'cwd="')
            self.assertEquals('version="'+__version__+'"', tags[3])
     