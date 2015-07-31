from __future__ import print_function, absolute_import, division

from argparse import Namespace
import datetime
import filecmp
import os
import re
from subprocess import call
import sys
from testfixtures import TempDirectory
import time
import unittest
from zither import __version__
import zither.zither as zither

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

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

        
def _create_file(path, filename, contents):
    filename = os.path.join(path, filename)
    with open(filename, 'w') as new_file:
        new_file.write(contents)
    return filename

#def _create_bam(dir, sam_contents, filename="test.sam"):
def _create_bam(dir, filename, sam_contents):
    test_input_sam = os.path.join(dir, filename)
    bam_filename = test_input_sam.replace(".sam", ".bam")
    sam_file = open(test_input_sam, "w")
    sam_file.write(sam_contents)
    sam_file.close()
    call(["samtools view -S -b "+ test_input_sam + " -o " + bam_filename + " 2>/dev/null "], shell=True)
    call(["samtools index " + bam_filename + " 2>/dev/null "], shell=True)
    return bam_filename


def _get_zither_metaheader(lines):
    for line in lines:
        if line.startswith("##zither=<"):
            return line
    return None
    
class ZitherBaseTestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.simple_forward_read = _BamFlag.PAIRED + _BamFlag.PROPER_PAIR + _BamFlag.MREVERSE + _BamFlag.READ1
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
    
    
class BamReaderTestCase(ZitherBaseTestCase):
    def test_equals(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr1	LN:10
readA	99	chr1	10	0	1M	=	105	0	A	>
'''
        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            bam_A = _create_bam(tmp_path, "sample_A.sam", sam_contents)
            bam_B = _create_bam(tmp_path, "sample_B.sam", sam_contents)

            base = zither._BamReader(bam_A)
            base_equivalent = zither._BamReader(bam_A)
            self.assertEquals(base, base_equivalent)
            different_file = zither._BamReader(bam_B)
            self.assertNotEquals(base, different_file)
 
    def testHash(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr1	LN:10
readA	99	chr1	10	0	1M	=	105	0	A	>
'''
        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            bam_A = _create_bam(tmp_path, "sample_A.sam", sam_contents)

            base = zither._BamReader(bam_A)
            base_equivalent = zither._BamReader(bam_A)
            self.assertEquals(base.__hash__(), base_equivalent.__hash__())
            
            reader_set = set()
            reader_set.add(base)
            reader_set.add(base_equivalent)
            self.assertEquals(1, len(reader_set))


    def test_get_depth_and_alt_freq(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	5M	=	105	0	AAAAA	>>>>>
readNameB	99	chr10	6	0	5M	=	106	0	CCCCC	>>>>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)
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
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)
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
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)
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
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)
            self.assertEquals((1,'1.0'), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="a"))
            self.assertEquals((1,'1.0'), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A"))

    def test_get_depth_and_alt_freq_skipsFreqForIndels(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	99	chr10	5	0	3M	=	105	0	AAA	>>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)
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
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)
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
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)
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
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)
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
'''.format(self.simple_forward_read + _BamFlag.DUP)
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A")[0])
            
    def test_get_depth_and_alt_freq_failedQualityCounted(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	{}	chr10	5	0	5M	=	105	0	AAAAA	>>>>>
'''.format(self.simple_forward_read + _BamFlag.QCFAIL)
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A")[0])
            
    def test_get_depth_and_alt_freq_secondaryReadsCounted(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
readNameA	{}	chr10	5	0	5M	=	105	0	AAAAA	>>>>>
'''.format(self.simple_forward_read + _BamFlag.SECONDARY)
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, "test.sam", sam_contents)
            reader = zither._BamReader(input_bam)
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A")[0])

class ZitherTestCase(ZitherBaseTestCase):
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
            zither._create_vcf(input_vcf, sample_reader_dict)
        
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
             
            sample_reader_dict = {"sample_A": zither._BamReader(bam_A), "sample_B": zither._BamReader(bam_B)}
            zither._create_vcf(input_vcf, sample_reader_dict)
        
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
            zither._create_vcf(input_vcf, sample_reader_dict)
        
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
            zither._create_vcf(input_vcf, sample_reader_dict)

        
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
            bam_A = _create_bam(tmp_path, "sample_A.sam", sam_contents)

            sample_reader_dict = {"sample_A": zither._BamReader(bam_A)}
            zither._create_vcf(input_vcf, sample_reader_dict)
        
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
            actual_sample_reader_mapping = zither._build_reader_dict(sample_bam_mapping)
            self.assertEquals(["sA", "sB"],  sorted(actual_sample_reader_mapping.keys()))
            self.assertEquals(zither._BamReader(bam_A),  actual_sample_reader_mapping["sA"])
            self.assertEquals(zither._BamReader(bam_B),  actual_sample_reader_mapping["sB"])
    
    def test_get_sample_bam_strategy(self):
        input_vcf_contents = \
'''##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_A
chr1	10	.	A	C	.	.	.	GT	0/1
'''
        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            vcf_filename = _create_file(tmp_path, "foo.vcf", input_vcf_contents)
            
            args = Namespace(input_vcf=vcf_filename, mapping_file=None)
            actual_strategy = zither._get_sample_bam_strategy(args)
            self.assertIsInstance(actual_strategy, zither._MatchingNameStrategy)
            
            args = Namespace(input_vcf=vcf_filename, mapping_file="foo.txt")
            actual_strategy = zither._get_sample_bam_strategy(args)
            self.assertIsInstance(actual_strategy, zither._MappingFileStrategy)
    
    def test_main_matching_names(self):
        args = ["zither", "/ccmb/BioinfCore/SoftwareDev/projects/zither/examples/matching_names/input.vcf"]
        expected_vcf_filename ="/ccmb/BioinfCore/SoftwareDev/projects/zither/examples/matching_names/expected_output.vcf"
        with open(expected_vcf_filename, 'r') as expected_vcf:
            expected_vcf_contents = expected_vcf.read()
            zither.main(args)
            actual_output_lines = self.stdout.getvalue()    
            self._compare_lines(expected_vcf_contents, actual_output_lines)
    
class MatchingNameStrategyTestCase(ZitherBaseTestCase):
    def test_build_dict_from_matching_bams(self):
        actual_mapping = zither._MatchingNameStrategy(["sA", "sB"], "/foo/bar/input.vcf").build_sample_bam_mapping()        
        self.assertEquals(["sA", "sB"],  sorted(actual_mapping.keys()))
        self.assertEquals("/foo/bar/sA.bam",  actual_mapping["sA"])
        self.assertEquals("/foo/bar/sB.bam",  actual_mapping["sB"])

    def test_build_dict_from_matching_bams_emptySampleList(self):
        actual_mapping = zither._MatchingNameStrategy([], "/foo/bar/input.vcf").build_sample_bam_mapping()        
        self.assertEquals([],  actual_mapping.keys())

class MappingFileStrategyTestCase(ZitherBaseTestCase):      
    def test_build_dict_from_mapping_file(self):
        mapping_file_contents = \
'''sA	/foo/bar/sA.bam
sB	/foo/bar/sB.bam
'''
        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            mapping_file = _create_file(tmp_path, "mapping_file.txt", mapping_file_contents)

            actual_mapping = zither._MappingFileStrategy(mapping_file).build_sample_bam_mapping()
            self.assertEquals(["sA", "sB"],  sorted(actual_mapping.keys()))
            self.assertEquals("/foo/bar/sA.bam",  actual_mapping["sA"])
            self.assertEquals("/foo/bar/sB.bam",  actual_mapping["sB"])

    def test_build_dict_from_mapping_file_invalidMappingFile(self):
        mapping_file_contents = \
'''sA
sB	/foo/bar/sB.bam
'''
        with TempDirectory() as tmp_dir:
            tmp_path = tmp_dir.path
            mapping_file = _create_file(tmp_path, "mapping_file.txt", mapping_file_contents)
            
            strategy = zither._MappingFileStrategy(mapping_file)
            
            self.assertRaises(ValueError,
                strategy.build_sample_bam_mapping)
