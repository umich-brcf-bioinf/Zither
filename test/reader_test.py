from __future__ import print_function, absolute_import, division

import bvc.reader as bvc_reader
import os
from testfixtures import TempDirectory
import unittest
from subprocess import call

def _create_bam(dir, sam_contents):
    test_input_sam = os.path.join(dir, "test.sam")
    sam_file = open(test_input_sam, "w")
    sam_file.write(sam_contents)
    sam_file.close()
    test_input_bam = os.path.join(dir, "test.bam")
    call(["samtools view -S -b "+ test_input_sam + " -o " + test_input_bam], shell=True)
    call(["samtools index " + test_input_bam], shell=True)
    return test_input_bam

class ReaderTestCase(unittest.TestCase):
    def test_get_depth_and_alt_freq(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
HWI-D00709:21:C6FJGANXX:1:2214:16769:34260	99	chr10	5	0	5M	=	105	0	AAAAA	>>>>>
HWI-D00709:21:C6FJGANXX:1:2113:10102:49540	99	chr10	6	0	5M	=	106	0	CCCCC	>>>>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = bvc_reader.Reader(input_bam)
            self.assertEquals((0, 0), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=4, ref="C", alt="A"))
            self.assertEquals((1, 1), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A"))
            self.assertEquals((2, .5), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=6, ref="C", alt="A"))

    def test_get_depth_and_alt_freq_ignoresDeletions(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
HWI-D00709:21:C6FJGANXX:1:2214:16769:34260	99	chr10	5	0	3M	=	105	0	AAA	>>>
HWI-D00709:21:C6FJGANXX:1:2113:10102:49540	99	chr10	5	0	1M1D1M	=	105	0	CG	>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = bvc_reader.Reader(input_bam)
            self.assertEquals(2, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A")[0])
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=6, ref="C", alt="A")[0])
            self.assertEquals(2, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=7, ref="C", alt="A")[0])

    def test_get_depth_and_alt_freq_ignoresSkippedRef(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
HWI-D00709:21:C6FJGANXX:1:2214:16769:34260	99	chr10	5	0	3M	=	105	0	AAA	>>>
HWI-D00709:21:C6FJGANXX:1:2113:10102:49540	99	chr10	5	0	1M1N1M	=	105	0	CG	>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = bvc_reader.Reader(input_bam)
            self.assertEquals(2, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A")[0])
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=6, ref="C", alt="A")[0])
            self.assertEquals(2, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=7, ref="C", alt="A")[0])

    def test_get_depth_and_alt_freq_ignoresCase(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
HWI-D00709:21:C6FJGANXX:1:2214:16769:34260	99	chr10	5	0	3M	=	105	0	AAA	>>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = bvc_reader.Reader(input_bam)
            self.assertEquals((1,1), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="a"))
            self.assertEquals((1,1), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A"))

    def test_get_depth_and_alt_freq_skipsFreqForIndels(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
HWI-D00709:21:C6FJGANXX:1:2214:16769:34260	99	chr10	5	0	3M	=	105	0	AAA	>>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = bvc_reader.Reader(input_bam)
            self.assertEquals((1,0), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="A", alt="ATA"))
            self.assertEquals((1,0), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="ATA", alt="A"))
 
    def test_get_depth_and_alt_freq_ignoresInsertions(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
HWI-D00709:21:C6FJGANXX:1:2214:16769:34260	99	chr10	5	0	8M	=	105	0	ACGTACGT	>>>>>>>>
HWI-D00709:21:C6FJGANXX:1:2214:16769:34260	99	chr10	5	0	4M2I4M	=	105	0	ACGTTTACGT	>>>>>>>>>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = bvc_reader.Reader(input_bam)
            self.assertEquals((2,1), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A"))
            self.assertEquals((2,1), reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=9, ref="C", alt="A"))

    def test_get_depth_and_alt_freq_overlappingReadPairsCountedSeparately(self):
        sam_contents = \
'''@HD	VN:1.4	GO:none	SO:coordinate
@SQ	SN:chr10	LN:135534747
HWI-D00709:21:C6FJGANXX:1:2214:16769:34260	99	chr10	5	0	4M	=	105	0	AAAA	>>>>
HWI-D00709:21:C6FJGANXX:1:2214:16769:34260	147	chr10	7	0	4M	=	105	0	TTTT	>>>>
'''
        with TempDirectory() as tmp_dir:    
            input_bam = _create_bam(tmp_dir.path, sam_contents)
            reader = bvc_reader.Reader(input_bam)
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=5, ref="C", alt="A")[0])
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=6, ref="C", alt="A")[0])
            self.assertEquals(2, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=7, ref="C", alt="A")[0])
            self.assertEquals(2, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=8, ref="C", alt="A")[0])
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=9, ref="C", alt="A")[0])
            self.assertEquals(1, reader.get_depth_and_alt_freq(chrom="chr10", position_one_based=10, ref="C", alt="A")[0])

            