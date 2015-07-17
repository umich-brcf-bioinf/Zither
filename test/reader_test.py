from __future__ import print_function, absolute_import, division

import bvc.reader as bvc_reader
import unittest

class ReaderTestCase(unittest.TestCase):
    def test_get_depth(self):
        sam_file_name = "test.sam"
        reader = bvc_reader.Reader(sam_file_name)
        self.assertEquals(0, reader.get_depth(chrom="chr10", position=4))
        self.assertEquals(1, reader.get_depth(chrom="chr10", position=5))
        self.assertEquals(2, reader.get_depth(chrom="chr10", position=6))
