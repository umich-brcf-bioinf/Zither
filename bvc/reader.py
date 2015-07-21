from __future__ import print_function, absolute_import, division
import pysam

class Reader(object):
    def __init__(self, sam_file_name):
        self.samfile = pysam.AlignmentFile(sam_file_name, "rb" )

    def get_depth_and_alt_freq(self, chrom, position_one_based, ref, alt):
        BASE_INDEX = {'A':0, 'C':1, 'G':2, 'T':3}        
        alt = alt.upper()
        AF = 0
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
                AF = variant_count/total_depth
        except KeyError:
            pass #AF remains 0

        return (total_depth, AF)
