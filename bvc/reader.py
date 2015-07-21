from __future__ import print_function, absolute_import, division
import pysam

class Reader(object):
    def __init__(self, sam_file_name):
        self.samfile = pysam.AlignmentFile(sam_file_name, "rb" )

    def get_depth_and_alt_freq(self, chrom, position_one_based, ref, alt):
        alt = alt.upper()
        total_depth = 0
        total_AF = 0
        AF = 0
        position_zero_based = position_one_based - 1
        
        for pileupcolumn in self.samfile.pileup(chrom, position_zero_based, position_one_based):
            if (pileupcolumn.reference_pos == position_zero_based):
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        total_depth += 1

                        if len(ref)==1 and pileupread.alignment.query_sequence[pileupread.query_position].upper() == alt:
                            total_AF += 1

                break
        if total_depth:
            AF = total_AF/total_depth
        return (total_depth, AF)
