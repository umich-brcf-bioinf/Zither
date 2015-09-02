====== 
Zither
======

Command-line tool to pull raw depths and alt freqs from BAM file(s) based on an existing VCF, writing output as new VCF to stdout.

.. image:: https://travis-ci.org/umich-brcf-bioinf/Zither.svg?branch=develop
    :target: https://travis-ci.org/umich-brcf-bioinf/Zither
    :alt: Build Status

.. image:: https://coveralls.io/repos/umich-brcf-bioinf/zither/badge.svg?branch=develop&service=github
    :target: https://coveralls.io/github/umich-brcf-bioinf/zither?branch=develop
    :alt: Coverage Status



The official repository is at:

https://github.com/umich-brcf-bioinf/Zither

----------
Quickstart
----------

Read a single BAM file
======================

   $ zither --bam examples/explicit_bam/Sample_X.bam examples/explicit_bam/input.vcf > output.vcf 

Given a VCF and a BAM file, read positions in the input VCF and corresponding pileups 
from Sample_X.bam.


Read a set of matched VCF sample names and BAM files
====================================================

   $ zither examples/matching_names/input.vcf > output.vcf 

Given a VCF and a collection of BAM files whose file names 
match the VCF sample names, reads positions from the 
input VCF and corresponding BAM pileups.


Explicitly map VCF sample names to BAM files
====================================================

   $ zither --mapping_file=examples/mapping_files/mapping_file.txt examples/mapping_files/input.vcf > output.vcf 

Given a VCF, a collection of BAMs, and a file that maps sample names to BAM paths,
reads positions from the input VCF and corresponding pileups 
from BAM files names. 

The mapping file is a tab-separated text file where each line has a sample 
name and the path to the corresponding BAM file. Paths to BAM files can be 
absolute or relative; relative paths are resolved relative to the directory 
that contains the mapping file. 

====

Email bfx-zither@umich.edu for support and questions.

UM BRCF Bioinformatics Core 
