Installing Zither
=================
Zither has been tested with Python 2.7 and Python 3.4 on \*nix and OSX.

Prerequisites
-------------
* pysam (0.8.3)  
* nosetests, testfixtures (3.0.2) are required for running
  automated tests
* Note that pip installs all required libraries; see [Installing] below.

Installing
----------
You can install from source from github:

``$ pip install git+https://github.com/umich-brcf-bioinf/Zither``

If you don't have root permissions, you can install locally:

``$ pip install git+https://github.com/umich-brcf-bioinf/Zither --user``

Following the pip install, you may need to adjust your path settings to include home/.local/bin. 


If you already have pysam installed, you can also clone from github and run directly from the source like so:

``$ git clone https://github.com/umich-brcf-bioinf/Zither``

``$ Zither/zither-runner.py --bam input.bam input.vcf > zither.vcf``
  

