Installing Zither
==================
Zither has been tested with Python 2.7 on \*nix.

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
