from setuptools import setup
from setuptools import find_packages
import os
import zither

def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), 'r') as filename:
        return filename.read()

setup(name='zither',
      version=zither.__version__,
      description=('Command-line tool to pull raw depths and alt freqs from '
                   'BAM file(s) based on an existing VCF, writing output as '
                   'new VCF to stdout.'),
      long_description=(read('README.rst') + '\n\n' +
                        read('CHANGELOG.rst') + '\n\n' +
                        read('AUTHORS.rst')),
      url='https://github.com/umich-brcf-bioinf/Zither',
      author='University of Michigan Bioinformatics Core',
      author_email='bfx-zither@umich.edu',
      license='Apache',
      packages=find_packages(exclude=['test*']),
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: Apache Software License',
                   'Operating System :: Unix',
                   'Programming Language :: Python :: 2.7',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'],
      keywords='VCF bioinformatic exome-seq DNA-seq variant-call-format BAM',
      install_requires=['pysam'],
      entry_points={'console_scripts': ['zither=zither.zither:main']},
      test_suite='nose.collector',
      tests_require=['nose', 'testfixtures', 'pysam'],
      zip_safe=False)
