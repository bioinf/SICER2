#!/usr/bin/env/ python

from distutils.core import setup, Extension


def main():
    setup(name="broadPeaks",
          version="1.0",
          description="Tool for ChiP-seq analysis to find broad peaks",
          author='Bioinformatics institute',
          author_email='dmkr0001@gmail.com; prikaziuk@gmail.com; predeus@gmail.com',
          url='https://github.com/BroadPeaksBioinf/BroadPeaks',
          package_dir={'BroadPeaks1': './BroadPeaks1'},
          packages=['BroadPeaks1'],
          scripts=['BroadPeaks1/broadPeaks'])

if __name__ == '__main__':
    main()
