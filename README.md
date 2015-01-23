bntools
=======

Tools for BioNano data analysis

How to install
--------------

Download source from github:

    git clone http://github.com/yanlinlin82/bntools

or

    wget -N https://github.com/yanlinlin82/bntools/archive/master.zip
    unzip master.zip && mv bntools-master/ bntools/

Then build it:

    cd bntools/ && make

Quick start
-----------

1. Nick reference genome sequences (in FASTA format) to generate nick site map
   (in Tab-Separated-Values format).

       bntools nick hg38.fa -o hg38.tsv.gz

2. Map molecules onto the reference map.

       bntools map hg38.tsv.gz input.bnx.gz

3. Convert maps between formats (txt/tsv/bnx/cmap).

       bntools view input.bnx.gz -f cmap

Others
------

1. Gzip compression input/output files are supported. If you specify an output
   filename with ".gz" suffix, bntools will save it with gzip compression.

2. To build with debug info, try:

    make DEBUG=1
