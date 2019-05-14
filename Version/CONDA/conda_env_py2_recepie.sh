#!/bin/bash
## conda environment name

source /data/shamsaddinisha/conda/etc/profile.d/conda.sh

ENV_NAME_PY2=ATAC_Seq_py2

conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

conda activate $ENV_NAME_PY2

conda install bc=1.06.0 cython=0.28.5 boost=1.68.0 blas=1.1.0 openblas=0.3.3 ncurses=6.1.0 pandas=0.23.4 perl=5.26.2.0 libgfortran=3.0.0 pigz=2.3.4 zlib=1.2.11 graphviz=2.38.0 python-levenshtein=0.12.0 jinja2=2.10 gsl=2.4.0 imagemagick=7.0 matplotlib=2.2.0 -y
conda install r=3.2 r-bitops=1.0 r-catools=1.17 -y
conda install r-snow=0.3 r-snowfall=1.84 bioconductor-rsamtools bioconductor-deseq2 r-spp=1.13 ghostscript=9.18.0 trim-galore=0.5.0 fastqc=0.11.7 afterqc=0.9.7 qualimap=2.2.2a multiqc=1.6.0 pysam=0.15.1 pybedtools=0.7.10 pyfaidx=0.5.3 -y
conda install wget=1.19.5 nomkl=3.0.0 numpy=1.14 numpy-base=1.14 scipy=1.1.0 six=1.11.0 python-dateutil=2.7.3 -y
conda install ucsc-bedgraphtobigwig=366 ucsc-fetchchromsizes=366 ucsc-wigtobigwig=366 ucsc-bigwiginfo=366 ucsc-bedclip=366 ucsc-bedtobigbed=366 ucsc-bigbedtobed=366 ucsc-bigwigtowig=366 ucsc-twobittofa=366 ucsc-fatotwobit=366 ucsc-bigwiginfo=366 ucsc-bigwigtobedgraph=366  ucsc-liftover=366 ucsc-bigWigAverageOverBed=366 -y
conda install ucsc-bigwigmerge=366 ucsc-bigbedinfo=366 ucsc-bigbedsummary=366 -y
conda install bedtools=2.27.1 macs2=2.1.1.20160309 cutadapt=1.18.0 bedops=2.4.35 bowtie=1.2.2 bowtie2=2.3.4.2 bwa=0.7.17 samtools=1.9.0 sambamba=0.6.6 picard=2.18.7 metaseq=0.5 preseq=2.0.3 -y

conda deactivate





echo "Installing ${ENV_NAME_PY2} has been successfully done."