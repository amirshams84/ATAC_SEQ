#!/bin/bash
## conda environment name

source /data/shamsaddinisha/conda/etc/profile.d/conda.sh

ENV_NAME_PY3=ATAC_Seq_py3

conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda


conda activate $ENV_NAME_PY3

conda install bc=1.06.0 cython=0.28.5 boost=1.68.0 blas=1.1.0 openblas=0.3.3 ncurses=6.1.0 pandas=0.23.4 perl=5.26.2.0 libgfortran=3.0.0 pigz=2.3.4 zlib=1.2.11 graphviz=2.38.0 python-levenshtein=0.12.0 jinja2=2.10 gsl=2.4.0 imagemagick=7.0 matplotlib=2.2.0 -y
conda install wget=1.19.5 nomkl=3.0.0 numpy=1.14 numpy-base=1.14 scipy=1.1.0 six=1.11.0 python-dateutil=2.7.3 -y
conda install pysam=0.15.1 pybedtools=0.7.10 pyfaidx=0.5.3 -y
conda install perl-threaded=5.22.0 homer=4.9.1 meme=5.0.2 bedtools=2.27.1 idr=2.0.4.2 -y

conda deactivate






echo "Installing ${ENV_NAME_PY3} has been successfully done."