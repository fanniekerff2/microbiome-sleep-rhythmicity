#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "$1/2-filter-$2/taxa-filtered-asv-sequences.qza" \
  --output-dir "$1/5-phylogeny-$2"