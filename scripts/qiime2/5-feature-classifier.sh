#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime feature-classifier classify-sklearn \
  --i-classifier "$3" \
  --i-reads "$1/2-filter-$2/filtered-asv-sequences.qza" \
  --o-classification "$1/4-taxonomy-$2/taxonomy.qza"

qiime metadata tabulate \
  --m-input-file "$1/4-taxonomy-$2/taxonomy.qza" \
  --o-visualization "$1/4-taxonomy-$2/taxonomy.qzv"