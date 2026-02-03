#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime taxa barplot \
  --i-table "$1/2-filter-$3/taxa-filtered-feature-table.qza" \
  --i-taxonomy "$1/4-taxonomy-$3/taxonomy.qza" \
  --m-metadata-file "$2" \
  --o-visualization "$1/4-taxonomy-$3/taxa-barplot.qzv"