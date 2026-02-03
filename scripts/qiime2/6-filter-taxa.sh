#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime taxa filter-table \
  --i-table "$1/2-filter-$3/filtered-feature-table.qza" \
  --i-taxonomy "$1/4-taxonomy-$3/taxonomy.qza" \
  --p-mode contains \
  --p-include p__ \
  --p-exclude 'p__;,Chloroplast,Mitochondria' \
  --o-filtered-table "$1/2-filter-$3/taxa-filtered-feature-table.qza"

qiime feature-table summarize \
  --i-table "$1/2-filter-$3/taxa-filtered-feature-table.qza" \
  --m-sample-metadata-file "$2" \
  --o-visualization "$1/2-filter-$3/taxa-filtered-feature-table.qzv"

qiime feature-table filter-seqs \
  --i-data "$1/2-filter-$3/filtered-asv-sequences.qza" \
  --i-table "$1/2-filter-$3/taxa-filtered-feature-table.qza" \
  --o-filtered-data "$1/2-filter-$3/taxa-filtered-asv-sequences.qza"

qiime feature-table tabulate-seqs \
  --i-data "$1/2-filter-$3/taxa-filtered-asv-sequences.qza" \
  --o-visualization "$1/2-filter-$3/taxa-filtered-asv-sequences.qzv"