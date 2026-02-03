#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime feature-table filter-features \
  --i-table "$1/2-filter-$2/temp-filt-feature-table.qza" \
  --p-min-samples $4 \
  --p-min-frequency $5 \
  --o-filtered-table "$1/2-filter-$2/filtered-feature-table.qza"

qiime feature-table summarize \
  --i-table "$1/2-filter-$2/filtered-feature-table.qza" \
  --m-sample-metadata-file "$3" \
  --o-visualization "$1/2-filter-$2/filtered-feature-table.qzv"

qiime feature-table filter-seqs \
  --i-data "$1/1-denoise/dada2-asv-sequences.qza" \
  --i-table "$1/2-filter-$2/filtered-feature-table.qza" \
  --o-filtered-data "$1/2-filter-$2/filtered-asv-sequences.qza"

qiime feature-table tabulate-seqs \
  --i-data "$1/2-filter-$2/filtered-asv-sequences.qza" \
  --o-visualization "$1/2-filter-$2/filtered-asv-sequences.qzv"