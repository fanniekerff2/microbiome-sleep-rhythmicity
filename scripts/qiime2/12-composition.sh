#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime taxa collapse \
  --i-table "$1/2-filter-$4/taxa-filtered-feature-table.qza" \
  --i-taxonomy "$1/4-taxonomy-$4/taxonomy.qza" \
  --p-level $2 \
  --o-collapsed-table "$1/9-comp-$4/$3-feature-table.qza"

qiime feature-table summarize \
  --i-table "$1/9-comp-$4/$3-feature-table.qza" \
  --m-sample-metadata-file "$5" \
  --o-visualization "$1/9-comp-$4/$3-feature-table.qzv"