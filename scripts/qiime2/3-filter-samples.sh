#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime feature-table filter-samples \
    --i-table "$1/1-denoise/dada2-feature-table.qza" \
    --p-min-frequency $2 \
    --m-metadata-file "$3" \
    --o-filtered-table "$1/2-filter-$2/temp-filt-feature-table.qza"

qiime feature-table summarize \
  --i-table "$1/2-filter-$2/temp-filt-feature-table.qza" \
  --m-sample-metadata-file "$3" \
  --o-visualization "$1/2-filter-$2/temp-filt-feature-table.qzv"