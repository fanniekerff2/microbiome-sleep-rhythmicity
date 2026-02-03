#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime diversity alpha-rarefaction \
  --i-table "$1/1-denoise/dada2-feature-table.qza" \
  --m-metadata-file "$2" \
  --p-max-depth $3 \
  --p-steps $4 \
  --o-visualization "$1/1.2-rarefaction/alpha-rarefaction-plot-$4.qzv"