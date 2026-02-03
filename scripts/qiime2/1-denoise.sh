#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$5/demux-paired-end.qza" \
  --p-trunc-len-f $2 \
  --p-trunc-len-r $3 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --o-representative-sequences "$1/1-denoise/dada2-asv-sequences.qza" \
  --o-table "$1/1-denoise/dada2-feature-table.qza" \
  --o-denoising-stats "$1/1-denoise/dada2-stats.qza"

qiime metadata tabulate \
  --m-input-file "$1/1-denoise/dada2-stats.qza" \
  --o-visualization "$1/1-denoise/dada2-stats.qzv"

qiime feature-table summarize \
  --i-table "$1/1-denoise/dada2-feature-table.qza" \
  --m-sample-metadata-file "$4" \
  --o-visualization "$1/1-denoise/dada2-feature-table.qzv"

qiime feature-table tabulate-seqs \
  --i-data "$1/1-denoise/dada2-asv-sequences.qza" \
  --o-visualization "$1/1-denoise/dada2-asv-sequences.qzv"