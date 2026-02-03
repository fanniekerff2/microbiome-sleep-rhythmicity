#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$1/seq_data/16S-GB" \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path "$1/data-sensitive/seq_data/demux-paired-end.qza"