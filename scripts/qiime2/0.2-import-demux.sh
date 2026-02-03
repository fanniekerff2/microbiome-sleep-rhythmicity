#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime demux summarize \
    --i-data "$1/demux-paired-end.qza" \
    --o-visualization "$2/0-import-demux/demux-paired-end.qzv"