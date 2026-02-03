#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime feature-table relative-frequency \
  --i-table "$1/$3-$2/$4-feature-table.qza" \
  --o-relative-frequency-table "$1/8-rel-freq-$2/$5-feature-table.qza"