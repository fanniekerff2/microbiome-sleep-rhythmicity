#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10
module load eth_proxy

qiime tools import \
      --type NCBIAccessionIDs \
      --input-path "$1/raw/project_id.tsv" \
      --output-path "$1/project_id.qza"

qiime fondue get-sequences \
      --i-accession-ids "$1/project_id.qza" \
      --p-email my@mail.ch \
      --o-single-reads "$2/single-reads.qza" \
      --o-paired-reads "$2/demux-paired-end.qza" \
      --o-failed-runs "$2/failed-ids.qza" \
      --p-n-jobs 5 \
      --p-retries 20 \
      --verbose