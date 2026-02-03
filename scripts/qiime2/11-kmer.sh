#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime kmerizer seqs-to-kmers \
    --i-sequences "$1/2-filter-$2/taxa-filtered-asv-sequences.qza" \
    --i-table "$1/2-filter-$2/taxa-filtered-feature-table.qza" \
    --p-kmer-size $3 \
    --p-tfidf False \
    --p-max-features $4 \
    --o-kmer-table "$1/7-kmer-results-$2/kmer-table-$4.qza"
    
qiime feature-table relative-frequency \
  --i-table "$1/7-kmer-results-$2/kmer-table-$4.qza" \
  --o-relative-frequency-table "$1/7-kmer-results-$2/rel-freq-kmer-table-$4.qza"


qiime kmerizer seqs-to-kmers \
    --i-sequences "$1/2-filter-$2/taxa-filtered-asv-sequences.qza" \
    --i-table "$1/2-filter-$2/taxa-filtered-feature-table.qza" \
    --p-kmer-size $3 \
    --p-tfidf True \
    --o-kmer-table "$1/7-kmer-results-$2/tfidf-kmer-table-$4.qza" \
    --p-max-features $4

qiime feature-table relative-frequency \
  --i-table "$1/7-kmer-results-$2/tfidf-kmer-table-$4.qza" \
  --o-relative-frequency-table "$1/7-kmer-results-$2/rel-freq-tfidf-kmer-table-$4.qza"