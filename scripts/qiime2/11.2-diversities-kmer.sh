#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime kmerizer core-metrics \
    --i-sequences "$1/2-filter-$3/taxa-filtered-asv-sequences.qza" \
    --i-table "$1/2-filter-$3/taxa-filtered-feature-table.qza" \
    --p-sampling-depth $5 \
    --m-metadata-file "$2" \
    --p-kmer-size $4 \
    --p-no-with-replacement \
    --o-rarefied-table "$1/7-kmer-results-$3/diversity-kmer/rarefied-table.qza" \
    --o-kmer-table "$1/7-kmer-results-$3/diversity-kmer/kmer-table.qza" \
    --o-observed-features-vector "$1/7-kmer-results-$3/diversity-kmer/alpha-diversities/observed-features.qza" \
    --o-shannon-vector "$1/7-kmer-results-$3/diversity-kmer/alpha-diversities/shannon.qza" \
    --o-jaccard-distance-matrix "$1/7-kmer-results-$3/diversity-kmer/beta-div-distance-matrices/jaccard.qza" \
    --o-bray-curtis-distance-matrix "$1/7-kmer-results-$3/diversity-kmer/beta-div-distance-matrices/braycurtis.qza" \
    --o-jaccard-pcoa-results "$1/7-kmer-results-$3/diversity-kmer/pcoas-beta-div-metrics/jaccard-pcoa-results.qza" \
    --o-bray-curtis-pcoa-results "$1/7-kmer-results-$3/diversity-kmer/pcoas-beta-div-metrics/bray-curtis-pcoa-results.qza" \
    --o-scatterplot "$1/7-kmer-results-$3/diversity-kmer/scatterplot.qzv"