#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime boots core-metrics \
    --i-table "$1/2-filter-$3/taxa-filtered-feature-table.qza" \
    --i-phylogeny "$1/5-phylogeny-$3/rooted_tree.qza" \
    --m-metadata-file "$2" \
    --p-sampling-depth $5 \
    --p-n $4 \
    --p-no-replacement \
    --p-alpha-average-method median \
    --p-beta-average-method medoid \
    --o-resampled-tables "$1/6-diversity-$3/resampled-tables/" \
    --o-alpha-diversities "$1/6-diversity-$3/alpha-diversities/" \
    --o-distance-matrices "$1/6-diversity-$3/beta-div-distance-matrices/" \
    --o-pcoas "$1/6-diversity-$3/pcoas-beta-div-metrics/" \
    --o-emperor-plots "$1/6-diversity-$3/emperor-plots-beta-div-metrics/"

qiime diversity alpha-group-significance \
  --i-alpha-diversity "$1/6-diversity-$3/alpha-diversities/observed_features.qza" \
  --m-metadata-file "$2" \
  --o-visualization "$1/6-diversity-$3/alpha-diversities/observed_features.qzv"

qiime diversity alpha-group-significance \
  --i-alpha-diversity "$1/6-diversity-$3/alpha-diversities/shannon.qza" \
  --m-metadata-file "$2" \
  --o-visualization "$1/6-diversity-$3/alpha-diversities/shannon.qzv"

qiime diversity alpha-group-significance \
  --i-alpha-diversity "$1/6-diversity-$3/alpha-diversities/pielou_e.qza" \
  --m-metadata-file "$2" \
  --o-visualization "$1/6-diversity-$3/alpha-diversities/pielou_e.qzv"

qiime diversity alpha-group-significance \
  --i-alpha-diversity "$1/6-diversity-$3/alpha-diversities/faith_pd.qza" \
  --m-metadata-file "$2" \
  --o-visualization "$1/6-diversity-$3/alpha-diversities/faith_pd.qzv"