#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2025.7

# Filter features to those that are present at least 25 times in at least 4 samples
qiime feature-table filter-features \
  --i-table "$1/$4-$3/$5-feature-table.qza" \
  --p-min-frequency 25 \
  --p-min-samples 4 \
  --o-filtered-table "$1/10-diff-abund-$3/DA-$5-feature-table.qza"

# Run ANCOM-BC2 to identify differentially abundant taxa
qiime composition ancombc2 \
  --i-table "$1/10-diff-abund-$3/DA-$5-feature-table.qza" \
  --m-metadata-file "$2" \
  --p-fixed-effects-formula 'melatonin + age_days + sex' \
  --p-random-effects-formula '(1 | infant_id)' \
  --p-p-adjust-method $6 \
  --p-prevalence-cutoff $7 \
  --o-ancombc2-output "$1/10-diff-abund-$3/$5-ancombc_results-$6-$7.qza"

# Generate a table of these same values for all taxa
qiime composition tabulate \
  --i-data "$1/10-diff-abund-$3/$5-ancombc_results-$6-$7.qza" \
  --o-visualization "$1/10-diff-abund-$3/$5-ancombc_results-$6-$7.qzv"
    
# Generate a visualization of the ANCOM-BC2 results
qiime composition ancombc2-visualizer \
  --i-data "$1/10-diff-abund-$3/$5-ancombc_results-$6-$7.qza" \
  --o-visualization "$1/10-diff-abund-$3/$5-ancombc2_visualization-$6-$7.qzv"