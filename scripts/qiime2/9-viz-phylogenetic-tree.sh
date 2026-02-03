#!/bin/bash

source ${HOME}/.bashrc
conda activate empress

empress tree-plot \
    --i-tree "$1/5-phylogeny-$3/rooted_tree.qza" \
    --m-feature-metadata-file "$1/4-taxonomy-$3/taxonomy.qza" \
    --o-visualization "$1/5-phylogeny-$3/rooted_tree-viz.qzv"

empress community-plot \
    --i-tree "$1/5-phylogeny-$3/rooted_tree.qza" \
    --i-feature-table "$1/2-filter-$3/taxa-filtered-feature-table.qza" \
    --m-sample-metadata-file "$2" \
    --m-feature-metadata-file "$1/4-taxonomy-$3/taxonomy.qza" \
    --o-visualization "$1/5-phylogeny-$3/community-tree-viz.qzv"