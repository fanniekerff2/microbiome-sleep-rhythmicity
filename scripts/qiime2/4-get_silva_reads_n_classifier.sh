#!/bin/bash

source ${HOME}/.bashrc
conda activate qiime2-amplicon-2024.10

qiime rescript reverse-transcribe \
    --i-rna-sequences "$1/silva-138.1-ssu-nr99-rna-seqs.qza" \
    --o-dna-sequences "$1/silva-138.1-ssu-nr99-seqs.qza"

qiime rescript cull-seqs \
    --i-sequences "$1/silva-138.1-ssu-nr99-seqs.qza" \
    --o-clean-sequences "$1/silva-138.1-ssu-nr99-seqs-cleaned.qza"

qiime rescript filter-seqs-length-by-taxon \
    --i-sequences "$1/silva-138.1-ssu-nr99-seqs-cleaned.qza" \
    --i-taxonomy "$1/silva-138.1-ssu-nr99-tax.qza" \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs "$1/silva-138.1-ssu-nr99-seqs-filt.qza" \
    --o-discarded-seqs "$1/silva-138.1-ssu-nr99-seqs-discard.qza"

qiime rescript dereplicate \
    --i-sequences "$1/silva-138.1-ssu-nr99-seqs-filt.qza" \
    --i-taxa "$1/silva-138.1-ssu-nr99-tax.qza" \
    --p-mode 'uniq' \
    --o-dereplicated-sequences "$1/silva-138.1-ssu-nr99-seqs-derep-uniq.qza" \
    --o-dereplicated-taxa "$1/silva-138.1-ssu-nr99-tax-derep-uniq.qza"

qiime feature-classifier extract-reads \
    --i-sequences "$1/silva-138.1-ssu-nr99-seqs-derep-uniq.qza" \
    --p-f-primer GTGYCAGCMGCCGCGGTAA \
    --p-r-primer GGACTACNVGGGTWTCTAAT \
    --p-read-orientation 'forward' \
    --o-reads "$1/silva-138.1-ssu-nr99-seqs-515f-806r.qza"

qiime rescript dereplicate \
    --i-sequences "$1/silva-138.1-ssu-nr99-seqs-515f-806r.qza" \
    --i-taxa "$1/silva-138.1-ssu-nr99-tax-derep-uniq.qza" \
    --p-mode 'uniq' \
    --o-dereplicated-sequences "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq.qza" \
    --o-dereplicated-taxa  "$1/silva-138.1-ssu-nr99-tax-515f-806r-derep-uniq.qza"
# ! are the `o-dereplicated-sequences` the SILVA reference v4 reads that can be used for closed reference clustering?

qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$1/silva-138.1-ssu-nr99-seqs-515f-806r-uniq.qza" \
    --i-reference-taxonomy "$1/silva-138.1-ssu-nr99-tax-515f-806r-derep-uniq.qza" \
    --o-classifier "$1/silva-138.1-ssu-nr99-515f-806r-classifier.qza"

# remove unneeded files
rm -r "$1/silva-138.1-ssu-nr99-seqs.qza"
rm -r "$1/silva-138.1-ssu-nr99-seqs-cleaned.qza"
rm -r "$1/silva-138.1-ssu-nr99-seqs-discard.qza"
rm -r "$1/silva-138.1-ssu-nr99-seqs-filt.qza"
rm -r "$1/silva-138.1-ssu-nr99-seqs-derep-uniq.qza"
rm -r "$1/silva-138.1-ssu-nr99-tax-derep-uniq.qza"
rm -r "$1/silva-138.1-ssu-nr99-seqs-515f-806r.qza"