#! /usr/bin/bash

# To map the region of the O. biroi reference genome into which the transgene was inserted, we started by sequencing
# the whole genome of a transgenic individual using short read whole genome sequencing.
# We aligned that data to the sequence of the transgene insert as a reference genome.
# We loaded that aligned data in the integrative genomics viewer, and identified reads 
# that mapped to the beginning or the end of the transgene insert (and thus likely also aligned to the region of the
# O. biroi endogenous genome into which the transgene was inserted).

# After Identifying these "junction reads" we aligned the whole genome DNA sequencing data to the species' reference genome
# and queried those alignments by the names of each junction read using "samtools view"
# This revealed the chromosomal location to which each of these reads aligned in the species' reference genome, and we confirmed that they all mapped in the same location
for READ_NAME in A00815:199:HC3G5DRXY:2:2123:24885:3537; do
	samtools view /Path/To/DNAseq/Alignment/To/Species/Reference/SAMPLE_Trimmed_aligned_sorted_merged_dedup.bam | grep -m 2 $READ_NAME
done

# Next we created a consensus sequence of these reads with the R package 'msa' using CLUSTAL 2.1. (See ./MSA.R)
# The part of this consensus sequence that was identical to the transgene insert was removed,
# and the remaining sequence was assumed to be the sequence of the insertion site.

# We used BLAST on the NCBI website to identify the insertion site. 

