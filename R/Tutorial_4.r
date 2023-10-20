# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "clusterProfiler", "org.Hs.eg.db"))
install.packages("seqinr")

library(Biostrings)
library(seqinr)
library(clusterProfiler)
library(org.Hs.eg.db)

# Define and analyze a protein sequence
protein_sequence <- "MGPAGDVRQERIKWVRRSYW"
protein_obj <- AAString(protein_sequence)

# Simple analysis like finding motifs
motif <- "GDV"
if (matchPattern(motif, protein_obj)) {
  print(paste("Motif", motif, "found in sequence."))
} else {
  print(paste("Motif", motif, "not found in sequence."))
}

# Calculate molecular weight of the protein
mol_weight <- seqinr::mw(protein_sequence, seqtype = "protein", letters = AA_ALPHABET)
print(paste("The molecular weight of the protein is", mol_weight))

# Parse sequences from a FASTA file (assuming 'sequences.fasta' is available in the working directory)
fasta_file <- "sequences.fasta"
sequences <- readAAStringSet(fasta_file)

# Displaying the sequences (optional step for verification)
print(sequences)

# Enrichment Analysis - assuming you have a gene list from some source; here, we use a placeholder list
gene_list <- c("7157", "1956", "672", "324", "7124", "2475", "207")  # Replace with your actual list

# Performing GO enrichment analysis
ego <- enrichGO(gene         = gene_list,
                universe     = keys(org.Hs.eg.db, keytype = "ENTREZID"),
                keyType      = 'ENTREZID',
                OrgDb        = org.Hs.eg.db,
                ont          = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable     = TRUE)

# Viewing and visualizing the results
print(ego)
barplot(ego, showCategory=10)
dotplot(ego, showCategory=10)
