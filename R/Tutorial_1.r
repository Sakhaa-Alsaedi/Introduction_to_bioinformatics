#######Task 1#######

# Open the file for reading
file <- file('gene.fna', 'r')

# Initialize an empty string to store the gene sequences
gene_sequences <- ""

# Read each line in the file
while (length(line <- readLines(file, n = 1)) > 0) {
  # Check if the line starts with '>'
  if (startsWith(line, '>')) {
    # This is a header line, you can skip it or process it as needed
    next
  } else {
    # Append the gene sequence to the gene_sequences string
    gene_sequences <- paste(gene_sequences, line, sep = "")
  }
}

# Close the file
close(file)

# Print the entire gene sequence
cat(gene_sequences, "\n")

#######Task 2#######

find_ori_location <- function(genome) {
  ori <- "ATACAATA"
  
  for (i in 1:(nchar(genome) - nchar(ori) + 1)) {
    if (substring(genome, i, i + nchar(ori) - 1) == ori) {
      return(i)
    }
  }
  
  return(-1)  # Return -1 if ori is not found in the genome
}

ori_location <- find_ori_location(gene_sequences)

if (ori_location != -1) {
  cat("The origin of replication (ori) is found at position", ori_location, "\n")
} else {
  cat("The origin of replication (ori) was not found in the genome.\n")
}

#######Task 3#######

# Define the DNA sequence
dna_sequence <- gene_sequences

# Initialize an empty string to store the RNA sequence
rna_sequence <- ""

# Transcribe DNA to RNA by replacing 'T' with 'U'
for (base in strsplit(dna_sequence, "")[[1]]) {
  if (base == "T") {
    rna_sequence <- paste(rna_sequence, "U", sep = "")
  } else {
    rna_sequence <- paste(rna_sequence, base, sep = "")
  }
}

# Print the RNA sequence
cat(rna_sequence, "\n")

#######Task 4#######

# Load the library for regular expressions
# First, you need to install it using this command: install.packages('stringr')
library(stringr)

# Define the DNA sequence
dna_sequence <- gene_sequences

# Define the motif you want to search for
motif <- "TACGT"

# Use regular expressions to search for the motif
matches <- str_locate_all(dna_sequence, motif)

# Initialize a list to store the positions of motif matches
match_positions <- list()

# Iterate over the matches and record their positions
for (match in matches[[1]]) {
  start <- match[1]
  end <- match[2]
  match_positions <- append(match_positions, list(c(start, end)))
}

# Print the positions of motif matches
if (length(match_positions) > 0) {
  cat(paste("The motif '", motif, "' was found at the following positions:\n", sep = ""))
  for (match in match_positions) {
    cat("Start:", match[1], ", End:", match[2], "\n")
  }
} else {
  cat(paste("The motif '", motif, "' was not found in the DNA sequence.\n", sep = ""))
}

#######Task 5#######

# Define the DNA sequence
dna_sequence <- gene_sequences

# Calculate the GC content
gc_count <- sum(str_count(dna_sequence, "G") + str_count(dna_sequence, "C"))
sequence_length <- nchar(dna_sequence)
gc_content <- (gc_count / sequence_length) * 100

# Print the GC content
cat("GC Content: ", sprintf("%.2f%%", gc_content), "\n")

#######Task 6#######

# Define a named list for RNA to amino acid translation
rna_codon_table <- list(
  "AUG" = "M",  # Start codon
  "UUU" = "F",  "UUC" = "F",
  "UUA" = "L",  "UUG" = "L",
  "UCU" = "S",  "UCC" = "S",  "UCA" = "S",  "UCG" = "S",
  "UAU" = "Y",  "UAC" = "Y",
  "UGU" = "C",  "UGC" = "C",
  "UGG" = "W",
  "CUU" = "L",  "CUC" = "L",  "CUA" = "L",  "CUG" = "L",
  "CCU" = "P",  "CCC" = "P",  "CCA" = "P",  "CCG" = "P",
  "CAU" = "H",  "CAC" = "H",
  "CAA" = "Q",  "CAG" = "Q",
  "CGU" = "R",  "CGC" = "R",  "CGA" = "R",  "CGG" = "R",
  "AUU" = "I",  "AUC" = "I",  "AUA" = "I",
  "ACU" = "T",  "ACC" = "T",  "ACA" = "T",  "ACG" = "T",
  "AAU" = "N",  "AAC" = "N",
  "AAA" = "K",  "AAG" = "K",
  "AGU" = "S",  "AGC" = "S",
  "AGA" = "R",  "AGG" = "R",
  "GUU" = "V",  "GUC" = "V",  "GUA" = "V",  "GUG" = "V",
  "GCU" = "A",  "GCC" = "A",  "GCA" = "A",  "GCG" = "A",
  "GAU" = "D",  "GAC" = "D",
  "GAA" = "E",  "GAG" = "E",
  "GGU" = "G",  "GGC" = "G",  "GGA" = "G",  "GGG" = "G",
  "UAA" = "*",  "UAG" = "*",  "UGA" = "*"  # Stop codons
)

# Initialize an empty string to store the protein sequence
protein_sequence <- ""

# Translate the RNA sequence into a protein sequence
for (i in seq(1, nchar(rna_sequence), by = 3)) {
  codon <- substr(rna_sequence, i, i + 2)
  amino_acid <- rna_codon_table[[codon]]
  if (is.null(amino_acid)) {
    amino_acid <- "X"  # 'X' for unknown codons
  }
  protein_sequence <- paste0(protein_sequence, amino_acid)
}

# Print the protein sequence
cat("Protein Sequence:", protein_sequence, "\n")
