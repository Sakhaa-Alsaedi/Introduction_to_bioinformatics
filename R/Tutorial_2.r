#######Task 1#######

global_alignment <- function(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-2) {
  m <- nchar(seq1)
  n <- nchar(seq2)
  
  score_matrix <- matrix(0, nrow=m+1, ncol=n+1)
  
  # Initialization
  for(i in 1:(m+1)) {
    score_matrix[i, 1] <- gap_penalty * (i-1)
  }
  
  for(j in 1:(n+1)) {
    score_matrix[1, j] <- gap_penalty * (j-1)
  }
  
  # Fill the score matrix
  for(i in 2:(m+1)) {
    for(j in 2:(n+1)) {
      match <- score_matrix[i-1, j-1] + ifelse(substr(seq1, i-1, i-1) == substr(seq2, j-1, j-1), match_score, mismatch_score)
      delete <- score_matrix[i-1, j] + gap_penalty
      insert <- score_matrix[i, j-1] + gap_penalty
      
      score_matrix[i, j] <- max(match, delete, insert)
    }
  }
  
  # Traceback and compute the alignment 
  align1 <- ""
  align2 <- ""
  i <- m+1
  j <- n+1
  while (i > 1 || j > 1) {
    current_score <- score_matrix[i, j]
    if (i > 1 && score_matrix[i, j] == score_matrix[i-1, j] + gap_penalty) {
      align1 <- paste0(substr(seq1, i-1, i-1), align1)
      align2 <- paste0("-", align2)
      i <- i - 1
    } else if (j > 1 && score_matrix[i, j] == score_matrix[i, j-1] + gap_penalty) {
      align1 <- paste0("-", align1)
      align2 <- paste0(substr(seq2, j-1, j-1), align2)
      j <- j - 1
    } else {
      align1 <- paste0(substr(seq1, i-1, i-1), align1)
      align2 <- paste0(substr(seq2, j-1, j-1), align2)
      i <- i - 1
      j <- j - 1
    }
  }
  
  return(list(align1, align2))
}

# Test the function
seq1 <- "ATCG"
seq2 <- "TCA"
alignment <- global_alignment(seq1, seq2)
print(alignment[[1]])
print(alignment[[2]])

#######Task 2#######

calculate_allele_frequencies <- function(sequences) {
  # Calculate allele frequencies at each position for a list of aligned sequences.
  
  # Initialize a list for each position
  num_sequences <- length(sequences)
  sequence_length <- nchar(sequences[[1]])
  allele_counts <- vector("list", sequence_length)
  
  # Initialize the lists within allele_counts
  for(i in 1:sequence_length) {
    allele_counts[[i]] <- list()
  }
  
  # Count alleles at each position
  for(seq in sequences) {
    for(i in 1:sequence_length) {
      nucleotide <- substr(seq, i, i)
      if(!is.null(allele_counts[[i]][[nucleotide]])) {
        allele_counts[[i]][[nucleotide]] <- allele_counts[[i]][[nucleotide]] + 1
      } else {
        allele_counts[[i]][[nucleotide]] <- 1
      }
    }
  }
  
  # Convert counts to frequencies
  allele_frequencies <- list()
  for(count in allele_counts) {
    frequencies <- list()
    for(allele in names(count)) {
      frequencies[[allele]] <- count[[allele]] / num_sequences
    }
    allele_frequencies <- append(allele_frequencies, list(frequencies))
  }
  
  return(allele_frequencies)
}

# Test the function
sequences <- c(
  "ATCG",
  "ATCG",
  "ATCA",
  "ATGG",
  "CTAG"
)
frequencies <- calculate_allele_frequencies(sequences)
print(frequencies)

#######Task 3#######

compute_genotype_frequencies <- function(p, q) {
  # Compute genotype frequencies using Hardy-Weinberg equilibrium.
  
  # Check if p and q sum up to 1
  if (!(0 <= p && p <= 1) || !(0 <= q && q <= 1) || abs(p + q - 1) > 1e-6) {
    cat("Allele frequencies (p and q) should sum up to 1.\n")
    return(NULL)
  }
  
  AA <- p^2
  Aa <- 2 * p * q
  aa <- q^2
  
  return(list(AA = AA, Aa = Aa, aa = aa))
}

# Test the function
p <- 0.75  # Frequency of allele A
q <- 0.25  # Frequency of allele a

frequencies <- compute_genotype_frequencies(p, q)

cat("Genotype frequencies based on Hardy-Weinberg equilibrium:\n")
cat(sprintf("AA (homozygous dominant): %.2f\n", frequencies$AA))
cat(sprintf("Aa (heterozygous): %.2f\n", frequencies$Aa))
cat(sprintf("aa (homozygous recessive): %.2f\n", frequencies$aa))

#######Task 4#######

generate_dotplot <- function(seq1, seq2) {
  # Generate a dot plot for two sequences.
  
  # Initialize a blank matrix
  matrix <- matrix(0, nrow=nchar(seq1), ncol=nchar(seq2))
  
  # Fill the matrix: 1 indicates a match, 0 indicates no match
  for (i in 1:nchar(seq1)) {
    for (j in 1:nchar(seq2)) {
      if (substr(seq1, i, i) == substr(seq2, j, j)) {
        matrix[i, j] <- 1
      }
    }
  }
  
  return(matrix)
}

plot_dotplot <- function(matrix, seq1, seq2) {
  # Visualize the dot plot matrix using base R.
  
  image(t(matrix)[, nrow(matrix):1], col=c("white", "black"), xaxt='n', yaxt='n')
  axis(1, at=1:nchar(seq2), labels=strsplit(seq2, '')[[1]])
  axis(2, at=1:nchar(seq1), labels=strsplit(seq1, '')[[1]], las=2)
}

# Test the function
seq1 <- "TACGTACGGT"
seq2 <- "TACGATCGGT"

matrix <- generate_dotplot(seq1, seq2)
plot_dotplot(matrix, seq1, seq2)

#######Task 5#######

# Simulated reference genome
reference <- "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"

# Simulated sequencing reads from two individuals
individual1 <- "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"  # Same as reference
individual2 <- "AGCTAGCTAGCTAGCTTGCTAGCTAGCTAGCTAGCT"  # One SNP at position 16 (T instead of A)

call_variants <- function(reference, individual) {
  variants <- list()
  for (i in 1:nchar(reference)) {
    if (substr(reference, i, i) != substr(individual, i, i)) {
      # Position, reference nucleotide, individual's nucleotide
      variants <- append(variants, list(list(i, substr(reference, i, i), substr(individual, i, i))))
    }
  }
  return(variants)
}

# Variant calling
variants1 <- call_variants(reference, individual1)
variants2 <- call_variants(reference, individual2)

cat("Variants for Individual 1:", variants1, "\n")
cat("Variants for Individual 2:", variants2, "\n")

#######Task 6#######

library(Matrix)
library(stats)

# Parameters
num_individuals <- 1000
num_snps <- 100
phenotype_snp_index <- 42  # This SNP is associated with the phenotype

# Set seed for reproducibility
set.seed(0)

# Generate synthetic data
genotypes <- matrix(sample(0:1, num_individuals * num_snps, replace=TRUE), nrow=num_individuals)
phenotypes <- rep(0, num_individuals)

# The individuals with the associated SNP will have a higher chance of showing the phenotype
for (i in 1:num_individuals) {
  if (genotypes[i, phenotype_snp_index] == 1) {
    phenotypes[i] <- sample(0:1, 1, prob=c(0.4, 0.6))  # 60% chance to show phenotype
  } else {
    phenotypes[i] <- sample(0:1, 1, prob=c(0.7, 0.3))  # 30% chance to show phenotype
  }
}

# GWAS Testing
p_values <- numeric(num_snps)
for (i in 1:num_snps) {
  table <- matrix(0, nrow=2, ncol=2)
  for (j in 1:num_individuals) {
    table[genotypes[j, i] + 1, phenotypes[j] + 1] <- table[genotypes[j, i] + 1, phenotypes[j] + 1] + 1
  }
  p_values[i] <- chisq.test(table)$p.value
}

# Find the SNP with the most significant association
min_p_value_index <- which.min(p_values)
cat(sprintf("The SNP most associated with the phenotype is SNP %d with p-value %e", min_p_value_index, p_values[min_p_value_index]))
