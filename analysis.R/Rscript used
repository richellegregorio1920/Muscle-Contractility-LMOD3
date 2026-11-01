
# 1. Load Required Libraries
suppressPackageStartupMessages({
  library(rentrez)
  library(ggplot2)
})

# 2. Establish Directory Structure
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

cat("[INFO] Initializing LMOD3 computational workflow...\n")

# 3. Retrieve Data from NCBI via rentrez
# Gene ID for human LMOD3 is 348932. RefSeq accession: NM_198271.5 -> NP_938012.2
gene_id <- "348932"
transcript_id <- "NM_198271.5"
protein_id <- "NP_938012.2"

cat(sprintf("[INFO] Fetching records for gene %s and protein %s...\n", gene_id, protein_id))

# Fetch protein sequence in FASTA format
prot_fetch <- entrez_fetch(db = "protein", id = protein_id, rettype = "fasta")

# Parse FASTA string into header and sequence lines
fasta_lines <- strsplit(prot_fetch, "\n")[[1]]
seq_header <- fasta_lines[1]
protein_seq <- paste(fasta_lines[-1], collapse = "")
protein_seq <- gsub("\\s+", "", protein_seq) # remove any stray whitespace

# Basic Metrics
total_length <- nchar(protein_seq)
cat(sprintf("[SUCCESS] Retrieved protein sequence successfully. Total length: %d amino acids.\n", total_length))

# 4. Amino Acid Composition Analysis
aa_vector <- strsplit(protein_seq, "")[[1]]
aa_counts <- table(aa_vector)
aa_df <- as.data.frame(aa_counts)
names(aa_df) <- c("AminoAcid", "Count")
aa_df$Percentage <- (aa_df$Count / total_length) * 100

# Save composition results to CSV
write.csv(aa_df, "results/amino_acid_composition.csv", row.names = FALSE)
cat("[INFO] Saved amino acid composition table to results/amino_acid_composition.csv\n")

# 5. Generate Figure 1: Amino-Acid Composition Bar Plot
p1 <- ggplot(aa_df, aes(x = reorder(AminoAcid, -Count), y = Count, fill = Count)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.85) +
  scale_fill_viridis_c(option = "cividis") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Figure 1: Amino-Acid Composition of Human Lmod3",
    subtitle = sprintf("RefSeq Protein: %s (Total Length: %d residues)", protein_id, total_length),
    x = "Amino Acid Residue",
    y = "Absolute Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

ggsave("figures/01_amino_acid_composition.png", plot = p1, width = 8, height = 5, dpi = 300)
cat("[INFO] Generated and saved Figure 1: figures/01_amino_acid_composition.png\n")

# 6. Generate Figure 2: Linear Variant Position Map
# Variant: p.Leu550Phe (occurs at amino acid position 550 out of 560)
variant_pos <- 550
variant_label <- "p.Leu550Phe (ClinVar ID: 374498)"

p2 <- ggplot() +
  # Protein backbone track
  geom_segment(aes(x = 1, xend = total_length, y = 1, yend = 1), linewidth = 4, color = "grey60") +
  # Variant annotation point
  geom_point(aes(x = variant_pos, y = 1), color = "firebrick", size = 5) +
  geom_text(aes(x = variant_pos, y = 1.25, label = variant_label), color = "firebrick", fontface = "bold", size = 3.5) +
  ylim(0.5, 1.6) +
  xlim(0, total_length + 20) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Figure 2: Linear Variant Position Map for LMOD3",
    subtitle = sprintf("Assigned Variant c.1648C>T (%s) mapped on %s backbone", "p.Leu550Phe", protein_id),
    x = "Amino Acid Position (Codon Coordinate)",
    y = ""
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("figures/02_variant_position.png", plot = p2, width = 8, height = 3.5, dpi = 300)
cat("[INFO] Generated and saved Figure 2: figures/02_variant_position.png\n")

cat("[SUCCESS] All analysis steps, data outputs, and figure generations completed successfully!\n")
