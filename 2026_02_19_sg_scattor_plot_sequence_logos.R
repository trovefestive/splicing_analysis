# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(Biostrings)
library(seqLogo)
library(BSgenome.Mmusculus.UCSC.mm10) # Change to Hsapiens if human

# ==========================================
# 1. DATA LOADING & PRE-PROCESSING
# ==========================================

# Robust PSI calculation function to handle "NA" strings and commas
calc_mean_psi <- function(psi_col) {
  sapply(strsplit(as.character(psi_col), ","), function(x) {
    nums <- suppressWarnings(as.numeric(x))
    if (all(is.na(nums))) return(NA)
    mean(nums, na.rm = TRUE)
  })
}

# Load rMATS JCEC output
df <- read.table("SE.MATS.JCEC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Clean and compute means
df_processed <- df %>%
  mutate(
    PSI_WT    = calc_mean_psi(IncLevel1),
    PSI_Q157R = calc_mean_psi(IncLevel2),
    # Filter out non-standard chromosomes to avoid Biostrings errors later
    StandardChr = grepl("^chr[0-9XYM]+$", chr)
  ) %>%
  filter(StandardChr == TRUE)

# Define significance based on Mutant (Q157R) behavior relative to WT
# If Diff > 0: WT > Q157R -> Exon is SKIPPED by mutation
# If Diff < 0: Q157R > WT -> Exon is INCLUDED by mutation
plot_data <- df_processed %>%
  mutate(
    Status = case_when(
      FDR <= 0.05 & IncLevelDifference >= 0.1  ~ "Exon Skipping",
      FDR <= 0.05 & IncLevelDifference <= -0.1 ~ "Exon Inclusion",
      TRUE ~ "Not Significant"
    )
  )

# ==========================================
# 2. GENERATE SCATTER PLOT
# ==========================================

plot_sig <- plot_data %>% filter(Status != "Not Significant")
n_skip <- sum(plot_sig$Status == "Exon Skipping")
n_inc  <- sum(plot_sig$Status == "Exon Inclusion")

scatter_p <- ggplot(plot_sig, aes(x = PSI_WT, y = PSI_Q157R)) +
  geom_point(aes(color = Status), size = 1.5, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Exon Skipping" = "#D7191C", "Exon Inclusion" = "#2731F5")) +
  labs(x = "WT PSI", y = expression(U2AF1^{Q157R} ~ "PSI"), title = "Differential Splicing (SE)") +
  annotate("text", x = 0.05, y = 0.95, label = paste0("Skipping: ", n_skip), hjust = 0, fontface = "bold") +
  annotate("text", x = 0.95, y = 0.05, label = paste0("Inclusion: ", n_inc), hjust = 1, fontface = "bold") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", panel.border = element_rect(fill = NA, color = "black"))

print(scatter_p)

# ==========================================
# 3. SEQUENCE LOGO ANALYSIS (3' Splice Site)
# ==========================================

# Function to extract 3'ss (-20nt intron to +3nt exon)
get_3ss_seqs <- function(sub_df) {
  # pos is the junction: exonStart_0base for +, exonEnd for -
  starts <- ifelse(sub_df$strand == "+", sub_df$exonStart_0base - 19, sub_df$exonEnd - 2)
  ends   <- ifelse(sub_df$strand == "+", sub_df$exonStart_0base + 3, sub_df$exonEnd + 20)
  
  getSeq(BSgenome.Mmusculus.UCSC.mm10, names = sub_df$chr, 
         start = starts, end = ends, strand = sub_df$strand)
}

# Extract sequences for the two groups
seqs_skip <- get_3ss_seqs(subset(plot_data, Status == "Exon Skipping"))
seqs_inc  <- get_3ss_seqs(subset(plot_data, Status == "Exon Inclusion"))

# Create PWMs
pwm_skip <- makePWM(consensusMatrix(seqs_skip, as.prob = TRUE)[1:4,])
pwm_inc  <- makePWM(consensusMatrix(seqs_inc, as.prob = TRUE)[1:4,])

# Plot Logos
# Position 21 is the +1 position of the exon
par(mfrow = c(2, 1))
seqLogo(pwm_skip)
grid::grid.text("Exons Skipped in Q157R (WT Preference)", x = 0.5, y = 0.95, gp = grid::gpar(fontsize=10, font=2))
seqLogo(pwm_inc)
grid::grid.text("Exons Included in Q157R (Mutant Preference)", x = 0.5, y = 0.45, gp = grid::gpar(fontsize=10, font=2))

# ==========================================
# 4. SAVE OUTPUTS
# ==========================================
# ggsave("Scatter_WT_vs_Q157R.pdf", scatter_p, width = 5, height = 5)