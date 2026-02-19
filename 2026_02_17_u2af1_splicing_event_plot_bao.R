library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Load your rMATS Skipped Exon (SE) data
df1 <- read.csv("mo_q157r_3s_new_bao.csv", header = TRUE, stringsAsFactors = FALSE)


# ---------- Compute PSI ----------
df <- df %>%
  mutate(
    PSI_WT = calc_mean_psi(WT),
    PSI_Q157R    = calc_mean_psi(KI)
  )

# ---------- Define significant events (REVISED) ----------
plot_data <- df %>%
  mutate(
    Status = case_when(
      # If WT > Q157R, the mutant causes skipping
      FDR <= 0.05 & IncLevelDifference >= 0.1  ~ "Exon Skipping", 
      
      # If Q157R > WT, the mutant causes inclusion
      FDR <= 0.05 & IncLevelDifference <= -0.1 ~ "Exon Inclusion", 
      
      TRUE ~ "Not Significant"
    )
  )

plot_sig <- plot_data %>%
  filter(Status != "Not Significant")

# ---------- Counts for annotation ----------
n_skip <- sum(plot_sig$Status == "Exon Skipping", na.rm = TRUE)
n_inc  <- sum(plot_sig$Status == "Exon Inclusion", na.rm = TRUE)

# ---------- Manuscript-quality scatter ----------
p <- ggplot(plot_sig, aes(x = PSI_Q157R, y = PSI_WT)) +
  geom_point(aes(color = Status), size = 2, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  scale_color_manual(values = c(
    "Exon Skipping"  = "#D7191C",
    "Exon Inclusion" = "#2731F5"
  )) +
  labs(
    x = expression(U2AF1^{Q157R} ~ "PSI"),
    y = "WT PSI",
    color = NULL
  ) +
  annotate("text", x = 0.05, y = 0.95,
           label = paste0("Exon Skipping (n = ", n_skip, ")"),
           hjust = 0, vjust = 1, size = 4, fontface = "bold") +
  annotate("text", x = 0.95, y = 0.05,
           label = paste0("Exon Inclusion (n = ", n_inc, ")"),
           hjust = 1, vjust = 0, size = 4, fontface = "bold") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.line = element_line(color = "black"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8)
  )

p


