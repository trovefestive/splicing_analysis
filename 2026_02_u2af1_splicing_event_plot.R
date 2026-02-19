library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Load your rMATS Skipped Exon (SE) data
df <- read.table("SE.MATS.JCEC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 2. Function to calculate the mean PSI from comma-separated strings
calc_mean_psi <- function(psi_str) {
  sapply(strsplit(as.character(psi_str), ","), function(x) {
    val <- as.numeric(x)
    if (all(is.na(val))) return(NA)
    mean(val, na.rm = TRUE)
  })
}

# 3. Process the data for plotting
plot_sig <- plot_data %>%
  filter(Status != "Not Significant")


ggplot(plot_sig) +
  geom_point(
    aes(x = PSI_Q157R, y = PSI_WT, color = Status),
    size = 1.5
  ) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  scale_color_manual(values = c(
    "Exon Skipping"  = "blue",
    "Exon Inclusion" = "red"
  )) +
  labs(
    x = expression("U2AF1"^"Q157R"),
    y = "WT"
  ) +
  annotate(
    "text", x = 0.15, y = 0.92,
    label = paste0("Exon Skipping\nn=", sum(plot_sig$Status == "Exon Skipping")),
    color = "blue", fontface = "bold", hjust = 0
  ) +
  annotate(
    "text", x = 0.95, y = 0.08,
    label = paste0("Exon Inclusion\nn=", sum(plot_sig$Status == "Exon Inclusion")),
    color = "red", fontface = "bold", hjust = 1
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(color = "black"),
    panel.border = element_rect(fill = NA, color = "black")
  )



## unlabeled plot
ggplot(plot_sig) +
  geom_point(
    aes(x = PSI_Q157R, y = PSI_WT, color = Status),
    size = 1.5
  ) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  scale_color_manual(values = c(
    "Exon Skipping"  = "blue",
    "Exon Inclusion" = "red"
  )) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    axis.line = element_line(color = "black"),
    panel.border = element_rect(fill = NA, color = "black")
  )

