library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Load the summary data
# Ensuring we read the tab-separated format from rMATS
summary_df <- read.table("2026_02_19_summary_sg_analysis.txt", header = TRUE, sep = "\t")

# 2. Reshape data for plotting
# We want to compare Sample 1 Higher vs Sample 2 Higher for JCEC
plot_data <- summary_df %>%
  select(EventType, 
         "WT Higher (Skipping)" = SigEventsJCECSample1HigherInclusion, 
         "Q157R Higher (Inclusion)" = SigEventsJCECSample2HigherInclusion) %>%
  pivot_longer(cols = -EventType, names_to = "Direction", values_to = "Count")

# 3. Create the Manuscript-Ready Bar Plot
p_summary <- ggplot(plot_data, aes(x = EventType, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.8), 
            vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("WT Higher (Skipping)" = "#D7191C", 
                               "Q157R Higher (Inclusion)" = "#2731F5")) +
  labs(
    title = "Significant Splicing Events (JCEC)",
    subtitle = "FDR < 0.05, |delta PSI| > 0.1",
    x = "Alternative Splicing Event Type",
    y = "Number of Significant Events",
    fill = "Direction of Change"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90")
  )

# Display plot
print(p_summary)

# 4. Optional: Save for publication
# ggsave("rMATS_Summary_Barplot.pdf", p_summary, width = 7, height = 5)