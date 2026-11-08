#Package setup
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

#Data output directory setup
args <- commandArgs(trailingOnly = TRUE); out_dir <- if (length(args)) args[1] else "out"
f <- readRDS(file.path(out_dir, "prediction_rows.rds"))

# Sets the figure colours and bar spacing
pal <- c("NPi 4 or higher" = "#86b6ef", "NPi 3 to 3.9" = "#2a78d6", "NPi below 3 (abnormal)" = "#0d366b")
ink <- "#1f2328"; muted <- "#5b6168"; dg <- position_dodge(width = 0.8)

# Calculates observed 24 h risk for each NPi band within each admission GCS group
d <- f %>%
  mutate(npi = factor(cut(npi_min_c, c(-Inf, 3, 4, Inf), right = FALSE, labels = rev(names(pal))), levels = names(pal)),
         gcs = cut(admit_gcs, c(2, 8, 12, 15), labels = c("GCS 3-8", "GCS 9-12", "GCS 13-15"))) %>%
  group_by(gcs, npi) %>%
  summarise(risk = mean(y == "yes"), checks = n(), pts = n_distinct(patient_id), .groups = "drop")

# Draws the bar chart with labels, title and caption
p <- ggplot(d, aes(gcs, risk, fill = npi)) +
  geom_col(position = dg, width = 0.74, colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * risk)), position = dg, vjust = -0.45, size = 4.2, colour = ink, fontface = "bold") +
  scale_fill_manual(values = pal, name = "Current worst-eye NPi") +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%"), limits = c(0, 0.62),
                     breaks = seq(0, 0.6, 0.1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "A low NPi signals deterioration in the next 24 hours, whatever the admission GCS",
       subtitle = "Share of pupil checks followed by neurologic deterioration within 24 h",
       x = "Admission Glasgow Coma Scale", y = NULL,
       caption = paste0("Simulated Neuro ICU cohort: 592 patients, 23,793 pupil checks before any first deterioration. ",
                        "NPi = Neurological Pupil index, worse of the two eyes,\nlast valid value. ",
                        "Deterioration = GCS drop of 2 or more, new pupil asymmetry, or emergent osmotherapy or surgery. ",
                        "Observed proportions, no model.")) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", colour = ink, size = 15.5),
        plot.subtitle = element_text(colour = muted, margin = margin(b = 10)),
        plot.caption = element_text(colour = muted, size = 8.8, hjust = 0, margin = margin(t = 12)),
        plot.title.position = "plot", plot.caption.position = "plot",
        axis.text = element_text(colour = ink, size = 12), axis.title.x = element_text(colour = muted, margin = margin(t = 8)),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "#e6e6e3", linewidth = 0.4),
        legend.position = "top", legend.justification = "left", legend.title = element_text(colour = muted, size = 11),
        legend.text = element_text(colour = ink, size = 11), plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(16, 20, 12, 16))

# Saves the figure as PNG and PDF and prints the numbers behind it
for (ext in c("png", "pdf")) ggsave(file.path(out_dir, paste0("figure_npi_24h_risk.", ext)), p, width = 10, height = 6.4, dpi = 200)
print(d)
