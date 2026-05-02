# Figure 1: broad isotope metadata plot and labeled supplement.

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggrepel)

data_file <- file.path("data", "Isotope_metadata.tsv")
main_file <- file.path("figures", "20260420_isotope_metadata_main.pdf")
supplement_file <- file.path("figures", "20260420_isotope_metadata_supplement_labeled.pdf")
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

symbiosis_colors <- c(
  "Chemosymbiotic" = "#2F6BFF",
  "Non-symbiotic" = "#8A94A6",
  "Photosymbiotic" = "#27A36B"
)

read_and_clean_metadata <- function(path) {
  read_tsv(path, col_types = cols(.default = "c")) %>%
    mutate(
      Delta_13C_average = suppressWarnings(as.numeric(Delta_13C_average)),
      Delta_13C_SD = suppressWarnings(as.numeric(Delta_13C_SD)),
      Species = str_replace_all(coalesce(Species, ""), "[^A-Za-z]", ""),
      Species = na_if(str_trim(Species), ""),
      Tissue = na_if(str_trim(coalesce(Tissue, "")), ""),
      Life_Stage = na_if(str_trim(coalesce(Life_Stage, "")), ""),
      Category = if_else(Phylum == "Mollusca", "Mollusca", "Cnidaria"),
      Category = factor(Category, levels = c("Mollusca", "Cnidaria")),
      Symbiosis = recode(
        Symbiotic,
        "Y" = "Photosymbiotic",
        "N" = "Non-symbiotic",
        "C" = "Chemosymbiotic",
        .default = "Unknown"
      ),
      label = case_when(
        !is.na(Genus) & !is.na(Species) & !is.na(Life_Stage) & !is.na(Tissue) ~
          paste(Genus, Species, Life_Stage, Tissue, sep = " | "),
        !is.na(Genus) & !is.na(Species) & !is.na(Tissue) ~
          paste(Genus, Species, Tissue, sep = " | "),
        !is.na(Genus) & !is.na(Species) ~ paste(Genus, Species),
        !is.na(Genus) ~ Genus,
        TRUE ~ "Unspecified"
      )
    ) %>%
    filter(!is.na(Delta_13C_average)) %>%
    filter(Phylum != "Porifera") %>%
    filter(!(Genus %in% c("Pinguitellina", "Circe"))) %>%
    arrange(Category, Symbiosis, Delta_13C_average) %>%
    mutate(
      x_jitter = c(
        seq(0.72, 1.92, length.out = sum(Category == "Mollusca")),
        seq(2.18, 2.46, length.out = sum(Category == "Cnidaria"))
      )
    )
}

plot_metadata <- function(data, show_labels = FALSE) {
  p <- ggplot(data, aes(x = x_jitter, y = Delta_13C_average, color = Symbiosis)) +
    geom_hline(yintercept = seq(-40, -10, by = 5), color = "#E9E1D5", linewidth = 0.3) +
    geom_errorbar(
      aes(ymin = Delta_13C_average - Delta_13C_SD, ymax = Delta_13C_average + Delta_13C_SD),
      width = 0,
      linewidth = 0.4,
      alpha = 0.55,
      na.rm = TRUE
    ) +
    geom_point(size = 3.1, stroke = 0, alpha = 0.9) +
    scale_x_continuous(
      limits = c(0.55, 2.62),
      breaks = c(1.32, 2.32),
      labels = c("Mollusca", "Cnidaria"),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    scale_y_continuous(
      limits = c(-40, -10),
      breaks = seq(-40, -10, by = 5),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_color_manual(values = symbiosis_colors, drop = FALSE) +
    labs(x = NULL, y = expression(paste(delta^{13}, "C (\u2030)")), color = NULL) +
    theme_classic(base_size = 15) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.line = element_line(color = "#2B2624", linewidth = 0.5),
      axis.ticks = element_line(color = "#2B2624", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      axis.title.y = element_text(size = 18, margin = margin(r = 12), color = "#2B2624"),
      axis.text = element_text(size = 14, color = "#2B2624"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = element_text(size = 13),
      legend.title = element_blank(),
      plot.margin = margin(12, 18, 12, 12)
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(shape = 16, size = 4, alpha = 1)))

  if (show_labels) {
    p <- p +
      geom_text_repel(
        aes(label = label),
        size = 3,
        color = "#2B2624",
        box.padding = 0.55,
        point.padding = 0.3,
        segment.color = "#8E8376",
        segment.alpha = 0.7,
        seed = 42,
        max.overlaps = Inf,
        min.segment.length = 0,
        force = 3.5,
        force_pull = 0.2,
        max.time = 15,
        max.iter = 50000
      )
  }

  p
}

plot_data <- read_and_clean_metadata(data_file)

pdf(main_file, width = 10.8, height = 6.8)
print(plot_metadata(plot_data, show_labels = FALSE))
dev.off()

pdf(supplement_file, width = 14.5, height = 9.5)
print(plot_metadata(plot_data, show_labels = TRUE))
dev.off()
