# Figure 2: cardiid-focused isotope comparison plot.

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

data_file <- file.path("data", "Isotope_metadata.tsv")
wa_file <- file.path("data", "20250205_WA_isotope_data_long.csv")
output_file <- file.path("figures", "20260417_metadata_alt3_cardiid_focus.pdf")
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

symbiosis_colors <- c(
  "Chemosymbiotic" = "#2F6BFF",
  "Non-symbiotic" = "#8A94A6",
  "Photosymbiotic" = "#27A36B"
)

metadata <- read_tsv(data_file, col_types = cols(.default = "c")) %>%
  mutate(
    Delta_13C_average = suppressWarnings(as.numeric(Delta_13C_average)),
    Delta_13C_SD = suppressWarnings(as.numeric(Delta_13C_SD)),
    Species = str_replace_all(coalesce(Species, ""), "[^A-Za-z]", ""),
    Species = na_if(str_trim(Species), ""),
    Tissue = na_if(str_trim(coalesce(Tissue, "")), ""),
    Category = if_else(Phylum == "Mollusca", "Mollusca", "Cnidaria"),
    Category = factor(Category, levels = c("Cnidaria", "Mollusca")),
    Symbiosis = recode(
      Symbiotic,
      "Y" = "Photosymbiotic",
      "N" = "Non-symbiotic",
      "C" = "Chemosymbiotic",
      .default = "Unknown"
    )
  ) %>%
  filter(!is.na(Delta_13C_average)) %>%
  filter(Phylum != "Porifera") %>%
  filter(!(Genus %in% c("Pinguitellina", "Circe")))

wa_fragum_unedo <- read_csv(wa_file, show_col_types = FALSE) %>%
  filter(
    Genus == "Fragum",
    species == "unedo",
    Data == "Delta_13C",
    tissue %in% c("Foot", "Mantle", "Gill")
  ) %>%
  mutate(
    Tissue = tissue,
    delta13c = as.numeric(Value),
    Symbiosis = recode(
      Symbiosis,
      "Y" = "Photosymbiotic",
      "N" = "Non-symbiotic",
      "C" = "Chemosymbiotic",
      .default = "Unknown"
    )
  ) %>%
  group_by(Symbiosis, Tissue) %>%
  summarise(
    Delta_13C_average = mean(delta13c, na.rm = TRUE),
    Delta_13C_SD = sd(delta13c, na.rm = TRUE),
    Genus = "Fragum",
    Species = "unedo",
    Family = "Cardiidae",
    Category = "Mollusca",
    .groups = "drop"
  )

plot_data <- metadata %>%
  filter((Genus == "Fragum" | Family == "Cardiidae") & !(Genus == "Fragum" & Species == "unedo")) %>%
  bind_rows(wa_fragum_unedo) %>%
  mutate(
    taxon = case_when(
      Genus == "Fragum" & Species == "unedo" ~ "Fragum unedo",
      Genus == "Fragum" & Species == "scruposum" ~ "Fragum whitleyi",
      !is.na(Species) ~ paste(Genus, Species),
      TRUE ~ Genus
    ),
    tissue_group = case_when(
      str_to_lower(Tissue) == "foot" ~ "foot",
      str_to_lower(Tissue) == "mantle" ~ "mantle",
      str_to_lower(Tissue) == "gill" ~ "gill",
      TRUE ~ "other"
    )
  ) %>%
  group_by(taxon, Symbiosis, tissue_group) %>%
  summarise(
    Delta_13C_average = mean(Delta_13C_average, na.rm = TRUE),
    Delta_13C_SD = if (all(is.na(Delta_13C_SD))) NA_real_ else mean(Delta_13C_SD, na.rm = TRUE),
    Genus = first(Genus),
    Species = first(Species),
    Family = first(Family),
    Tissue = first(Tissue),
    .groups = "drop"
  ) %>%
  group_by(taxon) %>%
  mutate(
    species_mean = mean(Delta_13C_average, na.rm = TRUE),
    tissue_group = factor(tissue_group, levels = c("foot", "mantle", "gill", "other")),
    tissue_offset = case_when(
      n_distinct(tissue_group[!is.na(tissue_group)]) <= 1 ~ 0,
      taxon == "Fragum unedo" & tissue_group == "foot" ~ -0.28,
      taxon == "Fragum unedo" & tissue_group == "mantle" ~ 0,
      taxon == "Fragum unedo" & tissue_group == "gill" ~ 0.28,
      taxon == "Fragum unedo" & tissue_group == "other" ~ 0.42,
      tissue_group == "foot" ~ -0.18,
      tissue_group == "mantle" ~ 0,
      tissue_group == "gill" ~ 0.18,
      TRUE ~ 0.32
    ),
    replicate_offset = if (n() > 1) seq(-0.05, 0.05, length.out = n()) else rep(0, n())
  ) %>%
  ungroup() %>%
  mutate(
    taxon = factor(taxon, levels = rev(unique(taxon[order(species_mean)]))),
    taxon_base = as.numeric(taxon),
    y_pos = taxon_base + tissue_offset + replicate_offset,
    taxon_label = paste0("italic('", as.character(taxon), "')")
  )

plot_obj <- ggplot(plot_data, aes(x = Delta_13C_average, y = y_pos, color = Symbiosis)) +
  geom_errorbar(
    aes(xmin = Delta_13C_average - Delta_13C_SD, xmax = Delta_13C_average + Delta_13C_SD),
    orientation = "y",
    width = 0.16,
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  geom_point(aes(shape = tissue_group), size = 3.2, alpha = 0.95, stroke = 0.9, fill = "white") +
  scale_y_continuous(
    breaks = seq_along(levels(plot_data$taxon)),
    labels = function(x) {
      label_tbl <- unique(plot_data[order(plot_data$taxon_base), c("taxon_base", "taxon_label")])
      parse(text = label_tbl$taxon_label[match(x, label_tbl$taxon_base)])
    },
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_color_manual(values = symbiosis_colors, drop = FALSE) +
  scale_shape_manual(
    values = c("foot" = 8, "mantle" = 18, "gill" = 21, "other" = 17),
    breaks = c("foot", "mantle", "gill", "other"),
    labels = c("foot", "mantle", "gill", "other")
  ) +
  labs(
    x = expression(paste(delta^{13}, "C (\u2030)")),
    y = NULL,
    color = NULL,
    shape = "Tissue"
  ) +
  theme_classic(base_size = 15) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "#2B2624", linewidth = 0.5),
    axis.ticks = element_line(color = "#2B2624", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title.x = element_text(size = 18, margin = margin(t = 10), color = "#2B2624"),
    axis.text.x = element_text(size = 14, color = "#2B2624"),
    axis.text.y = element_text(size = 13, color = "#2B2624"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.margin = margin(12, 18, 12, 12)
  ) +
  guides(
    color = guide_legend(order = 1, nrow = 1, byrow = TRUE, override.aes = list(shape = 16, size = 4, alpha = 1)),
    shape = guide_legend(order = 2, nrow = 1, byrow = TRUE)
  )

pdf(output_file, width = 9.8, height = 7.2)
print(plot_obj)
dev.off()
