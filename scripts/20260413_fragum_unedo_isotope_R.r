# Model fitting and plotting for the Fragum unedo isotope dataset.

library(lme4)
library(lmerTest)
library(MuMIn)
library(broom.mixed)

data_file <- file.path("data", "2026_Fragum_unedo_isotope.csv")
output_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

data <- read.csv(data_file)
data$size <- factor(data$size, levels = c("small", "medium", "big"))
data$tissue <- factor(data$tissue, levels = c("foot", "mantle", "gill"))

model_a <- lm(Delta_13C ~ height_mm, data = data)
model_b <- lm(Delta_13C ~ size, data = data)
model_c <- lm(Delta_13C ~ tissue, data = data)
model_d <- lmer(Delta_13C ~ height_mm + (1 | tissue), data = data, REML = FALSE)

print(r.squaredGLMM(model_a))
print(r.squaredGLMM(model_b))
print(r.squaredGLMM(model_c))
print(r.squaredGLMM(model_d))
print(AIC(model_a, model_b, model_c, model_d))

print(broom.mixed::tidy(model_a))
print(broom.mixed::tidy(model_b))
print(broom.mixed::tidy(model_c))
print(broom.mixed::tidy(model_d))

print(anova(model_a))
print(anova(model_b))
print(anova(model_c))
print(anova(model_d))
print(VarCorr(model_d)$tissue[1])

diagnostics_file <- file.path(output_dir, "20260413_fragum_unedo_isotope_diagnostics.pdf")
pdf(diagnostics_file, width = 10, height = 5)
par(mfrow = c(1, 2))
plot(model_a, which = 1, col = "#1d5c5199", pch = 18)
plot(model_a, which = 2, col = "#a64d7999", pch = 18)
dev.off()

boxplot_file <- file.path(output_dir, "20260413_fragum_unedo_isotope_size_boxplot.pdf")
pdf(boxplot_file, width = 7, height = 5)
boxplot(
  Delta_13C ~ size,
  data = data,
  main = "Carbon Isotopes by Sample Category",
  xlab = "Size Category",
  ylab = expression(paste(delta^{13}, "C (\u2030)")),
  col = c("#1d5c5199", "#a64d7999", "#96734d99"),
  ylim = c(-22, -15)
)
stripchart(
  Delta_13C ~ size,
  data = data,
  vertical = TRUE,
  method = "jitter",
  add = TRUE,
  pch = 18,
  cex = 2,
  col = "#00000088"
)
dev.off()

scatter_file <- file.path(output_dir, "20260413_fragum_unedo_isotope_height_scatter.pdf")
pdf(scatter_file, width = 7, height = 5)
plot_colors <- c("#E35336", "#FFD3AC", "#9988A1")[data$size]
plot_shapes <- c(8, 18, 21)[data$tissue]
plot(
  data$height_mm,
  data$Delta_13C,
  main = "Carbon Isotopes vs Sample Height",
  xlab = "Shell Height (mm)",
  ylab = expression(paste(delta^{13}, "C (\u2030)")),
  cex = 2,
  col = plot_colors,
  pch = plot_shapes,
  xlim = c(0, 40),
  xaxt = "n",
  ylim = c(-22, -15)
)
axis(side = 1, at = seq(0, 40, by = 5))
abline(model_a, col = "black", lwd = 2, lty = 3)
legend("bottomright", legend = levels(data$tissue), pch = c(8, 18, 21), col = "black", bty = "n")
intercept <- round(coef(model_a)[1], 2)
slope <- round(coef(model_a)[2], 2)
r_squared <- round(summary(model_a)$r.squared, 2)
label <- bquote(R^2 == .(r_squared) ~ ";" ~ y == .(slope) * x - .(abs(intercept)))
text(x = 26, y = -19.5, labels = label, pos = 4, cex = 1.1)
dev.off()
