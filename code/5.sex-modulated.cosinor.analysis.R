library(lmerTest)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(qvalue)
library(CircStats)

load("data/Human_annual_metabolomics.RData")

source("code/supporting.functions.R")

# parameter definition
period <- 365

w <- 2 * pi / period

seasonal.metabolites <- readRDS("results/human_amide_lmm_cosinor_annual.rds")
rownames(seasonal.metabolites) <- seasonal.metabolites$unique_id

SOMs_amide <- subset(seasonal.metabolites, q.value < 0.05)


################################################################################
#   metabolome data fitting with sex difference #
################################################################################


features <- human_amide_BGI[, sample_metadata$Sample]
rownames(features) <- human_amide_BGI$unique_id

## outlier removal
norm.data <- apply(features, 1, function(x) log2(x) |> scale()) |> t()
sum(abs(norm.data) > 3)
# 10127
features[abs(norm.data) > 3] <- NA

PCA <- prcomp(norm.data)


PCA.m <- prcomp(norm.data[, sample_metadata$Sex == "male"])

PCA.f <- prcomp(norm.data[, sample_metadata$Sex == "female"])


ggplot(data = subset(sample_metadata, Sex == "male")) +
  geom_point(aes(
    x = PCA.m$rotation[, 1], y = PCA.m$rotation[, 2],
    color = Subject
  ), size = 3) +
  scale_colour_manual(values = subject.color(1:6)) +
  theme_classic() +
  theme(legend.position = "bottom", validate = TRUE) +
  ggplot(data = subset(sample_metadata, Sex == "female")) +
  geom_point(aes(
    x = PCA.f$rotation[, 1], y = PCA.f$rotation[, 2],
    color = Subject
  ), size = 3) +
  scale_colour_manual(values = subject.color(7:11)) +
  theme_classic() +
  theme(legend.position = "bottom", validate = TRUE)

ggplot(data = sample_metadata) +
  geom_point(aes(
    x = PCA$rotation[, 1], y = PCA$rotation[, 2],
    color = Subject
  ), size = 3) +
  scale_x_continuous(name = "PC1") +
  scale_y_continuous(name = "PC2") +
  scale_colour_manual(values = subject.color(1:11)) +
  theme_classic() +
  theme(legend.position = "bottom", validate = TRUE)


# fitting
lmm_cosinor.sm.result <- data.frame(
  row.names = human_amide_BGI$unique_id
)


for (m in rownames(features)) {
  y <- features[m, ] |>
    log2() |>
    unlist()
  suppressMessages({
    fit.cosinor.sm <- lmer(y ~ (sin.a + cos.a) * Sex + (1 | Subject), data = sample_metadata, REML = F)
  })

  para <- fixef(fit.cosinor.sm)
  para["p.value"] <- contest(fit.cosinor.sm, diag(6)[c(2:3, 5:6), ])[1, "Pr(>F)"]

  lmm_cosinor.sm.result[m, names(para)] <- para
}

## FDR correction

lmm_cosinor.sm.result$q.value <- qvalue(lmm_cosinor.sm.result$p.value)$qvalues


## parameter calculation

lmm_cosinor.sm.result$amplitude.male <-
  with(lmm_cosinor.sm.result, amplitude(cos.a, sin.a) * 2)

lmm_cosinor.sm.result$phase.male <-
  with(lmm_cosinor.sm.result, phase(cos.a, sin.a))

lmm_cosinor.sm.result$mesor.male <-
  with(lmm_cosinor.sm.result, `(Intercept)`)

## parameter calculation
lmm_cosinor.sm.result$amplitude.female <-
  with(lmm_cosinor.sm.result, amplitude(cos.a + `cos.a:Sexfemale`, sin.a + `sin.a:Sexfemale`) * 2)

lmm_cosinor.sm.result$phase.female <-
  with(lmm_cosinor.sm.result, phase(cos.a + `cos.a:Sexfemale`, sin.a + `sin.a:Sexfemale`))

lmm_cosinor.sm.result$mesor.female <-
  with(lmm_cosinor.sm.result, `(Intercept)` + Sexfemale)


# sort by FDR
lmm_cosinor.sm.result <- lmm_cosinor.sm.result[seasonal.metabolites$unique_id, ]

# export parameter table
merge(human_amide_BGI[, 1:13], lmm_cosinor.sm.result, by.x = "unique_id", by.y = 0) |>
  write.csv("results/human_metabolome_amide_mixed-effects-cosinor_sex-modulated-annual.csv", na = "", row.names = F)

saveRDS(
  merge(human_amide_BGI[, 1:13], lmm_cosinor.sm.result, by.x = "unique_id", by.y = 0),
  "results/human_amide_lmm_cosinor_sex-modulated-annual.rds"
)

sum(lmm_cosinor.sm.result$q.value < 0.05)




sm.seasonal.metabolites <- readRDS("results/human_amide_lmm_cosinor_sex-modulated-annual.rds")
rownames(sm.seasonal.metabolites) <- sm.seasonal.metabolites$unique_id


q.cutoff <- 0.05


smSOMs_amide <- subset(sm.seasonal.metabolites, q.value < q.cutoff)
smSOMs_amide <- smSOMs_amide[order(smSOMs_amide$phase.male), ]

# sex differences
################################################################################

ids <- c(
  "Amide_POS_6.028_243.05294", # trp
  "Amide_NEG_7.138_214.04883", # Glycerophosphoethanolamine
  "Amide_POS_1.036_385.30602", # N−palmitoyl glutamin
  "Amide_NEG_2.165_207.0511", # hydroxy adipic acid
  "Amide_POS_2.108_220.11827", # pantothenic acid
  "Amide_NEG_2.662_309.20722" # EpKODE
)

smSOM.representatives <- smSOMs_amide[ids, ]

smSOM.representatives$curated_name <-
  c(
    "tryptophan",
    "glycerophosphoethanolamine",
    "N−palmitoyl glutamine",
    "2-hydroxyadipic acid",
    "pantothenic acid",
    "12(13)Ep-9-KODE"
  )

representative.sm.plots <- plot_metabolite.sm(smSOM.representatives, centered = T)

representative.sm.plots[[2]] <- representative.sm.plots[[2]] +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
representative.sm.plots[[4]] <- representative.sm.plots[[4]] +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
representative.sm.plots[[6]] <- representative.sm.plots[[6]] +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


pdf("figures/sex_difference/representative_smSOMs.pdf", width = 9, height = 4)
representative.sm.plots |>
  wrap_plots(ncol = 3, byrow = F)
dev.off()


representative.sm.sep.plots <- plot_metabolite.sm.sep(smSOM.representatives, centered = T)
pdf("figures/sex_difference/representative_sdSOMs_sep.pdf", width = 18, height = 4)
representative.sm.sep.plots |>
  wrap_plots(ncol = 6, byrow = T)
dev.off()

bins <- 25

count.m <- circular.SW.count(smSOMs_amide$phase.male, bins = bins, period = period)
count.f <- circular.SW.count(smSOMs_amide$phase.female, bins = bins, period = period)

sum(count.m$count)
sum(count.f$count)

scale.seq <- seq(10, 70, 20)

rose_male <- ggplot(data = count.m) +
  geom_hline(yintercept = scale.seq |> sqrt(), color = "grey80") +
  geom_vline(xintercept = c(0), color = "grey20") +
  geom_bar(aes(x = center * w, y = sqrt(count)),
    stat = "identity",
    fill = circular.phase.color(seq(from = 0, to = period, length.out = bins))
  ) +
  scale_y_continuous(
    name = "number of metabolites", limit = c(0, 9),
    breaks = scale.seq |> sqrt(), labels = scale.seq
  ) +
  coord_radial(expand = F, start = 0, end = 2 * pi, r.axis.inside = F, clip = "off") +
  scale_x_continuous(
    name = "phase (day)",
    limits = c(0, 2 * pi),
    breaks = first.month * w,
    labels = names(first.month)
  ) +
  theme_annual() +
  labs(title = "male") +
  theme(plot.title = element_text(size = 12, color = "blue"))

rose_female <- ggplot(data = count.f) +
  geom_hline(yintercept = scale.seq |> sqrt(), color = "grey80") +
  geom_vline(xintercept = c(0), color = "grey20") +
  geom_bar(aes(x = center * w, y = sqrt(count)),
    stat = "identity",
    fill = circular.phase.color(seq(from = 0, to = period, length.out = bins))
  ) +
  scale_y_continuous(
    name = "number of metabolites", limit = c(0, 9),
    breaks = scale.seq |> sqrt(), labels = scale.seq
  ) +
  coord_radial(expand = F, start = 0, end = 2 * pi, r.axis.inside = F, clip = "off") +
  scale_x_continuous(
    name = "phase (day)",
    limits = c(0, 2 * pi),
    breaks = first.month * w,
    labels = names(first.month)
  ) +
  theme_annual() +
  labs(title = "female") +
  theme(plot.title = element_text(size = 12, color = "red"))

rose_male + rose_female


pdf("Figures/sex_difference/rose_plot_SOMs_M.pdf", width = 3, height = 3)
print(rose_male)
dev.off()


pdf("Figures/sex_difference/rose_plot_SOMs_F.pdf", width = 3, height = 3)
print(rose_female)
dev.off()


pdf("Figures/sex_difference/rose_plot_SOMs_MF.pdf", width = 3, height = 6)
print(rose_male / rose_female)
dev.off()



r.test(smSOMs_amide$phase.male * w)
kuiper(smSOMs_amide$phase.male * w)

r.test(smSOMs_amide$phase.female * w)
kuiper(smSOMs_amide$phase.female * w)

rao.homogeneity(list(smSOMs_amide$phase.male * w, smSOMs_amide$phase.female * w))

watson.two(smSOMs_amide$phase.male * w, smSOMs_amide$phase.female * w)

phase.diff <- circular.diff(smSOMs_amide$phase.female, smSOMs_amide$phase.male)

r.test(phase.diff * w)

mu <- vm.ml(phase.diff * w, bias = T)
set.seed(777)
mu.ci <- vm.bootstrap.ci(phase.diff * w, bias = T)

# mean lag
mu$mu / w
mu.ci$mu.ci / w

plt_SD_phase <- ggplot(data = data.frame(phase.diff, row.names = smSOMs_amide$unique_id)) +
  geom_hline(yintercept = seq(10, 50, 10), color = "grey80") +
  geom_histogram(aes(x = phase.diff * w), bins = 12 * 3, fill = "palevioletred") +
  geom_vline(xintercept = c(0), color = "grey20") +
  annotate("pointrange",
    x = mu$mu, ,
    xmin = mu.ci$mu.ci[1],
    xmax = mu.ci$mu.ci[2],
    y = 55, color = "black", size = 0.05, linewidth = 0.75
  ) +
  scale_y_continuous(
    name = "number of metabolites", limit = c(0, 60),
    breaks = seq(10, 50, 10), transform = "sqrt"
  ) +
  coord_radial(expand = F, start = -1 * pi, end = pi, r.axis.inside = T, clip = "off") +
  scale_x_continuous(
    name = "phase difference (female - male, day)",
    limits = c(-1 * pi, pi),
    breaks = c(-1 * seq(0, period / 2 - 10, 30), seq(0, period / 2 - 10, 30), period / 2) * w,
    labels = c(-1 * seq(0, period / 2 - 10, 30), seq(0, period / 2 - 10, 30), period / 2)
  ) +
  theme_annual() +
  labs(title = "difference") +
  theme(plot.title = element_text(size = 12, color = "palevioletred"))

plt_SD_phase

write.csv(plt_SD_phase@data, "Figures/sex_difference/phase_sex-difference.csv")


pdf("Figures/sex_difference/phase_sex-difference.pdf", width = 3, height = 3)
plt_SD_phase
dev.off()

pdf("Figures/sex_difference/rose_plot_SOMs_MFdiff.pdf", width = 9, height = 3)
print(rose_male + rose_female + plt_SD_phase)
dev.off()


plt_SD_amplitude <- ggplot(data = smSOMs_amide) +
  geom_point(aes(
    x = amplitude.male,
    y = amplitude.female,
    color = (amplitude.female - amplitude.male)
  ), size = 2) +
  scale_color_gradient2(
    name = "difference",
    low = "#3970eb", mid = "#AAAAAA", high = "#eb7093"
  ) +
  scale_x_continuous(name = "male amplitude", limits = c(0, 5.5), expand = 0) +
  scale_y_continuous(name = "female amplitude", limits = c(0, 5.5), expand = 0) +
  geom_abline(
    slope = 1 |> rep(3),
    intercept = c(-log2(1.5), 0, log2(1.5)),
    color = "darkblue", linetype = c(2, 1, 2)
  ) +
  theme_annual() +
  theme(legend.position = c(0.9, 0.3), legend.key.size = unit(0.4, "cm"))

plt_SD_amplitude +
  geom_text_repel(
    aes(x = amplitude.male, y = amplitude.female, label = curated_name),
    data = smSOM.representatives[3:4, ],
    size = 3,
    box.padding = 5,
    point.padding = 0,
    min.segment.length = 0,
    force = 10,
    max.overlaps = Inf
  )

pdf("Figures/sex_difference/amplitude_sex-difference.pdf", width = 3, height = 3)
plt_SD_amplitude
dev.off()


pdf("Figures/sex_difference/amplitude_sex-difference_annotation.pdf", width = 3, height = 3)
plt_SD_amplitude +
  geom_text_repel(
    aes(x = amplitude.male, y = amplitude.female, label = curated_name),
    data = smSOM.representatives[3:4, ],
    size = 3,
    box.padding = 0.8,
    point.padding = 0,
    min.segment.length = 0,
    force = 0.5,
    max.overlaps = Inf
  )
dev.off()



ids <- smSOMs_amide$unique_id

features <- human_amide_BGI[, sample_metadata$Sample]
rownames(features) <- human_amide_BGI$unique_id
##### Heatmap drawing #####
sd.cutoff <- 3

mat <- matrix(nrow = length(ids))
metadata <- sample_metadata[sample_metadata$Sex == "male", ]
mat <- matrix(0, ncol = 12, nrow = length(ids))
for (i in unique(metadata$Subject)) {
  tmp <- (features[ids, grep(i, colnames(features))])
  tmp <- log2(tmp)
  tmp <- tmp - (apply(tmp, 1, mean))
  mat <- mat + tmp
}
mat.m <- mat / 6


mat <- matrix(nrow = length(ids))
metadata <- sample_metadata[sample_metadata$Sex == "female", ]
mat <- matrix(0, ncol = 12, nrow = length(ids))
for (i in unique(metadata$Subject)) {
  tmp <- (features[ids, grep(i, colnames(features))])
  tmp <- log2(tmp)
  tmp <- tmp - (apply(tmp, 1, mean))
  mat <- mat + tmp
}
mat.f <- mat / 5

mat <- cbind(mat.m, mat.f)

mat <- mat / apply(mat, 1, sd, na.rm = T)

mat[abs(mat) > sd.cutoff & !is.na(mat)] <-
  sd.cutoff * sign(mat[abs(mat) > sd.cutoff & !is.na(mat)])

mat.m <- mat[smSOMs_amide[order(smSOMs_amide$phase.male), ]$unique_id, 1:12]
mat.f <- mat[smSOMs_amide[order(smSOMs_amide$phase.female), ]$unique_id, 13:24]

df.row <- data.frame(
  "common" = ifelse(smSOMs_amide$unique_id %in% SOMs_amide$unique_id, "common", "specific"),
  "Ion mode" = smSOMs_amide$Ion.Mode,
  "RT" = smSOMs_amide$RT,
  "m/z" = smSOMs_amide$CalcMz,
  row.names = smSOMs_amide$unique_id
)

colnames(mat.m) <- colnames(mat.f) <- month.abb
colors <- colorRampPalette(
  colors = c("#0050FF", "#FFF", "#FF5000")
)(300 + 1)


pdf("figures/sex_difference/heatmap_sdSOM_male.pdf", width = 4, height = 3)
pheatmap::pheatmap(mat.m,
  col = colors,
  annotation_row = df.row[rownames(mat.m), ],
  annotation_colors = list(common = c(common = "grey20", specific = "grey90")),
  cluster_rows = F,
  cluster_cols = F,
  legend = T,
  angle_col = 90,
  scale = "none",
  na_col = "#FFFFFFFF",
  labels_row = "",
  show_rownames = F,
  show_colnames = T,
  legend_breaks = c(-3, 0, 3),
  legend_labels = c("< -3", "0", "> 3"),
  fontsize = 8
)

dev.off()

pdf("figures/sex_difference/heatmap_sdSOM_female.pdf", width = 4, height = 3)
pheatmap::pheatmap(mat.f,
  col = colors,
  annotation_row = df.row[rownames(mat.f), ],
  annotation_colors = list(common = c(common = "grey20", specific = "grey90")),
  cluster_rows = F,
  cluster_cols = F,
  legend = T,
  angle_col = 90,
  scale = "none",
  na_col = "#FFFFFFFF",
  labels_row = "",
  show_rownames = F,
  show_colnames = T,
  legend_breaks = c(-3, 0, 3),
  legend_labels = c("< -3", "0", "> 3"),
  fontsize = 8
)

dev.off()




library(eulerr)

pdf("figures/sex_difference/sdSOM-SOM_venn.pdf", width = 5, height = 4)

euler(c(
  additive = 66,
  interaction = 58,
  "additive&interaction" = 250
)) |> plot(
  quantities = TRUE,
  fills = list(fill = c("#EEEEEE", "palevioletred"), alpha = 0.25),
  labels = list(col = "black")
)

dev.off()
