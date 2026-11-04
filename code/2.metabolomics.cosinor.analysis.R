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

table(interaction(human_amide_BGI$Ion.Mode, human_amide_BGI$Level))


################################################################################
#   metabolome data fitting  #
################################################################################


features <- human_amide_BGI[, sample_metadata$Sample]
rownames(features) <- human_amide_BGI$unique_id

## outlier removal
norm.data <- apply(features, 1, function(x) log2(x) |> scale()) |> t()
sum(abs(norm.data) > 3)
# 10127
sum(abs(norm.data) > 3) / prod(dim(norm.data)) * 100 # %
# 0.899199 % were set as outliers

features[abs(norm.data) > 3] <- NA

# fitting
lmm_cosinor.result <- data.frame(row.names = human_amide_BGI$unique_id)

for (m in rownames(features)) {
  y <- features[m, ] |>
    log2() |>
    unlist()
  suppressMessages({
    fit.cosinor <- lmer(y ~ sin.a + cos.a + Sex + (1 | Subject),
      data = sample_metadata, REML = F
    )
  })
  para <- fixef(fit.cosinor)
  para["p.value"] <- contest(fit.cosinor, diag(4)[c(2:3), ])[1, "Pr(>F)"]
  lmm_cosinor.result[m, names(para)] <- para
}

## FDR correction

lmm_cosinor.result$q.value <- qvalue(lmm_cosinor.result$p.value)$qvalues

## parameter calculation
lmm_cosinor.result$amplitude <- with(lmm_cosinor.result, amplitude(cos.a, sin.a) * 2)
lmm_cosinor.result$phase <- with(lmm_cosinor.result, phase(cos.a, sin.a))

lmm_cosinor.result$mesor.male <- with(lmm_cosinor.result, `(Intercept)`)
lmm_cosinor.result$mesor.female <- with(lmm_cosinor.result, `(Intercept)` + Sexfemale)

# sort by FDR
lmm_cosinor.result <- lmm_cosinor.result[order(lmm_cosinor.result$q.value), ]

# export parameter table
merge(human_amide_BGI[, 1:13], lmm_cosinor.result, by.x = "unique_id", by.y = 0) |>
  write.csv("results/human_metabolome_amide_mixed-effects-cosinor_annual.csv", na = "", row.names = F)

saveRDS(
  merge(human_amide_BGI[, 1:13], lmm_cosinor.result, by.x = "unique_id", by.y = 0),
  "results/human_amide_lmm_cosinor_annual.rds"
)

sum(lmm_cosinor.result$q.value < 0.05)


seasonal.metabolites <- readRDS("results/human_Amide_lmm_cosinor_annual.rds")
row.names(seasonal.metabolites) <- seasonal.metabolites$unique_id


q.cutoff <- 0.05

SOMs_amide <- subset(seasonal.metabolites, q.value < q.cutoff)
SOMs_amide <- SOMs_amide[order(SOMs_amide$phase), ]

range(SOMs_amide$amplitude)

2^range(SOMs_amide$amplitude)

features <- human_amide_BGI[, sample_metadata$Sample]
rownames(features) <- human_amide_BGI$unique_id

ids <- SOMs_amide$unique_id

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

mat.m <- mat.m / apply(mat.m, 1, sd)
mat.f <- mat.f / apply(mat.f, 1, sd)

mat.m[abs(mat.m) > sd.cutoff] <- sd.cutoff * sign(mat.m[abs(mat.m) > sd.cutoff])
mat.f[abs(mat.f) > sd.cutoff] <- sd.cutoff * sign(mat.f[abs(mat.f) > sd.cutoff])

colnames(mat.m) <- colnames(mat.f) <- month.abb
colors <- colorRampPalette(
  colors = c("#0050FF", "#FFF", "#FF5000")
)(300 + 1)


pdf("figures/SOMs/heatmap_SOMs.pdf", width = 4, height = 3)
pheatmap::pheatmap(cbind(mat.m, " " = NA, mat.f),
  col = colors,
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

write.csv(cbind(mat.m, mat.f), "figures/SOMs/heatmap_SOMs.csv")



##### rose plot for peak phase distributions #####

bins <- 25

# unimodality : T
r.test(SOMs_amide$phase * w)
# non-uniformity : T
kuiper(SOMs_amide$phase * w)
# non-uniformity : T
watson(SOMs_amide$phase * w)
# non-uniformity : T
rao.spacing(SOMs_amide$phase * w)

(vm.ml(SOMs_amide$phase * w)$mu / w) %% period


count <- circular.SW.count(SOMs_amide$phase, bins = bins, period = period)
scale.seq <- seq(10, 70, 20)

rose_phase <- ggplot(data = count) +
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
  scale_x_continuous(
    name = "phase (day)",
    limits = c(0, 2 * pi),
    breaks = sort(month.label %% period) * w,
    labels = names(sort(month.label %% period))
  ) +
  coord_radial(expand = F, start = 0, end = 2 * pi, r.axis.inside = F, clip = "off") +
  theme_classic()


rose_phase

write.csv(count, "figures/SOMs/phase_rose_plot.csv")


pdf("Figures/SOMs/phase_rose_plot.pdf", width = 4, height = 4)
rose_phase
dev.off()


library(ggrepel)

#############

volcano_plt <- ggplot(data = seasonal.metabolites) +
  geom_point(aes(x = amplitude, y = -log10(q.value)), color = "grey90") +
  geom_abline(
    slope = 0, intercept = -log10(0.05),
    color = "black", linetype = "dotted"
  ) +
  geom_point(
    aes(x = amplitude, y = -log10(q.value), color = factor(phase)),
    data = SOMs_amide[order(SOMs_amide$unique_id %in% SOM.representatives$unique_id), ],
    size = 2
  ) +
  scale_color_manual(values = circular.phase.color(SOMs_amide$phase)) +
  scale_x_continuous(name = "log2 amplitude (peak-to-trough)", limits = c(0, 5), expand = c(0, 0)) +
  scale_y_continuous(name = "-log10( q-value ) ", limits = c(0, 40), expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = "n")


pdf("Figures/SOMs/volcano_annotated.pdf", width = 4, height = 4)

volcano_plt +
  geom_text_repel(
    aes(x = amplitude, y = -log10(q.value), label = curated_name),
    data = subset(SOM.representatives, q.value < 1e-3 | amplitude > 0.5),
    size = 3,
    box.padding = 1.1,
    point.padding = 0,
    min.segment.length = 0,
    force = 100,
    max.overlaps = Inf
  )

dev.off()


features <- human_amide_BGI[, sample_metadata$Sample]
rownames(features) <- human_amide_BGI$unique_id

## outlier removal
norm.data <- apply(features, 1, function(x) log2(x) |> scale()) |> t()
sum(abs(norm.data) > 3)
# 10127
features[abs(norm.data) > 3] <- NA

SOM.representatives <- SOMs_amide[
  c(
    "Amide_POS_2.496_150.05861", # methionine
    "Amide_POS_13.705_241.03072", # cystine
    "Amide_NEG_13.553_296.03808", # cys-ss-cys-gly
    "Amide_NEG_13.473_353.05969", # putative gly-cys-ss-cys-gly
    "Amide_POS_8.431_173.09247", # gly-pro
    "Amide_POS_1.317_121.03981", # urocanate
    # "Amide_NEG_14.026_159.11394", #N(6)−Methyllysine

    "Amide_NEG_9.818_135.03001", # threonic acid
    "Amide_NEG_9.207_72.99319", # glyoxylate,

    "Amide_NEG_13.989_427.00752", # IDP
    "Amide_POS_3.732_312.13013", # dimethylguanosine
    "Amide_NEG_13.743_331.04546", # dIMP
    "Amide_POS_3.947_298.11462", # methylguanosine

    "Amide_NEG_6.955_800.65404", # PC(20:0/P−18:0)
    "Amide_NEG_6.954_716.69172", # Cer(d18:2/29:0)
    "Amide_NEG_1.076_553.48298", # DG(14:0/17:0/0:0)[iso2]
    "Amide_POS_0.745_579.4988", # DG(16:1(9Z)/17:1(9Z)/0:0)[iso2]

    "Amide_POS_0.869_583.25515", # biliverdin
    "Amide_NEG_0.861_299.2015" # retinoic acid AtRA
  ),
]

SOM.representatives$curated_name <-
  c(
    "methionine",
    "cystine",
    "cys-ss-cys-gly",
    "gly-cys-ss-cys-gly",
    "gly-pro",
    "urocanic acid",
    "threonic acid",
    "glyoxylic acid",
    "IDP",
    "N2,N2-dimethylguanosine",
    "dIMP",
    "N2-methylguanosine",
    "PC(38:0)",
    "Cer(47:2)",
    "DG(31:0)",
    "DG(33:2)",
    "biliverdin",
    "all-trans retinoic acid"
  )


representative.plots <- plot_metabolite(SOM.representatives, centered = T)

representative.plots[[17]] <-
  representative.plots[[17]] + theme(axis.text.x = element_text(angle = 45, hjust = 1))

representative.plots[[18]] <-
  representative.plots[[18]] + theme(axis.text.x = element_text(angle = 45, hjust = 1))

lapply(representative.plots, function(x) {
  (x + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
})


pdf("figures/SOMs/representative_SOMs_plots.pdf", width = 6, height = 16)
representative.plots |>
  wrap_plots(ncol = 2, byrow = T) |>
  print()
dev.off()

merge(SOM.representatives, features, by.x = "unique_id", by.y = 0) |>
  write.csv("figures/SOMs/representative_SOMs_plots.csv", na = "", row.names = F)
