library(pls)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(tidyverse)


# load data
load("data/Human_annual_metabolomics.RData")

source("code/supporting.functions.R")

seasonal.metabolites <- readRDS("results/human_amide_lmm_cosinor_annual.rds")
rownames(seasonal.metabolites) <- seasonal.metabolites$unique_id

SOMs_amide <- subset(seasonal.metabolites, q.value < 0.05)

ids <- SOMs_amide$unique_id

plsr.pred <- data.frame(sample_metadata, row.names = sample_metadata$Sample)

plsr.models <- list()

ncomp <- 10
exprained.variance.threshold <- 0.01
accepted.error <- 30

RS.phase <- function(cos.p, sin.p, cos.a, sin.a) {
  error <- (circular.diff(phase(cos.a, sin.a), phase(cos.p, sin.p)))
  RS <- rsq(phase(cos.a, sin.a) + error, phase(cos.a, sin.a))
  return(RS)
}


################################################################################

for (subject.out in LETTERS[1:11]) {
  ## scaling data
  sex.train <- ifelse(grepl("[ABCDEF]", subject.out), "male", "female")

  train.metadata <- subset(sample_metadata, Subject != subject.out)
  test.metadata <- subset(sample_metadata, Subject == subject.out)

  train.features <- human_amide_BGI[, train.metadata$Sample] |> log2()
  test.features <- human_amide_BGI[, test.metadata$Sample] |> log2()

  rownames(train.features) <- rownames(test.features) <- human_amide_BGI$unique_id

  train.features <- train.features[ids, ]
  test.features <- test.features[ids, ]

  samples.m <- subset(train.metadata, Sex == "male")$Sample
  samples.f <- subset(train.metadata, Sex == "female")$Sample

  center.m <- apply(train.features[, samples.m], 1, mean, na.rm = T)
  center.f <- apply(train.features[, samples.f], 1, mean, na.rm = T)

  stdev.m <- apply(train.features[, subset(train.metadata, Sex == "male")$Sample], 1, sd, na.rm = T)
  stdev.f <- apply(train.features[, subset(train.metadata, Sex == "female")$Sample], 1, sd, na.rm = T)

  train.features[, subset(train.metadata, Sex == "male")$Sample] <-
    apply(train.features[, subset(train.metadata, Sex == "male")$Sample], 2, function(x) (x - center.m) / stdev.m)

  train.features[, subset(train.metadata, Sex == "female")$Sample] <-
    apply(train.features[, subset(train.metadata, Sex == "female")$Sample], 2, function(x) (x - center.f) / stdev.f)

  ## set inner LOOCV sample lists
  subject.set <- LETTERS[1:11]
  innerCV_list <-
    lapply(
      setdiff(subject.set, subject.out),
      function(x) {
        grep(x, train.metadata$Sample)
      }
    )

  ## dependent variables
  cos.a <- cos(train.metadata$day_of_year * w)
  sin.a <- sin(train.metadata$day_of_year * w)

  ## independent variables (predictor variables)
  p <- train.features[ids, ] |> t()

  ## perform PLSR
  pls.c <- plsr(cos.a ~ p, ncomp, validation = "CV", segments = innerCV_list, center = F, scale = F)
  pls.s <- plsr(sin.a ~ p, ncomp, validation = "CV", segments = innerCV_list, center = F, scale = F)

  ## extract LOOCV prediction
  CV.cos.pred <- as.list(as.data.frame(pls.c$validation$pred))
  CV.sin.pred <- as.list(as.data.frame(pls.s$validation$pred))

  ### optimal max components selection
  c.ncomp <- 1 + which(
    diff(unlist(lapply(CV.cos.pred, rsq, y = cos.a))) > exprained.variance.threshold
  ) |>
    max()
  s.ncomp <- 1 + which(
    diff(unlist(lapply(CV.sin.pred, rsq, y = sin.a))) > exprained.variance.threshold
  ) |>
    max()

  ###

  #### optimal LVs components searching
  ## evaluate optimal combination of the number of latent variables
  CV.error <- expand.grid("cos" = 1:c.ncomp, "sin" = 1:s.ncomp)
  CV.error$total <- CV.error$cos + CV.error$sin
  for (c in 1:c.ncomp) {
    for (s in 1:s.ncomp) {
      CV.error[CV.error$cos == c & CV.error$sin == s, "value"] <-
        RS.phase(
          CV.cos.pred[[c]], CV.sin.pred[[s]],
          cos.a, sin.a
        )
    }
  }

  opt <- CV.error[which.max(CV.error$value), 1:2] |>
    unlist()

  pls.c$candidate.ncomp <- c.ncomp
  pls.s$candidate.ncomp <- s.ncomp

  ## set optimal #LVs
  names(opt) <- c("cos", "sin")

  ### optimal LVs components determined
  pls.c$opt.ncomp <- opt.c <- opt[1]
  pls.s$opt.ncomp <- opt.s <- opt[2]

  ## extract selected PLS coefficients
  pls.c$opt.coefficients <- B.c <- pls.c$coefficients[, , opt.c]
  pls.s$opt.coefficients <- B.s <- pls.s$coefficients[, , opt.s]

  ### evaluate model by data of a left subject
  ## scaled by train data mean and standard deviation
  if (sex.train == "male") {
    q <- apply(test.features, 2, function(x) (x - center.m) / stdev.m) |> t()
  } else if (sex.train == "female") {
    q <- apply(test.features, 2, function(x) (x - center.f) / stdev.f) |> t()
  }

  ## input data
  q <- q[, ids]
  error.test <- c()
  ## estimated cos and sin values, phase, amplitude and actual sampling day
  for (sample in test.metadata$Sample) {
    c.coord <- sum(B.c * q[sample, ])
    s.coord <- sum(B.s * q[sample, ])

    plsr.pred[sample, "cos.coordinate"] <- c.coord
    plsr.pred[sample, "sin.coordinate"] <- s.coord

    plsr.pred[sample, "amplitude"] <- amplitude(c.coord, s.coord)
    plsr.pred[sample, "predicted"] <- phase(c.coord, s.coord)
    plsr.pred[sample, "actual"] <- subset(test.metadata, Sample == sample)$day_of_year
  }

  test.error <- plsr.pred[plsr.pred$Subject == subject.out, c("cos.coordinate", "sin.coordinate", "cos.a", "sin.a")]
  test.metric <- RS.phase(test.error[[1]], test.error[[2]], test.error[[3]], test.error[[4]])
  print(c(
    paste("model", subject.out),
    paste("optimal cos #components=", opt.c),
    paste("optimal sin #components=", opt.s),
    paste(
      round(sum(CV.error$value > test.metric) / nrow(CV.error) * 100, 1),
      "% CV combination better than test (potency of overfitting)"
    )
  ))

  ## save models and related parameters
  plsr.models[[subject.out]] <-
    list(
      "s.model" = pls.s, "c.model" = pls.c,
      "model.parameters" = list(
        train.metadata = train.metadata,
        scaling = list(
          male = list(center = center.m, scale = stdev.m),
          female = list(center = center.f, scale = stdev.f)
        ),
        RS.phase = RS.phase,
        exprained.variance.threshold = exprained.variance.threshold,
        LOOCV = CV.error
      ),
      "test.evaluation" = list(
        test.input = q,
        test.metadata = test.metadata,
        test.error = test.error,
        test.metric = test.metric
      )
    )
}

## errors between predicted and actual day
plsr.pred$error <- circular.diff(
  plsr.pred$predicted,
  plsr.pred$actual
)


saveRDS(plsr.models, "results/PLSR_calendar_models.rds")

write.csv(plsr.pred, "results/PLSR.LOO-prediction.csv", row.names = F)

lapply(plsr.models, function(model) {
  model$c.model$opt.ncomp
}) |>
  unlist() |>
  mean()
lapply(plsr.models, function(model) {
  model$s.model$opt.ncomp
}) |>
  unlist() |>
  mean()


## some indices for model performance
MdAE <- median(abs(plsr.pred$error))
IQR <- IQR(abs(plsr.pred$error))

RS.plsr <- rsq(plsr.pred$actual, plsr.pred$actual + plsr.pred$error)

accepted.ratio <- sum(abs(plsr.pred$error) < accepted.error) / nrow(plsr.pred) |> round(3)

print(paste(
  "RS=", RS.plsr |> round(3), "; MdAE ± IQR =", MdAE |> round(2), "±", IQR |> round(2), ";",
  round(accepted.ratio, 3) * 100, "% samples <", accepted.error, "day error"
))


## correlation between actual and estimated
plt_PLSR_RS <- ggplot(data = subset(plsr.pred)) +
  geom_point(aes(
    x = actual,
    y = predicted,
    color = Subject,
    fill = Subject,
    shape = Subject
  ), size = 2) +
  scale_Subject() +
  scale_x_continuous(name = "actual sampling day", limits = c(0, period), expand = c(0.01, 0.01)) +
  scale_y_continuous(name = "predicted sampling day", limits = c(0, period), expand = c(0.01, 0.01)) +
  annotate("text", x = period / 2, y = 360, label = paste("RS =", round(RS.plsr, 3)), size = 4) +
  geom_abline(intercept = 0, slope = 1, color = "grey20") +
  theme_annual()
plt_PLSR_RS

## check error distribution across the sampling month
plt_PLSR_monthly_error <-
  ggplot(data = subset(plsr.pred, Sex == "male")) +
  geom_point(aes(
    x = Month,
    y = error,
    shape = Subject,
    fill = Subject,
    color = Subject
  )) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20") +
  geom_line(aes(
    x = Month,
    y = error,
    color = Subject
  )) +
  scale_x_continuous(breaks = 0:12, name = "month") +
  scale_y_continuous(name = "error (estimated - actual day)", limits = c(-period / 2, period / 2)) +
  scale_Subject(1:6) +
  theme_annual() +
  theme(legend.position = "top") +
  ggplot(data = subset(plsr.pred, Sex == "female")) +
  geom_point(aes(
    x = Month,
    y = error,
    shape = Subject,
    fill = Subject,
    color = Subject
  )) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20") +
  geom_line(aes(
    x = Month,
    y = error,
    color = Subject
  )) +
  scale_x_continuous(breaks = 0:12, name = "month") +
  scale_y_continuous(name = "error (estimated - actual day)", limits = c(-period / 2, period / 2)) +
  scale_Subject(7:11) +
  theme_annual() +
  theme(legend.position = "top")

plt_PLSR_monthly_error


### predicted amplitude weakly correlated with prediction error
# probably including intrinsic body season variation

ggplot(data = plsr.pred) +
  geom_point(aes(x = amplitude, y = abs(error), color = Subject)) +
  scale_color_manual(values = subject.color(1:11)) +
  geom_hline(yintercept = 0, color = "grey20") +
  geom_hline(yintercept = MdAE, color = "grey20", linetype = "dotted") +
  geom_vline(xintercept = 1, color = "grey20", linetype = "dotted") +
  scale_y_continuous(
    limits = c(-10, period / 2), breaks = c(seq(0, 180, 30), MdAE),
    labels = c(seq(0, 180, 30), MdAE |> round(2)), name = "absolute error (day)"
  ) +
  theme_classic()

lmerTest::lmer(abs(error) ~ amplitude + (1 | Subject), data = plsr.pred) |>
  anova()


plt_PLSR_monthly_error_withMeanSD <-
  ggplot(data = plsr.pred) +
  stat_summary(aes(x = Month - 0.2, y = error),
    fun.data = "mean_sdl",
    size = 1, color = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20") +
  scale_x_continuous(breaks = 0:12, name = "month") +
  scale_y_continuous(name = "error (estimated - actual day)", limits = c(-period / 2, period / 2)) +
  # scale_size_continuous(limits = c(0,2.5))+
  theme_classic() +
  theme(legend.position = "top") +
  ggplot(data = subset(plsr.pred, Sex == "male")) +
  stat_summary(aes(x = Month - 0.2, y = error),
    fun.data = "mean_sdl",
    size = 1, color = "blue"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20") +
  scale_x_continuous(breaks = 0:12, name = "month") +
  scale_y_continuous(name = "error (estimated - actual day)", limits = c(-period / 2, period / 2)) +
  # scale_size_continuous(limits = c(0,2.5))+
  theme_classic() +
  theme(legend.position = "top") +
  ggplot(data = subset(plsr.pred, Sex == "female")) +
  stat_summary(aes(x = Month - 0.2, y = error),
    fun.data = "mean_sdl",
    size = 1, color = "red"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20") +
  scale_x_continuous(breaks = 0:12, name = "month") +
  scale_y_continuous(name = "error (estimated - actual day)", limits = c(-period / 2, period / 2)) +
  # scale_size_continuous(limits = c(0,2.5))+
  theme_classic() +
  theme(legend.position = "top")

plt_PLSR_monthly_error_withMeanSD


plt_PLSR_subject_error <-
  ggplot(data = plsr.pred) +
  geom_boxplot(aes(x = Subject, y = abs(error), fill = Subject), color = "black") +
  scale_fill_manual(values = subject.color()) +
  geom_hline(yintercept = 0, color = "grey20") +
  annotate("text", y = period / 2, x = 6, label = paste("MdAE + IQR =", round(MdAE, 2), "+", round(IQR, 2)), size = 4) +
  scale_y_continuous(
    limits = c(-10, period / 2), breaks = c(seq(0, 180, 30), MdAE),
    labels = c(seq(0, 180, 30), MdAE |> round(2)), name = "absolute error (day)"
  ) +
  theme_classic() +
  theme(legend.position = "none")

plt_PLSR_subject_error

## draw error cummulative distribution curve
x <- sort(abs(plsr.pred$error))

AUC.plsr <- AUC(x)


plt_PLSR_ECD <- plot_AUC(list("PLSR" = AUC.plsr), colors = "black")
plt_PLSR_ECD

pdf("figures/PLSR/PLSR_cummulative_error_distribution.pdf", width = 3.5, height = 3.5)
plt_PLSR_ECD
dev.off()

write.csv(AUC.plsr[1:2] |> as.data.frame(), "figures/PLSR/PLSR_cummulative_error_distribution.csv")

## calculate significance of PLS coefficents
# significance was defined as Z-statistics
# mean and sd of coefficients were calculated from 11 LOO PLSR models
# indicating contribution of each features and its robustness to estimation

plsr.coef.cos <- lapply(plsr.models, function(model) {
  model$c.model$opt.coefficients
}) |> as.data.frame()

plsr.coef.sin <- lapply(plsr.models, function(model) {
  model$s.model$opt.coefficients
}) |> as.data.frame()

mean.beta.cos <- plsr.coef.cos |>
  apply(1, function(x) mean(x))

mean.beta.sin <- plsr.coef.sin |>
  apply(1, function(x) mean(x))

## absolute mean of coefficients (for ranking metabolites)
sum.effect.size <- abs(mean.beta.sin) + abs(mean.beta.cos)

### scatter plot for significances
plt_significance <- ggplot(data = data.frame(
  mean.beta.cos, mean.beta.sin, sum.effect.size
)) +
  geom_point(aes(
    x = mean.beta.cos |> abs(),
    y = mean.beta.sin |> abs(),
    color = sum.effect.size
  ), size = 3) +
  scale_x_continuous(
    name = "mean absolute cos coefficients",
    limits = c(0, 0.07), expand = 0
  ) +
  scale_y_continuous(
    name = "mean absolute sin coefficients",
    limits = c(0, 0.07), expand = 0
  ) +
  scale_color_gradient(name = "sum", low = "grey70", high = "darkorange") +
  theme_classic() +
  theme(legend.position = c(0.9, 0.5))


plt_significance


pdf("figures/PLSR/coefficients_cos-sin_scatter.pdf", width = 3.5, height = 3.5)
plt_significance
dev.off()

plt_significance@data |>
  write.csv("figures/PLSR/coefficients_cos-sin_scatter.csv")


### heatmap of coefficients ###
pdf("figures/PLSR/mean-coefficients_heatmap.pdf", width = 2, height = 3)
pheatmap::pheatmap(
  as.matrix(data.frame(
    "β cos" = mean.beta.cos,
    "β sin" = mean.beta.sin
  )),
  cluster_rows = F, cluster_cols = F, na_col = "white", gaps_col = 1,
  color = colorRampPalette(c("darkgreen", "white", "darkmagenta"))(100 + 1),
  show_rownames = F, legend_breaks = seq(-10, 10, 2),
  legend = T, angle_col = 90, annotation_names_row = F,
  scale = "column"
)

dev.off()




### heatmap of coefficients ###
pdf("figures/PLSR/all-coefficients_heatmap.pdf", width = 6, height = 3)
pheatmap::pheatmap(
  as.matrix(cbind(as.data.frame(plsr.coef.cos),
    " " = NA,
    as.data.frame(plsr.coef.sin)
  )),
  cluster_rows = F, cluster_cols = F, na_col = "white",
  color = colorRampPalette(c("darkgreen", "white", "darkmagenta"))(100 + 1),
  show_rownames = F, legend_breaks = seq(-0.05, 0.05, 0.025),
  legend = T, angle_col = 90, annotation_names_row = F,
  scale = "none"
)

dev.off()


## save figures
pdf("figures/PLSR/PLSR_prediction_subject_error.pdf", width = 5, height = 4)
plt_PLSR_subject_error
dev.off()


pdf("figures/PLSR/PLSR_prediction_monthly_error.pdf", width = 8, height = 4)
plt_PLSR_monthly_error
dev.off()


pdf("figures/PLSR/PLSR_prediction-vs-actual.pdf", width = 4.5, height = 5)
plt_PLSR_RS
dev.off()


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

SOMs_amide.out <- merge(SOMs_amide, cbind(
  "PLSR.cos" = mean.beta.cos,
  "PLSR.sin" = mean.beta.sin,
  "PLSR.sum.effectsize" = sum.effect.size
),
by.x = "unique_id", by.y = 0
) |>
  merge(SOM.representatives[, c("unique_id", "curated_name")], all.x = T)


SOMs_amide.out <-
  SOMs_amide.out[, c(
    "unique_id", "Name", "curated_name", "Level", "Match", "sin.a", "cos.a",
    "p.value", "q.value", "amplitude", "phase",
    "mesor.male", "mesor.female",
    "PLSR.cos", "PLSR.sin", "PLSR.sum.effectsize"
  )]

colnames(SOMs_amide.out) <- c(
  "unique_id", "Name", "Curated_Name", "Level", "Match", "sin", "cos",
  "p.value", "q.value", "amplitude", "phase",
  "mesor.male", "mesor.female",
  "PLSR.cos", "PLSR.sin", "PLSR.sum.effectsize"
)

write.csv(SOMs_amide.out, "results/SOM.table.csv", row.names = F, na = "")
