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

LOO_lmm_cosinor <- readRDS("results/LOO_human_amide_lmm_cosinor_annual.rds")
LOO_lmm_cosinor <- lapply(LOO_lmm_cosinor, function(df) {
  df$unique_id <- row.names(df)
  df <- merge(human_amide_BGI[, 1:20], df, by = "unique_id")
  df
})

LOO_lmm_cosinor |>
  lapply(function(df) {
    subset(df, q.value < 0.05) |> nrow()
  }) |>
  unlist() |>
  min()


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
  SOMs <- subset(LOO_lmm_cosinor[[subject.out]], q.value < 0.05)
  ids <- SOMs$unique_id
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
  q <- q[test.metadata$Sample, ids]
  error.test <- c()
  ## estimated cos and sin values, phase, amplitude and actual sampling day
  c.coord <- q %*% B.c
  s.coord <- q %*% B.s

  plsr.pred[rownames(q), "cos.coordinate"] <- c.coord
  plsr.pred[rownames(q), "sin.coordinate"] <- s.coord

  plsr.pred[rownames(q), "amplitude"] <- amplitude(c.coord, s.coord)
  plsr.pred[rownames(q), "predicted"] <- phase(c.coord, s.coord)
  plsr.pred[rownames(q), "actual"] <- test.metadata$day_of_year


  test.error <- plsr.pred[plsr.pred$Subject == subject.out, c("cos.coordinate", "sin.coordinate", "cos.a", "sin.a")]
  test.metric <- RS.phase(test.error[[1]], test.error[[2]], test.error[[3]], test.error[[4]]) |>
    as.numeric()
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
  scale_x_continuous(
    name = "actual sampling day", limits = c(0, period),
    breaks = seq(0, 365, 30), expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    name = "predicted sampling day", limits = c(0, period),
    breaks = seq(0, 365, 30), expand = c(0.01, 0.01)
  ) +
  annotate("text",
    x = period / 2, y = 360,
    label = paste("RS =", round(RS.plsr, 3)), size = 4
  ) +
  geom_abline(intercept = 0, slope = 1, color = "grey20") +
  theme_annual() +
  coord_fixed()
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
plt_PLSR_ECD + coord_fixed()

pdf("figures/PLSR/PLSR_cummulative_error_distribution2.pdf", width = 3.5, height = 3.5)
plt_PLSR_ECD + coord_fixed()
dev.off()

write.csv(AUC.plsr[1:2] |> as.data.frame(), "figures/PLSR/PLSR_cummulative_error_distribution.csv")



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


## testing resolution of PLSR prediction across months

bs <- rbind(
  subset(plsr.pred, Month == 12),
  subset(plsr.pred, Month == 1),
  subset(plsr.pred, Month == 2)
)
bs$actual <- ifelse(bs$Month == 12, bs$actual - 365, bs$actual + 365)
bs$Month <- ifelse(bs$Month == 12, bs$Month - 12, bs$Month + 12)

bs <- rbind(plsr.pred, bs)
bs <- bs[order(bs$Subject), ]

lapply(0:11, function(n) {
  t.test(
    with(subset(bs, Month == n), actual + error),
    with(subset(bs, Month == (n + 1)), actual + error),
    paired = T
  )[["p.value"]]
}) |>
  p.adjust("holm") |>
  signif(3)


lapply(0:11, function(n) {
  t.test(
    with(subset(bs, Month == n), actual + error),
    with(subset(bs, Month == (n + 2)), actual + error),
    paired = T
  )[["p.value"]]
}) |>
  p.adjust("fdr") |>
  signif(2)
