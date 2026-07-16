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

plsr.models <- readRDS("results/PLSR_calendar_models.rds")


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


################################################################################


ids <- sort(sum.effect.size, decreasing = T) |> names()

topN <- c(seq(10, 50, 5), seq(60, 100, 10), seq(120, 200, 20), 250, 300)


set.name <- "plsr_selected"

error <- list()

for (n in topN) {
  timetable.pred <- data.frame()
  SOMs <- SOMs_amide[ids[1:n], ]

  for (subject.out in LETTERS[1:11]) {
    sex.train <- ifelse(grepl("[ABCDEF]", subject.out), "male", "female")
    test.metadata <- subset(sample_metadata, Subject == subject.out) |> as.data.frame()
    train.metadata <- subset(sample_metadata, Subject != subject.out & Sex == sex.train)
    train.features <- human_amide_BGI[, train.metadata$Sample] |> log2()
    test.features <- human_amide_BGI[, test.metadata$Sample] |> log2()
    rownames(train.features) <- rownames(test.features) <- human_amide_BGI$unique_id
    center <- apply(train.features, 1, mean, na.rm = T)
    stdev <- apply(train.features, 1, sd, na.rm = T)
    test.features[, test.metadata$Sample] <-
      apply(test.features, 2, function(x) (x - center) / stdev)
    for (month in test.metadata$Sample) {
      y <- test.features[SOMs$unique_id, month] |> unlist()
      ref.phase <- SOMs[, "phase"]

      fit.lm <- glm(y ~ sin(ref.phase * w) + cos(ref.phase * w))
      para <- coef(fit.lm)
      para["RS"] <- 1 - deviance(fit.lm) / deviance(glm(y ~ 1))
      para["amplitude"] <- amplitude(para["cos(ref.phase * w)"], para["sin(ref.phase * w)"])
      para["predicted"] <- phase(para["cos(ref.phase * w)"], para["sin(ref.phase * w)"])
      para["actual"] <- test.metadata[test.metadata$Sample == month, ]$day_of_year %% period
      para["error"] <- circular.diff(para["predicted"], para["actual"])
      test.metadata[test.metadata$Sample == month, names(para)] <- para
    }
    timetable.pred <- rbind(timetable.pred, test.metadata)
  }
  error[[paste0("top", nrow(SOMs))]] <- timetable.pred
}

################################################################################
colors <- colorRampPalette(
  colors = c("#0050FF", "grey50", "#FF5000")
)(length(topN))

timetable.errors <- do.call(rbind, error)
timetable.errors$top <-
  gsub("top([0-9]*)\\..*", "\\1", rownames(timetable.errors)) |>
  as.numeric() |>
  factor()

plsr.pred <- read.csv("results/PLSR.LOO-prediction.csv")
MdAE <- median(abs(plsr.pred$error))
RS.plsr <- with(plsr.pred, rsq(actual, actual + error))

plt_box_timetable_top <- ggplot(data = timetable.errors) +
  geom_hline(yintercept = MdAE, color = "grey20", linetype = "dotted") +
  geom_boxplot(aes(x = top, y = abs(error), fill = top), color = "black", size = 0.5) +
  geom_hline(yintercept = 0, color = "grey20") +
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    limits = c(-10, period / 2), breaks = c(seq(0, 180, 30), MdAE),
    labels = c(seq(0, 180, 30), MdAE |> round(2)), name = "absolute error (day)"
  ) +
  theme_annual() +
  theme(legend.position = "none") +
  labs(x = "top number of metabolites")
plt_box_timetable_top


write.csv(timetable.errors, "figures/timetable/timetable.errors.csv", row.names = F)

plt_RS_timetable_top <-
  group_by(timetable.errors, top) %>%
  summarise(
    "RS" = rsq(actual, actual + error)
  ) %>%
  ggplot() +
  geom_hline(yintercept = RS.plsr) +
  geom_vline(xintercept = c(20, 90), linetype = "dotted") +
  annotate("text",
    x = 12, y = RS.plsr + 0.005,
    label = expression("R"^2 * " of PLSR"),
    size = 4, hjust = 0, vjust = 0
  ) +
  geom_point(aes(x = top |> as.character() |> as.numeric(), y = RS, color = top), size = 3) +
  scale_color_manual(values = colors) +
  scale_x_continuous(
    name = "Top N metabolites", breaks = c(seq(20, 100, 20), seq(120, 200, 20), 250, 300),
    limits = c(0, 316), expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = expression("R"^2 * " of actual vs predicted sampling date"),
    limits = c(0.8, 1), expand = c(0, 0)
  ) +
  theme_annual(x.angle = 45) +
  theme(legend.position = "none")
plt_RS_timetable_top


pdf("figures/timetable/timetable_top-SOMs_RS.pdf", height = 3.5, width = 4.5)
plt_RS_timetable_top
dev.off()

write.csv(plt_RS_timetable_top@data, "figures/timetable/timetable_top-SOMs_RS.csv")


timetable.pred <- error$top20

timetable.pred$error |>
  abs() |>
  quantile()


AUC.timetable <- AUC(abs(timetable.pred$error))
plt_timetable_ECD <- plot_AUC(
  list("timetable" = AUC.timetable),
  colors = c("#3D66C0")
)

plt_timetable_ECD

with(
  timetable.pred,
  rsq(actual, actual + error)
)

SOMs_amide[ids[1:20], ] |> View()


pdf("figures/timetable/timetable_top20_cummulative_error_distribution.pdf", height = 3.5, width = 3.5)
plt_timetable_ECD
dev.off()


write.csv(AUC.timetable[1:2] |> as.data.frame(), "figures/timetable/timetable_top20_cummulative_error_distribution.csv")
