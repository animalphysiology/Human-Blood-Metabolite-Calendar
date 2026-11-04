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

LOO_lmm_cosinor <- LOO_lmm_cosinor |>
  lapply(function(df) {
    df$unique_id <- rownames(df)
    df
  })

plsr.models <- readRDS("results/PLSR_calendar_models.rds")

################################################################################

topN <- c(seq(10, 50, 5), seq(60, 100, 10), seq(120, 200, 20))


set.name <- "plsr_selected"

error <- list()

for (n in topN) {
  timetable.pred <- data.frame()
  for (subject.out in LETTERS[1:11]) {
    sum.effect.size <- with(
      plsr.models[[subject.out]],
      sqrt(abs(c.model$opt.coefficients)^2 + abs(s.model$opt.coefficients)^2)
    ) |> sort(decreasing = T)

    SOMs <- LOO_lmm_cosinor[[subject.out]] |> subset(q.value < 0.05)
    SOMs <- SOMs[names(sum.effect.size)[1:n], ]

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
      y <- test.features[rownames(SOMs), month] |> unlist()
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

plt_box_timetable_top <-
  ggplot(data = timetable.errors) +
  geom_hline(yintercept = MdAE, color = "grey20", linetype = "dotted") +
  geom_boxplot(aes(x = top, y = abs(error), fill = top), color = "black", size = 0.5) +
  geom_hline(yintercept = 0, color = "grey20") +
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    limits = c(-10, period / 2), breaks = c(seq(0, 180, 30), MdAE),
    labels = c(seq(0, 180, 30), MdAE |> round(2)), name = "absolute error (day)"
  ) +
  theme_annual(x.angle = 45) +
  theme(legend.position = "none") +
  labs(x = "top number of metabolites")


plt_RS_timetable_top <-
  group_by(timetable.errors, top) %>%
  summarise(
    "RS" = rsq(actual, actual + error)
  ) %>%
  ggplot() +
  geom_vline(xintercept = c(30), linetype = "dotted") +
  geom_point(aes(x = top |> as.character() |> as.numeric(), y = RS, color = top), size = 3) +
  scale_color_manual(values = colors) +
  scale_x_continuous(
    name = "Top N metabolites", breaks = c(seq(20, 100, 20), seq(120, 200, 20)),
    limits = c(0, 205), expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = expression("R"^2 * " of actual vs estimated"),
    limits = c(0.7, 1), expand = c(0, 0)
  ) +
  theme_annual(x.angle = 45) +
  theme(legend.position = "none")

plt_MdAE_timetable_top <-
  group_by(timetable.errors, top) %>%
  summarise(
    "MdAE" = abs(error) |> median()
  ) %>%
  ggplot() +
  geom_vline(xintercept = c(30), linetype = "dotted") +
  geom_point(aes(x = top |> as.character() |> as.numeric(), y = MdAE, color = top), size = 3) +
  scale_color_manual(values = colors) +
  scale_x_continuous(
    name = "Top N metabolites", breaks = c(seq(20, 100, 20), seq(120, 200, 20)),
    limits = c(0, 205), expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "MdAE of actual vs estimated",
    limits = c(0, 40), expand = c(0, 0)
  ) +
  theme_annual(x.angle = 45) +
  theme(legend.position = "none")

plt_RS_timetable_top + plt_MdAE_timetable_top # + plt_box_timetable_top



pdf("figures/timetable/timetable_top-SOMs_boxplot.pdf", height = 3.5, width = 4.5)
plt_box_timetable_top
dev.off()

write.csv(timetable.errors, "figures/timetable/timetable.errors.csv", row.names = F)


pdf("figures/timetable/timetable_top-SOMs_RS.pdf", height = 3.5, width = 4.5)
plt_RS_timetable_top

pdf("figures/timetable/timetable_top-SOMs_RS-MdAE.pdf", height = 3.5, width = 8)
plt_RS_timetable_top + plt_MdAE_timetable_top
dev.off()


write.csv(plt_RS_timetable_top@data, "figures/timetable/timetable_top-SOMs_RS.csv")


timetable.pred <- error$top30

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


pdf("figures/timetable/timetable_top30_cummulative_error_distribution.pdf", height = 3.5, width = 3.5)
plt_timetable_ECD + coord_fixed()
dev.off()


write.csv(AUC.timetable[1:2] |> as.data.frame(), "figures/timetable/timetable_top20_cummulative_error_distribution.csv")


n <- 30
subject.out <- "F"
sum.effect.size <- with(
  plsr.models[[subject.out]],
  abs(c.model$opt.coefficients) + abs(s.model$opt.coefficients)
) |> sort(decreasing = T)
SOMs <- LOO_lmm_cosinor[[subject.out]] |> subset(q.value < 0.05)
SOMs <- SOMs[names(sum.effect.size)[1:n], ]

pdf("figures/timetable/examples_Subject_F.pdf", width = 14, height = 3)
timetable(
  SOMs,
  subset(sample_metadata, Month %in% c(3, 6, 9, 12) & Subject == "A"),
  features
) |> wrap_plots(ncol = 4)
dev.off()


pdf("figures/timetable/examples_Subject_F.pdf", width = 6, height = 5)
timetable(
  SOMs,
  subset(sample_metadata, Month %in% c(3, 6, 9, 12) & Subject == "A"),
  features
) |> wrap_plots(ncol = 2)
dev.off()
