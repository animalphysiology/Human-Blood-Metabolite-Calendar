library(lmerTest)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(pheatmap)

source("code/supporting.functions.R")
load("data/Human_annual_metabolomics.RData")

################################################################################
#   Environmental data plotting #
################################################################################

ylim <- c(-5, 40)
AT_2024 <- ggplot(data = environmental_data_2024) +
  geom_line(aes(x = consecutive_days, y = AirTemperature), color = "#50A0FFAA") +
  geom_line(aes(x = consecutive_days, AirTemperature_monthly_runmean),
    color = "#50A0FF", linewidth = 1
  ) +
  scale_y_continuous(
    name = "Ambient temperature (°C)", limits = ylim, expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(90, 455),
    breaks = month.label,
    labels = names(month.label),
    expand = c(0, 0), name = "month",
    sec.axis =
      sec_axis(~.,
        breaks = c(seq(90, 330, 30), seq(365, 455, 30)),
        labels = c(seq(90, 330, 30), seq(365, 455, 30)) %% period,
        name = "day-of-year"
      )
  ) +
  annotate.summer.winter(ylim) +
  theme_annual() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.x.top = element_text(angle = 45, hjust = 0))


ylim <- c(9, 15)
DL_2024 <- ggplot(data = environmental_data_2024 |> na.omit()) +
  geom_line(aes(x = consecutive_days, y = Daylength), color = "#FFAA00", linewidth = 1) +
  scale_y_continuous(
    name = "Day length (hour)", limits = ylim, expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(90, 455),
    breaks = month.label,
    labels = names(month.label),
    expand = c(0, 0), name = "month",
    sec.axis =
      sec_axis(~.,
        breaks = c(seq(90, 330, 30), seq(365, 455, 30)),
        labels = c(seq(90, 330, 30), seq(365, 455, 30)) %% period,
        name = "day-of-year"
      )
  ) +
  annotate.summer.winter(ylim) +
  theme_annual() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.x.top = element_text(angle = 45, hjust = 0))

DL_2024 / AT_2024 # Fig1A


################################################################################
#   hormone data analysis #
################################################################################

# parameter definition
period <- 365
w <- 2 * pi / period

hormone.seasonality <- data.frame("hormone" = c("Prolactin", "Free T3"))

##### prolactin ######

fit.lmm.PRL <- lmer(
  log2(PRL) ~ cos.a + sin.a + Sex + (1 | Subject),
  data = sample_metadata, subset = log2(PRL) < log2(50), REML = F
)

fit.lmm.PRL


fixed <- fixef(fit.lmm.PRL)


A.PRL <- amplitude(fixed[2], fixed[3]) * 2
P.PRL <- phase(fixed[2], fixed[3])
M.PRL.m <- fixed[1]
M.PRL.f <- fixed[1] + fixed[4]

M.PRL <- mean(M.PRL.m, M.PRL.f)

contrast <- 2:3


significance.PRL <- contest(
  fit.lmm.PRL,
  diag(length(fixef(fit.lmm.PRL)))[contrast, ]
)[1, 6]


hormone.seasonality[1, c("cos", "sin")] <- fixed[2:3]
hormone.seasonality[1, c("amplitude", "phase")] <- c(A.PRL, P.PRL)
hormone.seasonality[1, c("mesor.male", "mesor.female")] <- c(M.PRL.m, M.PRL.f)
hormone.seasonality[1, "p.value"] <- significance.PRL


human_PRL <- merge(sample_metadata,
  coef(fit.lmm.PRL)$Subject[c(1, 4)],
  by.x = "Subject", by.y = 0
)

human_PRL$baseline <-
  ifelse(human_PRL$Sex == "male",
    human_PRL$`(Intercept)`,
    human_PRL$`(Intercept)` + human_PRL$Sexfemale
  )

human_PRL$Sexfemale <- NULL
human_PRL$`(Intercept)` <- NULL

human_PRL$relative.PRL <- log2(human_PRL$PRL) - human_PRL$baseline


ylim <- c(-1, 1)

plt_PRL <- ggplot(data = subset(human_PRL, PRL < 50)) +
  geom_vline(xintercept = P.PRL, color = "grey20", linewidth = 0.5) +
  annotate.summer.winter(ylim) +
  geom_point(aes(x = consecutive_days, y = relative.PRL, color = Subject, fill = Subject, shape = Subject), size = 1.5) +
  stat_function(
    fun = function(x) {
      A.PRL / 2 * cos((x - P.PRL) * w)
    },
    color = "black", linewidth = 0.8, n = period
  ) +
  scale_y_continuous(name = "relative concentration (log2 ng/mL)", limits = ylim, expand = c(0, 0)) +
  scale_x_continuous(
    name = "", limits = c(90, 455), expand = c(0, 0),
    breaks = month.label, labels = names(month.label)
  ) +
  scale_Subject(fill = T) +
  annotate("text",
    x = 100, y = max(ylim) * 0.9, hjust = 0, size = 4,
    label = paste("p =", signif(significance.PRL, 3))
  ) +
  theme_annual(legend = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Prolactin")

P.PRL

plt_PRL


##### Free T3


fit.lmm.FT3 <- lmer(
  log2(FT3) ~ (cos.a + sin.a) + Sex + (1 | Subject),
  data = sample_metadata, subset = log2(FT3) < log2(4.5), REML = F
)

fit.lmm.FT3

fixed <- fixef(fit.lmm.FT3)


A.FT3 <- amplitude(fixed[2], fixed[3]) * 2
P.FT3 <- phase(fixed[2], fixed[3])
M.FT3.m <- fixed[1] + fixed[4]
M.FT3.f <- fixed[1]

M.FT3 <- mean(M.FT3.m, M.FT3.f)


contrast <- 2:3

significance.FT3 <- contest(
  fit.lmm.FT3,
  diag(length(fixef(fit.lmm.FT3)))[contrast, ]
)[1, 6]


hormone.seasonality[2, c("cos", "sin")] <- fixed[2:3]
hormone.seasonality[2, c("amplitude", "phase")] <- c(A.FT3, P.FT3)
hormone.seasonality[2, c("mesor.male", "mesor.female")] <- c(M.FT3.m, M.FT3.f)
hormone.seasonality[2, "p.value"] <- significance.FT3


human_FT3 <- merge(sample_metadata,
  coef(fit.lmm.FT3)$Subject[c(1, 4)],
  by.x = "Subject", by.y = 0
)


human_FT3$baseline <-
  ifelse(human_FT3$Sex == "male",
    human_FT3$`(Intercept)`,
    human_FT3$`(Intercept)` + human_FT3$Sexfemale
  )

human_FT3$Sexfemale <- NULL
human_FT3$`(Intercept)` <- NULL

human_FT3$relative.FT3 <- log2(human_FT3$FT3) - human_FT3$baseline


ylim <- c(-0.25, 0.25)

plt_FT3 <- ggplot(data = subset(human_FT3, FT3 < 4.5)) +
  geom_vline(xintercept = P.FT3, color = "grey20", linewidth = 0.5) +
  annotate.summer.winter(ylim) +
  # geom_line(aes(x=consecutive_days, y=relative.FT3, color=Subject),size=0.1)+
  geom_point(aes(x = consecutive_days, y = relative.FT3, color = Subject, fill = Subject, shape = Subject), size = 1.5) +
  stat_function(
    fun = function(x) {
      A.FT3 / 2 * cos((x - P.FT3) * w)
    },
    color = "black", linewidth = 0.8, n = period
  ) +
  scale_y_continuous(name = "relative concentration (log2 pg/mL)", limits = ylim, expand = c(0, 0)) +
  scale_x_continuous(
    name = "", limits = c(90, 455), expand = c(0, 0),
    breaks = month.label, labels = names(month.label)
  ) +
  scale_Subject(fill = T) +
  annotate("text",
    x = 100, y = max(ylim) * 0.9, hjust = 0, size = 4,
    label = paste("p =", signif(significance.FT3, 3))
  ) +
  theme_annual(legend = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Free T3")

P.FT3

plt_FT3


##########

# Fig 1B
pdf("figures/hormone/PRL-FT3_cosinor.pdf", width = 8, height = 3)
plt_PRL + theme(legend.position = "none") +
  plt_FT3 +
  theme(legend.position = "right", legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.4, "cm")) +
  guides(color = guide_legend(ncol = 1, label.position = "right"))
dev.off()





################################################################################

Subject <- LETTERS[1:11]
Sex <- c(rep("male", 6), rep("female", 5)) |> factor(levels = c("male", "female"))
baseline.FT3 <- unique(human_FT3$baseline)
baseline.PRL <- unique(human_PRL$baseline)


hormone.seasonality[1, paste0("baseline.", LETTERS[1:11])] <- baseline.PRL
hormone.seasonality[2, paste0("baseline.", LETTERS[1:11])] <- baseline.FT3

write.csv(hormone.seasonality, "results/hormone.seasonality.csv", row.names = F)


ylim <- c(2, 5.5)

plt_baseline_PRL <- ggplot() +
  geom_point(aes(x = Sex, y = baseline.PRL, color = Subject, fill = Subject, shape = Subject),
    position = position_dodge(width = 0.5), size = 2
  ) +
  stat_summary(aes(x = Sex, y = baseline.PRL), fun.data = mean_sdl, color = "grey40", shape = 95, linewidth = 0.3) +
  scale_Subject(fill = T) +
  scale_y_continuous(
    name = "Prolactin baseline (log2 ng/mL)",
    limits = ylim, expand = c(0, 0)
  ) +
  theme_annual()


ylim <- c(1, 2)

plt_baseline_FT3 <- ggplot() +
  geom_point(aes(x = Sex, y = baseline.FT3, color = Subject, fill = Subject, shape = Subject),
    position = position_dodge(width = 0.5), size = 2
  ) +
  stat_summary(aes(x = Sex, y = baseline.FT3), fun.data = mean_sdl, color = "grey40", shape = 95, linewidth = 0.3) +
  scale_Subject(fill = T) +
  scale_y_continuous(
    name = "Free T3 baseline (log2 pg/mL)",
    limits = ylim, expand = c(0, 0)
  ) +
  theme_annual()


# fig 1C
pdf("figures/hormone/PRL_FT3_baseline.pdf", width = 5, height = 3)
plt_baseline_PRL + theme(legend.position = "n") + plt_baseline_FT3 +
  theme(
    legend.position = "right",
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.4, "cm")
  ) +
  guides(color = guide_legend(ncol = 1, label.position = "right"))
dev.off()
