library(readxl)
library(tidyr)
library(dplyr)

# frequently used parameter definition
period <- 365
w <- 2 * pi / period

rawdata_file_path <- "data/Human_Amide_2025.txt"
metadata_table_path <- "data/Supplementary_tables.xlsx"

################################################################################
#   load data #
################################################################################

# supplementary data table S1, S2
subject_metadata <- read_excel(metadata_table_path, sheet = 1, skip = 1)

sample_metadata <- read_excel(metadata_table_path, sheet = 2, skip = 1)

sample_metadata$sin.a <- sin(sample_metadata$consecutive_days * w)
sample_metadata$cos.a <- cos(sample_metadata$consecutive_days * w)

sample_metadata$Sex <- factor(sample_metadata$Sex, levels = c("male", "female"))

# remove units
colnames(sample_metadata)[8:9] <- c("PRL", "FT3")

## human_amide_BGI is raw metabolomics quantification table
## annotation part was provided as supplementary table S3
## metabolite abundances of each samples follows
human_amide_BGI <- read.delim(rawdata_file_path)

# make unique ID for metabolite feature
human_amide_BGI$unique_id <- paste(
  "Amide",
  human_amide_BGI$Ion.Mode,
  human_amide_BGI$RT,
  human_amide_BGI$CalcMz,
  sep = "_"
)

# extract abundance data
features <- human_amide_BGI[, sample_metadata$Sample]
rownames(features) <- human_amide_BGI$unique_id


feature.table <- human_amide_BGI[, 1:20]
feature.table$male.log2.mean <- apply(
  log2(features[, subset(sample_metadata, Sex == "male")$Sample]),
  1, mean
)
feature.table$male.log2.sd <- apply(
  log2(features[, subset(sample_metadata, Sex == "male")$Sample]),
  1, sd
)
feature.table$female.log2.mean <- apply(
  log2(features[, subset(sample_metadata, Sex == "female")$Sample]),
  1, mean
)
feature.table$female.log2.sd <- apply(
  log2(features[, subset(sample_metadata, Sex == "female")$Sample]),
  1, sd
)

rownames(feature.table) <- feature.table$unique_id

## table S3
write.csv(feature.table, "data/feature.table.csv", row.names = F, na = "")

subjects <- LETTERS[1:11]
males <- LETTERS[1:6]
females <- LETTERS[7:11]

### objects for plot
days.of.month <- c(
  "Jan" = 31, "Feb" = 28, "Mar" = 31, "Apr" = 30,
  "May" = 31, "Jun" = 30, "Jul" = 31, "Aug" = 31,
  "Sep" = 30, "Oct" = 31, "Nov" = 30, "Dec" = 31
)

mid.month <- cumsum(days.of.month) - days.of.month / 2

first.month <- cumsum(days.of.month) - days.of.month


center.of.sampling <- group_by(sample_metadata, Month) %>%
  summarise(center.month = mean(consecutive_days)) %>%
  .$center.month %>%
  unlist() %>%
  sort()

names(center.of.sampling) <- c(
  "Apr", "May", "Jun", "Jul", "Aug", "Sep",
  "Oct", "Nov", "Dec", "Jan", "Feb", "Mar"
)


month.label <- c(first.month[4:12], first.month[1:4] + 365)
names(month.label) <- month.abb[c(4:12, 1:4)]


Daylength_data <- read_excel("data/meteorological_data/2024/2024-2025_Daylength.xlsx")

meteorological_data <- list()
for (f in grep(".csv", list.files("data/meteorological_data/2024"), value = T)) {
  meteorological_data[[gsub("_data.csv", "", f)]] <-
    read.csv(paste0("data/meteorological_data/2024/", f), skip = 0, header = F, fileEncoding = "CP932")
}

meteorological_data <- lapply(meteorological_data, function(raw) {
  df <- raw[(nrow(raw) - 10943):nrow(raw), ]
  colnames(df) <- apply(raw[2:(nrow(raw) - 10945), ], 2, paste0, collapse = "")
  df$date <- as.POSIXct(df[[1]], format = "%Y/%m/%d %H:%M", tz = "Asia/Tokyo")
  return(df)
})

# sliding.windows (day)
window <- 30
library(caTools)

environmental_data_2024 <- data.frame(
  date = meteorological_data$`2024_2025_AirTempareture`$date,
  consecutive_days = as.numeric(as.POSIXct(meteorological_data$`2024_2025_AirTempareture`$date) - as.POSIXct("2024-01-01 00:00:00")) / 24
)

environmental_data_2024$AirTemperature <-
  meteorological_data$`2024_2025_AirTempareture`[[2]]

environmental_data_2024$AirTemperature_monthly_runmean <-
  runmean(meteorological_data$`2024_2025_AirTempareture`[[2]],
    k = 24 * window,
    align = "c", endrule = "m"
  )

environmental_data_2024$Daylength <-
  Daylength_data$Daylength[match(
    as.POSIXct(environmental_data_2024$date),
    as.POSIXct(Daylength_data$Date)
  )]

environmental_data_2024$SunshineDuration <-
  meteorological_data$`2024_2025_SunshineDuration`[[2]]

environmental_data_2024$SunshineDuration_monthly_runmean <-
  runmean(meteorological_data$`2024_2025_SunshineDuration`[[2]],
    k = 24 * window, align = "c", endrule = "m"
  )


environmental_data_2024[, -1] <- apply(environmental_data_2024[, -1], 2, as.numeric)

range(environmental_data_2024$AirTemperature, na.rm = T)
# [1]  -2.1 38.7

range(environmental_data_2024$Daylength, na.rm = T)
# [1]  9.783333 14.533333

cor.test(environmental_data_2024$Daylength,
  environmental_data_2024$SunshineDuration_monthly_runmean,
  use = "complete"
)
# Pearson's product-moment correlation
#
# data:  environmental_data_2024$Daylength and environmental_data_2024$SunshineDuration_monthly_runmean
# t = 10.178, df = 454, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3531743 0.5029447
# sample estimates:
#       cor
# 0.4310234

# remove intermediates
rm(list = c(
  "meteorological_data", "Daylength_data", "f", "window",
  "metadata_table_path", "rawdata_file_path"
))

################################################################################
#### functions for plotting data
################################################################################
library(ggplot2)

subject.color <- function(subject.set, a = "AA") {
  subject.color <- c(
    "#0000FF", # Blue
    "#4169E1", # Royal Blue
    "#191970", # Midnight Blue
    "#6A5ACD", # Slate Blue
    "#40E0D0", # Turquoise
    "#1E90FF", # Dodger Blue
    "#FF0000", # Red
    "#FF4500", # Orange Red
    "#C71585", # Medium Violet Red
    "#FF69B4", # Hot Pink
    "#A0522D", # Sienna
  ) |>
    paste0(a)

  names(subject.color) <- LETTERS[1:11]
  return(subject.color[subject.set])
}

subject.shape <- function(subject.set, fill = F) {
  if (fill) {
    subject.shape <- c(rep(c(21:23), each = 1, times = 4)[1:11], 24:25)
  } else {
    subject.shape <- c(rep(c(1, 5, 0), each = 1, times = 4)[1:11], 2, 6)
  }
  names(subject.shape) <- LETTERS[1:13]
  return(subject.shape[subject.set])
}

scale_Subject <- function(subject.set = LETTERS[1:13], alpha = "70", fill = T) {
  list(
    scale_color_manual(values = subject.color(subject.set, alpha)),
    scale_fill_manual(values = subject.color(subject.set, alpha)),
    scale_shape_manual(values = subject.shape(subject.set, fill = fill))
  )
}

annotate.summer.winter <- function(ylim, ratio = 0.02) {
  list(
    annotate("rect", ymax = ylim[1], ymin = ylim[1] + diff(ylim) * ratio, xmax = month.label[3], xmin = month.label[6], fill = "#FF800080"),
    annotate("rect", ymax = ylim[1], ymin = ylim[1] + diff(ylim) * ratio, xmax = month.label[9], xmin = month.label[12], fill = "#0080FF80")
  )
}

theme_annual <- function(x.angle = 0, legend = "top") {
  list(theme_classic(), theme(
    axis.text.x = element_text(
      angle = x.angle,
      hjust = ifelse(x.angle == 0, 0.5, 1),
      vjust = ifelse(x.angle == 0, 1, 1)
    ),
    plot.title = element_text(hjust = 0.5),
    legend.key.width = unit(0.25, "cm"),
    legend.position = legend,
    panel.background = element_blank(),
    legend.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank()
  ), if (legend == "right") {
    guides(color = guide_legend(ncol = 1, label.position = legend))
  } else if (legend == "top") {
    guides(color = guide_legend(nrow = 1, label.position = legend))
  } else {
    NULL
  })
}


plot_metabolite <- function(metabolites, centered = T) {
  plot.name <- metabolites$unique_id
  names(plot.name) <- metabolites$curated_name
  plot.list <- lapply(plot.name, function(metabolite) {
    relative.abundance <- unlist(features[metabolite, ]) |> log2()
    if (centered) {
      y <- lapply(LETTERS[1:11], function(i) {
        y <- relative.abundance[grep(i, names(relative.abundance))]
        return(y - mean(y, na.rm = T))
      }) |> unlist()
      relative.abundance <- y[names(relative.abundance)]
    }

    A <- metabolites[metabolite, "amplitude"]
    P <- metabolites[metabolite, "phase"]
    M.f <- metabolites[metabolite, "mesor.female"]
    M.m <- metabolites[metabolite, "mesor.male"]

    ylim <- c(min(relative.abundance, na.rm = T), max(relative.abundance, na.rm = T))
    ylim <- ylim + c(diff(ylim) * -0.1, diff(ylim) * 0.1)

    plt_p <- ggplot(data = data.frame(
      "timepoints" = sample_metadata$consecutive_days,
      "Subject" = sample_metadata$Subject,
      relative.abundance
    ) |> na.omit()) +
      annotate.summer.winter(ylim, 0.05) +
      geom_point(aes(
        x = timepoints, y = relative.abundance,
        color = Subject, fill = Subject, shape = Subject
      ), size = 1.5) +
      scale_x_continuous(
        name = "", expand = c(0, 0),
        limits = c(90, 455),
        breaks = month.label,
        labels = names(month.label)
      ) +
      geom_vline(xintercept = ifelse(P <= 90, P + 365, P), linewidth = .01, color = "grey20") +
      scale_y_continuous(name = "", limits = ylim, expand = c(0, 0)) +
      scale_Subject() +
      stat_function(
        fun = function(x) {
          A / 2 * cos((x - P) * w) + mean(relative.abundance, na.rm = T)
        },
        color = "black", linewidth = 0.8, n = period
      ) +
      theme_classic() +
      theme_annual() +
      theme(
        legend.position = "n",
        axis.text.x = element_blank()
      )
    return(plt_p + labs(title = metabolites[metabolite, "curated_name"]))
  })
  return(plot.list)
}


save.image(file = "data/Human_annual_metabolomics.RData")
