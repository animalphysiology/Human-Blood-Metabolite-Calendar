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


################################################################################
#   metabolome data fitting with LOO #
################################################################################


###### leave-one-subject-out ######

LOO_lmm_cosinor.result <- list()

for (subject.out in subjects) {
  # fitting

  print(subject.out)
  lmm_cosinor.result <- data.frame(row.names = human_amide_BGI$unique_id)
  sample_metadata_sub <- subset(sample_metadata, Subject != subject.out)

  features <- human_amide_BGI[, sample_metadata_sub$Sample]
  rownames(features) <- human_amide_BGI$unique_id

  ## outlier removal
  norm.data <- apply(features, 1, function(x) log2(x) |> scale()) |> t()
  sum(abs(norm.data) > 3)
  print(paste(
    sum(abs(norm.data) > 3) / prod(dim(norm.data)) * 100,
    "% value with |STDev| > 3 removed as outlier"
  ))

  features[abs(norm.data) > 3] <- NA


  for (m in rownames(features)) {
    y <- features[m, ] |>
      log2() |>
      unlist()
    suppressMessages({
      fit.cosinor <- lmer(y ~ sin.a + cos.a + Sex + (1 | Subject),
        data = sample_metadata_sub, REML = F
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


  LOO_lmm_cosinor.result[[subject.out]] <- lmm_cosinor.result
}

LOO_lmm_cosinor.result <- lapply(LOO_lmm_cosinor, function(df) {
  df$unique_id <- rownames(df)
  df
})


saveRDS(
  LOO_lmm_cosinor.result,
  "results/LOO_human_amide_lmm_cosinor_annual.rds"
)


LOO_SOMs_Amide <- LOO_lmm_cosinor.result |>
  lapply(function(df) {
    subset(df, q.value < 0.05)
  })

seasonal.metabolites <- readRDS("results/human_Amide_lmm_cosinor_annual.rds")
row.names(seasonal.metabolites) <- seasonal.metabolites$unique_id

SOMs_LOO_common <- seasonal.metabolites[Reduce(
  intersect,
  LOO_SOMs_Amide |>
    lapply(function(df) {
      rownames(df)
    })
), ]
