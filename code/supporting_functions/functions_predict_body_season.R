timetable <- function(SOMs, metadata, features, plot = T, period = 365, phase.col = "phase") {
  w <- 2 * pi / period
  timetable.pred <- data.frame(metadata)
  # rn <- intersect(SOMs$unique_id, rownames(features))
  # SOMs <- SOMs[rn, ]
  norm.data <- apply(features |> log2(), 1, scale) |> t()
  dimnames(norm.data) <- dimnames(features)
  norm.data <- norm.data[SOMs$unique_id, ]
  plots <- list()
  for (month in metadata$Sample) {
    y <- norm.data[, month] |> unlist()
    x <- SOMs[, phase.col]
    fit.lm <- glm(y ~ sin(x * w) + cos(x * w))
    para <- coef(fit.lm)
    para["RS"] <- 1 - deviance(fit.lm) / deviance(glm(y ~ 1))
    para["amplitude"] <- amplitude(para["cos(x * w)"], para["sin(x * w)"])
    para["predicted"] <- phase(para["cos(x * w)"], para["sin(x * w)"])
    para["actual"] <- metadata[metadata$Sample == month, ]$day_of_year
    para["error"] <- circular.diff(para["predicted"], para["actual"])

    timetable.pred[metadata$Sample == month, names(para)] <- para


    if (plot) {
      ylim <- c(min(y) |> floor(), max(y) |> ceiling())
      plots[[month]] <-
        ggplot(data = data.frame(y = y, x = x)) +
        geom_point(aes(x = x, y = y, color = x)) +
        geom_vline(xintercept = para["actual"], linetype = 2, color = "red") +
        geom_vline(xintercept = para["predicted"], linetype = 2, color = "black") +
        scale_color_gradientn(colors = circular.phase.color(1:365)) +
        stat_function(
          fun = eval(substitute(
            function(x) {
              A * cos((x - P) * W) + B
            },
            list(
              A = para["amplitude"],
              P = para["predicted"],
              W = w,
              B = para["(Intercept)"]
            )
          )),
          color = "black", linewidth = 0.8, n = period
        ) +
        scale_x_continuous(
          name = "reference peak phase of metabolites",
          breaks = month.label %% period,
          labels = (month.label %% period) |> names(),
          limits = c(0, period), expand = 0
        ) +
        scale_y_continuous(name = "relative abundance", limits = ylim, expand = 0) +
        theme_annual(x.angle = 45) +
        theme(legend.position = "n") +
        labs(title = paste(
          month.name[as.numeric(sub(".*_", "", month), "")],
          # "sample from subject", sub("_.*", "", month),
          ""
        ))
    }
  }
  if (plot) plots else timetable.pred
}
