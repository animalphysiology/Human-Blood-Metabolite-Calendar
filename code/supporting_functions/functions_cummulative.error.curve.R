library(ggplot2)
library(patchwork)


# calculate error cummlative curve and its AUC
AUC <- function(error, range = c(0, 365 / 2), accepted.error = 30) {
  MdAE <- median(abs(error))
  n <- length(error)
  accepted.ratio <- (sum(abs(error) < accepted.error)) / n
  x <- (sort(abs(error)) - range[1]) / range[2]
  y <- seq_along(error) / n
  x[x >= 1] <- 1
  x[x <= 0] <- 0
  if (sum(x == 0) < 1) {
    x <- c(0, x)
    y <- c(0, y)
  }
  if (sum(x == 1) < 1) {
    x <- c(x, 1)
    y <- c(y, 1)
  }
  AUC <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

  return(list(
    "x" = x, "y" = y, "AUC" = AUC, "MdAE" = MdAE,
    "accepted.error" = accepted.error,
    "accepted.ratio" = accepted.ratio,
    "interpolated.MdAE" = approx(y, x, xout = 0.5)$y * diff(range),
    "interpolated.accepted.ratio" = approx(x, y, xout = accepted.error / diff(range))$y,
    "interpolated.AUC" = sum(diff(x) * tail(y, -1))
  ))
}


plot_AUC <- function(AUC, colors, linetype = rep(1, length(AUC)), main = 1, text.size = 2.5, interpolate = F) {
  accepted.error <- AUC[[1]]$accepted.error
  accepted.ratio <- ifelse(interpolate, AUC[[1]]$interpolated.accepted.ratio, AUC[[1]]$accepted.ratio)
  MdAE <- ifelse(interpolate, AUC[[1]]$interpolated.MdAE, AUC[[1]]$MdAE)
  plt_AUC <- ggplot() +
    annotate("segment",
      x = 0, xend = 1,
      y = 0, yend = 1, color = "grey50", linewidth = 0.5
    ) +
    scale_x_continuous(
      name = "accepted errors (days)",
      limits = c(0, 1), expand = 0,
      breaks = seq(0, period / 2, 30) / period * 2,
      labels = paste("<", seq(0, period / 2, 30))
    ) +
    scale_y_continuous(name = "propotion of samples", limits = c(0, 1), expand = 0) +
    annotate("text", x = 0.65, y = 0.35, label = "MdAE", size = text.size) +
    annotate("text", x = 0.8, y = 0.35, label = "error < 30", size = text.size) +
    annotate("text", x = 0.95, y = 0.35, label = "AUC", size = text.size) +
    annotate("segment",
      x = accepted.error / period * 2, xend = accepted.error / period * 2,
      y = accepted.ratio, yend = 0, linetype = "dotted", color = "red"
    ) +
    annotate("segment",
      x = accepted.error / period * 2, xend = 0,
      y = accepted.ratio, yend = accepted.ratio,
      linetype = "dotted", color = "red"
    ) +
    annotate("segment",
      x = MdAE / period * 2, xend = MdAE / period * 2,
      y = 0.5, yend = 0, linetype = "dotted"
    ) +
    annotate("segment",
      x = MdAE / period * 2, xend = 0,
      y = 0.5, yend = 0.5, linetype = "dotted"
    ) +
    theme_annual() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  lineform <- ifelse(interpolate, "line", "step")
  add <- list()
  for (n in 1:length(AUC)) {
    AUC.a <- AUC[[n]]
    accepted.error <- AUC.a$accepted.error
    prefix <- if (interpolate) "interpolated." else ""
    accepted.ratio <- AUC.a[[paste0(prefix, "accepted.ratio")]]
    MdAE <- AUC.a[[paste0(prefix, "MdAE")]]
    AUC.value <- AUC.a[[paste0(prefix, "AUC")]]

    ypos <- 0.35 - 0.055 * n

    add <- c(add, list(
      annotate(lineform, x = AUC.a$x, y = AUC.a$y, color = colors[n], linetype = linetype[n]),
      annotate("text", x = 0.525, y = ypos, label = names(AUC)[n], size = text.size),
      annotate("text",
        x = 0.65, y = ypos,
        label = paste0(round(MdAE, 2)), size = text.size
      ),
      annotate("text",
        x = 0.8, y = ypos,
        label = paste0(round(accepted.ratio * 100, 2), "%"), size = text.size
      ),
      annotate("text", x = 0.95, y = ypos, label = round(AUC.value, 3), size = text.size),
      annotate("segment",
        x = 0.45, y = ypos,
        xend = 0.35, yend = ypos, color = colors[n],
        linetype = linetype[n], linewidth = 0.5
      )
    ))
  }

  plt_AUC <- plt_AUC + add
  return(plt_AUC)
}
