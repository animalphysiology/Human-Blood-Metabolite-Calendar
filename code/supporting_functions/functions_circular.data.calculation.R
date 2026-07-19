period <- 365
w <- 2 * pi / period

## function for min-max normalization to -1 to 1
min.max.normalize <- function(y, range = c(-1, 1)) {
  2 * (y - min(y)) / max(y - min(y)) - 1
}

## function for min-max normalization to 0 to 1
zero.one.normalize <- function(y) {
  (y - min(y)) / max(y - min(y))
}

## standard error
se <- function(x, na.rm = T) {
  n <- ifelse(na.rm, sum(!is.na(x)), length(x))
  se <- sd(x, na.rm = T) / sqrt(n)
  return(se)
}

amplitude <- function(cos, sin) {
  A <- sqrt(cos^2 + sin^2)
  return(A)
}

phase <- function(cos, sin, period = 365) {
  return((atan2(sin, cos) * period / (2 * pi)) %% period)
}

rsq <- function(x, y) {
  cor(x, y)^2
}

## calculate amplitude and phase from cosinor coefficients
## input : cosinor=c(intercept,sin.coef,cos.coef)
## output : mesor, sin, cos, amp, phase
cos.params <- function(cosinor, period = 365, suffix = "", prefix = "", omit.original = F) {
  out <- numeric()
  if (!omit.original) {
    out[paste0(prefix, "mesor", suffix)] <- cosinor[1]
    out[paste0(prefix, "sin", suffix)] <- cosinor[2]
    out[paste0(prefix, "cos", suffix)] <- cosinor[3]
  }
  out[paste0(prefix, "amplitude", suffix)] <- sqrt(sum((cosinor[c(2, 3)])^2)) * 2
  out[paste0(prefix, "phase", suffix)] <- (atan2(cosinor[2], cosinor[3]) * period / (2 * pi)) %% period
  return(out)
}


## function to compute the circular mean of angles
circular.mean <- function(angles, weights = NULL, period = 365) {
  if (is.null(weights)) {
    weights <- rep(1, length(angles))
  }
  sin.sum <- sum(sin(angles * 2 * pi / period) * weights)
  cos.sum <- sum(cos(angles * 2 * pi / period) * weights)
  mean.angle <- (atan2(sin.sum, cos.sum) * period / (2 * pi)) %% period
  return(mean.angle)
}


## function to compute the circular mean of angles
circular.var <- function(angles, weights = NULL, period = 365) {
  n <- length(angles)
  if (is.null(weights)) {
    weights <- rep(1, n)
  }
  sin.sum <- sum(sin(angles * 2 * pi / period) * weights) / n
  cos.sum <- sum(cos(angles * 2 * pi / period) * weights) / n
  var.angle <- sum(c(sin.sum, cos.sum)^2) |> sqrt()
  return((1 - var.angle) * period)
}


circular.sd <- function(angles, weights = NULL, period = 365) {
  n <- length(angles)
  if (is.null(weights)) {
    weights <- rep(1, n)
  }
  sin.sum <- sum(sin(angles * 2 * pi / period) * weights) / n
  cos.sum <- sum(cos(angles * 2 * pi / period) * weights) / n
  var.angle <- sum(c(sin.sum, cos.sum)^2) |> sqrt()
  return(sqrt(-2 * log(var.angle)) / (2 * pi) * period)
}


circular.two.mean <- function(angles1, angles2, period = 365) {
  mean.angles <- apply(data.frame(angles1, angles2), 1, circular.mean, period = period)
  return(mean.angles)
}

circular.two.sd <- function(angles1, angles2, period = 365) {
  var.angles <- apply(data.frame(angles1, angles2), 1, circular.sd, period = period)
  return(var.angles)
}


# calculate signed angle differences
circular.diff <- function(angles1, angles2, period = 365) {
  angles1 <- angles1 %% period
  angles2 <- angles2 %% period
  phase.diff <- (angles1 - angles2)
  phase.diff <- ifelse(phase.diff > period / 2, phase.diff - period, phase.diff)
  phase.diff <- ifelse(phase.diff < -1 * period / 2, phase.diff + period, phase.diff)
  return(phase.diff)
}

circular.cor <- function(angles1, angles2, period = 365) {
  w <- 2 * pi / period
  rm.na <- complete.cases(angles1, angles2)
  angles1 <- angles1[rm.na]
  angles2 <- angles2[rm.na]
  mean.angles1 <- circular.mean(angles1)
  mean.angles2 <- circular.mean(angles2)

  # Jammalamadaka–Sengupta correlation
  num <- sum(sin((angles1 - mean.angles1) * w) *
    sin(angles2 - mean.angles2) * w)
  den <- sqrt(sum(sin((angles1 - mean.angles1) * w)^2) *
    sum(sin((angles1 - mean.angles1) * w)^2))

  rho <- num / den
  return(rho)
}


RS.phase <- function(cos.p, sin.p, cos.a, sin.a) {
  error <- (circular.diff(phase(cos.a, sin.a), phase(cos.p, sin.p)))
  RS <- rsq(phase(cos.a, sin.a) + error, phase(cos.a, sin.a))
  return(RS)
}



# counting number of samples by phase within sliding windows
circular.SW.count <- function(phase, bins, period = 365, range = c(0, period)) {
  w.bins <- period / bins
  seq.bin <- seq(min(range), max(range), length.out = bins + 1)
  counts <- sapply(1:bins, function(i) {
    sum(phase >= seq.bin[i] & phase < seq.bin[i + 1])
  })
  seq.bin <- round(seq.bin, 10)
  phase.bins <- data.frame(
    row.names = paste0("[", seq.bin[1:bins], ",", seq.bin[1:bins + 1], "]"),
    "lower" = seq.bin[1:bins],
    "upper" = seq.bin[1:bins + 1],
    "center" = (seq.bin[1:bins] + seq.bin[1:bins + 1]) / 2,
    "count" = counts
  )
  return(phase.bins)
}


# color gradation of circular phase
circular.phase.color <-
  function(phase, period = 365) {
    w <- 2 * pi / period
    cols <- rgb(1 - (cos(phase * w) + 1) / 2, (sin(phase * w) + 1) / 8 + 0.5, (cos(phase * w) + 1) / 2)
    return(cols)
  }
