## ----decayGraph, echo = FALSE, fig.height = 3, fig.width = 6, fig.caption = "Exponential decay model with a decay rate, alpha, of -2.3, minCutoff = 1e-3 and distCutoff given by the vertical dotted lines."----
xx <- seq(0, 4, length = 100)
fn <- function(edgeDists, minCov = 1e-3, distCutoff = 2) {
    q1 <- (minCov - 1)/(exp(-2.3 * distCutoff) - 1)
    q2 <- 1 - q1
    ans <- (q2 + q1 * exp(-2.3 * edgeDists))
    ans[ans < minCov] <- 0
    return(ans)
}
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
plot(xx, fn(xx), type = "l", las = 1,
     xlab = "Distance",
     ylab = "Covariance")
abline(v = 2, lty = 2, lwd = 0.5)
plot(xx, fn(xx, distCutoff = 1), type = "l", las = 1,
     xlab = "Distance",
     ylab = "Covariance")
abline(v = 1, lty = 2, lwd = 0.5)

