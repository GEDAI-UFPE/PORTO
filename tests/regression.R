library(PORTO)
tol <- 1e-6
ok <- function(f) stopifnot(all(f$solver_status == 0L))
Y <- matrix(c(5, 8, 8, 5, 6, 6), 3, 2, byrow = TRUE)
R <- matrix(rep(c(1, 2), 3), 3, 2, byrow = TRUE)
u <- dea.order(Y = Y)
w <- dea.order(Y = Y, OY = R)
m <- dea.order(Y = Y, OY = R, Gamma = 0.5)
s <- dea.order(Y = Y, OY = R, Weak = FALSE)
for (f in list(u, w, m, s)) ok(f)
stopifnot(ncol(w$weights_input) == 1L, w$auto_unit_input,
          all(u$efficiency + tol >= w$efficiency),
          all(w$efficiency + tol >= m$efficiency),
          all(m$efficiency + tol >= s$efficiency),
          max(abs(s$virtual_output - s$surrogate_output)) < tol,
          max(abs(s$efficiency - dea.order(Y = Y, OY = R,
                                          Gamma = 1)$efficiency)) < tol)
# One-dimensional hand-verifiable BoD: x=1, y=(1,2,4).
f <- dea.order(Y = matrix(c(1, 2, 4), ncol = 1))
ok(f)
stopifnot(max(abs(f$efficiency - c(.25, .5, 1))) < tol)
# Reciprocal one-output model: higher inputs are worse.
f <- dea.order(X = matrix(c(1, 2, 4), ncol = 1), ORIENTATION = "in")
ok(f)
stopifnot(ncol(f$weights_output) == 1L,
          max(abs(f$efficiency - c(1, .5, .25))) < tol)
# Pooled ties across all five surrogate schemes.
T <- matrix(c(2, 2, 1, 1, 1.8, 2.1, 1.2, .9,
              2.2, 1.9, .8, 1.1), 3, 4, byrow = TRUE)
TR <- matrix(rep(c(1, 1, 3, 4), 3), 3, 4, byrow = TRUE)
for (method in c("ROC", "RS", "RR", "RE", "EW")) {
  f <- dea.order.ties(Y = T, OY = TR, SURROGATE = method, Weak = FALSE)
  ok(f)
  stopifnot(max(abs(rowSums(f$surrogate_output) - 1)) < tol,
            max(abs(f$virtual_output - f$surrogate_output)) < tol,
            max(abs(f$virtual_output[, 1] - f$virtual_output[, 2])) < tol)
}
f <- dea.order.ties(Y = T, OY = TR, Weak = FALSE)
stopifnot(max(abs(f$surrogate_output[1, ] - c(19/48,19/48,7/48,1/16))) < tol)
stopifnot(inherits(try(dea.order(Y = T, OY = TR), silent = TRUE), "try-error"),
          inherits(try(dea.order.ties(Y = T, OY = TR, TIES = "worst"), silent = TRUE), "try-error"))
# Directional example: strict ROC two criteria, data direction.
G <- matrix(c(5, 7, 6), ncol = 1); B <- matrix(c(4, 3, 5), ncol = 1)
RG <- matrix(1, 3, 1); RB <- matrix(2, 3, 1)
dw <- dir.dea.order(Y = G, B = B, OY = RG, OB = RB)
dm <- dir.dea.order(Y = G, B = B, OY = RG, OB = RB, Gamma = .5)
ds <- dir.dea.order(Y = G, B = B, OY = RG, OB = RB, Weak = FALSE)
for (f in list(dw, dm, ds)) ok(f)
stopifnot(all(dw$score + tol >= dm$score), all(dm$score + tol >= ds$score),
          max(abs(ds$virtual_share - ds$surrogate_share)) < tol,
          max(abs(ds$score - 1/(1 + ds$beta))) < tol)
# With both criteria tied, each directional share is one half.
dt <- dir.dea.order.ties(Y = G, B = B, OY = RG, OB = RG, Weak = FALSE)
ok(dt)
stopifnot(max(abs(dt$virtual_share - .5)) < tol)
# User-supplied directions and missing side.
dg <- dir.dea.order(Y = G, GY = 1)
ok(dg)
stopifnot(max(abs(dg$beta - c(2, 0, 1))) < tol)
# A zero-total radial virtual output is undefined rather than a false simplex.
f0 <- dea.order(Y = matrix(c(0, 1), ncol = 1))
ok(f0)
stopifnot(is.na(f0$virtual_output[1, 1]))
stopifnot(max(abs(ds$beta - c(.3625, 0, .225))) < tol)
# Tie-free implementations must agree with pooled implementation without ties.
tf <- dea.order.ties(Y = Y, OY = R, Gamma = .5)
ok(tf)
stopifnot(max(abs(tf$efficiency - m$efficiency)) < tol)
dtf <- dir.dea.order.ties(Y = G, B = B, OY = RG, OB = RB, Gamma = .5)
ok(dtf)
stopifnot(max(abs(dtf$score - dm$score)) < tol)
