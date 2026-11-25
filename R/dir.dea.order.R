# -----------------------------------------------------------------------------
# dir.dea.order.R
# Directional ordinal Benefit-of-the-Doubt model (NO tied priorities).
#
# Y  = desirable indicators (larger is better)
# B  = undesirable indicators (smaller is better)
# OY = ordinal priorities corresponding to Y
# OB = ordinal priorities corresponding to B
#
# The model follows the multiplier form of the directional BoD logic:
#   min beta_o = -u' y_o + p' b_o + omega
#   s.t.        -u' y_j + p' b_j + omega >= 0, for all j
#               u' g_y,o + p' g_b,o = 1
# with omega free under VRS and omega = 0 under CRS.
# The default direction is the data direction g_y,o = y_o, g_b,o = b_o.
# Hence directional virtual shares are u_r*y_ro and p_k*b_ko and sum to 1.
# The composite score is CI_o = 1 / (1 + beta_o).
#
# Weak = TRUE  : weak/benevolent ordering + least-priority surrogate floor.
# Weak = FALSE : strong/strict exact surrogate virtual shares.
#
# Use dir.dea.order.ties() when any row contains tied priorities.
# -----------------------------------------------------------------------------

dir.dea.order <- function(Y = NULL, B = NULL, OY = NULL, OB = NULL,
                          RTS = "vrs", SURROGATE = "ROC", Weak = TRUE,
                          ALPHA = 2, GY = NULL, GB = NULL, Gamma = NULL) {

  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("Package 'lpSolve' is required. Install it with install.packages('lpSolve').")
  }

  .base_surrogate <- function(s, method = "ROC", alpha = 2) {
    method <- toupper(method); r <- seq_len(s)
    if (method == "ROC") {
      w <- vapply(r, function(j) sum(1 / seq.int(j, s)) / s, numeric(1))
    } else if (method == "RS") {
      w <- 2 * (s + 1 - r) / (s * (s + 1))
    } else if (method == "RR") {
      w <- (1 / r) / sum(1 / r)
    } else if (method == "RE") {
      if (!is.finite(alpha) || alpha <= 0) stop("ALPHA must be positive for RE weights.")
      w <- exp(-alpha * log(r)) / sum(exp(-alpha * log(r)))
    } else if (method == "EW") {
      w <- rep(1 / s, s)
    } else stop("SURROGATE must be one of 'ROC', 'RS', 'RR', 'RE', or 'EW'.")
    as.numeric(w / sum(w))
  }

  .rank_info_no_ties <- function(R, method, alpha) {
    R <- as.matrix(R); n <- nrow(R); s <- ncol(R)
    if (!is.numeric(R) || any(!is.finite(R))) stop("Ranks must be finite numeric values.")
    if (all(R == 0)) return(list(active = FALSE, positions = matrix(0,n,s), tau = matrix(0,n,s)))
    if (any(!is.finite(R)) || any(R <= 0)) {
      stop("Rank matrices must be all zeros or contain positive finite ranks for every indicator.")
    }
    base <- .base_surrogate(s, method, alpha)
    pos <- tau <- matrix(NA_real_, n, s)
    for (o in seq_len(n)) {
      rr <- as.numeric(R[o, ])
      if (anyDuplicated(rr)) {
        stop(paste0("Tied ranks detected in row ", o, ". Use dir.dea.order.ties()."))
      }
      ord <- order(rr, seq_along(rr)); pp <- integer(s); pp[ord] <- seq_len(s)
      pos[o, ] <- pp; tau[o, ] <- base[pp]
    }
    list(active = TRUE, positions = pos, tau = tau)
  }

  .direction_row <- function(D, observed, o, label) {
    p <- ncol(observed)
    if (p == 0L) return(numeric(0))
    if (is.null(D)) return(as.numeric(observed[o, ]))
    if (is.vector(D) && !is.list(D)) {
      if (length(D) != p) stop(label, " must have length ", p, ".")
      z <- as.numeric(D)
    } else {
      D <- as.matrix(D)
      if (!all(dim(D) == dim(observed))) stop(label, " must have the same dimensions as its data matrix.")
      z <- as.numeric(D[o, ])
    }
    if (any(!is.finite(z)) || any(z < 0)) stop(label, " must be finite and nonnegative.")
    z
  }

  if (!is.logical(Weak) || length(Weak) != 1L || is.na(Weak)) {
    stop("Weak must be TRUE (weak/benevolent) or FALSE (strong/strict).")
  }
  if (is.null(Gamma)) Gamma <- if (Weak) 0 else 1
  if (!is.numeric(Gamma) || length(Gamma) != 1L ||
      !is.finite(Gamma) || Gamma < 0 || Gamma > 1) {
    stop("Gamma must be NULL or a finite number between 0 and 1.")
  }
  if (!Weak && Gamma != 1) stop("Weak=FALSE requires Gamma=1 (or NULL).")
  Weak <- Gamma < 1
  if (!is.character(SURROGATE) || length(SURROGATE) != 1L ||
      is.na(SURROGATE) || !toupper(SURROGATE) %in% c("ROC", "RS", "RR", "RE", "EW")) {
    stop("SURROGATE must be ROC, RS, RR, RE or EW.")
  }
  if (!is.numeric(ALPHA) || length(ALPHA) != 1L || !is.finite(ALPHA) || ALPHA <= 0) {
    stop("ALPHA must be a finite positive number.")
  }
  RTS <- tolower(RTS)
  if (!RTS %in% c("crs", "vrs")) stop("RTS must be 'crs' or 'vrs'.")
  if (is.null(Y) && is.null(B)) stop("At least one of Y (desirable) or B (undesirable) must be supplied.")

  if (!is.null(Y)) Y <- as.matrix(Y)
  if (!is.null(B)) B <- as.matrix(B)
  n <- if (!is.null(Y)) nrow(Y) else nrow(B)
  if (is.null(Y)) Y <- matrix(numeric(0), nrow = n, ncol = 0L)
  if (is.null(B)) B <- matrix(numeric(0), nrow = n, ncol = 0L)
  if (nrow(Y) != nrow(B)) stop("Y and B must have the same number of rows (DMUs).")
  storage.mode(Y) <- "double"; storage.mode(B) <- "double"
  if (any(!is.finite(Y)) || any(!is.finite(B))) stop("Y and B must contain only finite values.")
  if (any(Y < 0) || any(B < 0)) stop("Directional BoD requires nonnegative desirable and undesirable indicators.")

  ng <- ncol(Y); nb <- ncol(B); s <- ng + nb
  if (n < 1L || s < 1L) stop("At least one DMU and one indicator are required.")
  if (is.null(OY)) OY <- matrix(0, n, ng)
  if (is.null(OB)) OB <- matrix(0, n, nb)
  OY <- as.matrix(OY); OB <- as.matrix(OB)
  if (!all(dim(OY) == c(n, ng))) stop("OY must have the same dimensions as Y.")
  if (!all(dim(OB) == c(n, nb))) stop("OB must have the same dimensions as B.")

  R <- cbind(OY, OB)
  rr <- .rank_info_no_ties(R, SURROGATE, ALPHA)

  beta <- score <- rep(NA_real_, n)
  solver_status <- rep(NA_integer_, n)
  weights_desirable <- matrix(NA_real_, n, ng, dimnames = list(rownames(Y), colnames(Y)))
  weights_undesirable <- matrix(NA_real_, n, nb, dimnames = list(rownames(B), colnames(B)))
  virtual_share <- matrix(NA_real_, n, s)
  surrogate_share <- rr$tau
  for (o in seq_len(n)) {
    gy <- .direction_row(GY, Y, o, "GY")
    gb <- .direction_row(GB, B, o, "GB")
    if (sum(gy) + sum(gb) <= 0) stop("The directional vector must contain at least one positive component for every DMU.")

    # Variables: u (desirable), p (undesirable), and omega+ / omega- under VRS.
    if (RTS == "crs") {
      nv <- ng + nb; op <- om <- NA_integer_
    } else {
      nv <- ng + nb + 2L; op <- ng + nb + 1L; om <- ng + nb + 2L
    }

    obj <- numeric(nv)
    if (ng) obj[seq_len(ng)] <- -Y[o, ]
    if (nb) obj[ng + seq_len(nb)] <- B[o, ]
    if (RTS == "vrs") { obj[op] <- 1; obj[om] <- -1 }

    A <- list(); cdir <- character(); crhs <- numeric()
    add <- function(row, sense, rhs) {
      A[[length(A) + 1L]] <<- row
      cdir <<- c(cdir, sense)
      crhs <<- c(crhs, rhs)
    }

    # Common directional frontier.
    for (j in seq_len(n)) {
      row <- numeric(nv)
      if (ng) row[seq_len(ng)] <- -Y[j, ]
      if (nb) row[ng + seq_len(nb)] <- B[j, ]
      if (RTS == "vrs") { row[op] <- 1; row[om] <- -1 }
      add(row, ">=", 0)
    }

    # Directional normalization u'g_y + p'g_b = 1.
    row <- numeric(nv)
    if (ng) row[seq_len(ng)] <- gy
    if (nb) row[ng + seq_len(nb)] <- gb
    add(row, "=", 1)

    if (rr$active) {
      pos <- rr$positions[o, ]; tau <- rr$tau[o, ]; ord <- order(pos)
      contrib <- vector("list", s)
      if (ng) for (r in seq_len(ng)) {
        z <- numeric(nv); z[r] <- gy[r]; contrib[[r]] <- z
      }
      if (nb) for (k in seq_len(nb)) {
        z <- numeric(nv); z[ng + k] <- gb[k]; contrib[[ng + k]] <- z
      }

      if (Weak) {
        if (s >= 2L) for (h in seq_len(s - 1L)) {
          a <- ord[h]; b <- ord[h + 1L]
          add(contrib[[a]] - contrib[[b]], ">=", Gamma * (tau[a] - tau[b]))
        }
        least <- ord[s]
        add(contrib[[least]], ">=", tau[least])
      } else {
        # Normalization already forces the shares to sum to one; s-1 exact
        # share equalities imply the final one automatically.
        if (s >= 2L) for (r in seq_len(s - 1L)) add(contrib[[r]], "=", tau[r])
      }
    }

    fit <- lpSolve::lp(direction = "min", objective.in = obj,
                       const.mat = do.call(rbind, A),
                       const.dir = cdir, const.rhs = crhs)
    solver_status[o] <- fit$status

    if (fit$status == 0L) {
      beta[o] <- fit$objval
      # Numerical tolerance: beta should be >= 0 because the evaluated DMU is
      # itself one of the frontier constraints.
      if (beta[o] < 0 && beta[o] > -1e-8) beta[o] <- 0
      score[o] <- if (beta[o] >= 0) 1 / (1 + beta[o]) else NA_real_
      if (ng) weights_desirable[o, ] <- fit$solution[seq_len(ng)]
      if (nb) weights_undesirable[o, ] <- fit$solution[ng + seq_len(nb)]
      sh <- numeric(s)
      if (ng) sh[seq_len(ng)] <- gy * fit$solution[seq_len(ng)]
      if (nb) sh[ng + seq_len(nb)] <- gb * fit$solution[ng + seq_len(nb)]
      virtual_share[o, ] <- sh
    }
  }

  # Safe column names after computation.
  yn <- if (ng) { z <- colnames(Y); if (is.null(z)) paste0("Y",seq_len(ng)) else z } else character(0)
  bn <- if (nb) { z <- colnames(B); if (is.null(z)) paste0("B",seq_len(nb)) else z } else character(0)
  rownames(virtual_share) <- rownames(surrogate_share) <- if (ng) rownames(Y) else rownames(B)
  colnames(virtual_share) <- c(yn, bn)
  colnames(surrogate_share) <- c(yn, bn)

  if (any(solver_status != 0L)) {
    warning(sum(solver_status != 0L), " LP(s) did not return an optimal solution; inspect solver_status.")
  }

  list(
    efficiency = score,
    score = score,
    beta = beta,
    weights_desirable = weights_desirable,
    weights_undesirable = weights_undesirable,
    virtual_share = virtual_share,
    surrogate_share = surrogate_share,
    ordinal_positions = rr$positions,
    solver_status = solver_status,
    preference = if (Gamma == 0) "weak" else if (Gamma == 1) "strict" else "moderate",
    Gamma = Gamma,
    Weak = Weak,
    RTS = RTS,
    SURROGATE = toupper(SURROGATE),
    ALPHA = ALPHA,
    direction = if (is.null(GY) && is.null(GB)) "data direction" else "user supplied"
  )
}
