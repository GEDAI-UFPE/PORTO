# -----------------------------------------------------------------------------
# dea.order.R
# Corrected ordinal DEA / Benefit-of-the-Doubt multiplier model (no ties).
#
# Main corrections relative to the former routine:
#   1) If X = NULL and Y is supplied, a SINGLE unit input is created (n x 1),
#      as required by the standard BoD formulation.
#   2) Weak = TRUE (default) implements benevolent/weak ordinal preferences:
#        virtual contribution(rank 1) >= ... >= virtual contribution(rank s)
#      plus a surrogate lower bound for the least-preferred indicator.
#   3) Weak = FALSE implements the strong/strict formulation:
#      each virtual share is fixed to its surrogate share.
#   4) VRS uses a truly free intercept represented as omega+ - omega-, because
#      lpSolve variables are nonnegative by default.
#   5) This function is for rankings WITHOUT ties. Use dea.order.ties() when
#      a row of OX or OY contains tied priorities.
#
# Gamma: optional moderate preference intensity (0 to 1).
# Surrogate schemes: ROC, RS, RR, RE, EW.
# RE uses ALPHA = 2 by default.
# -----------------------------------------------------------------------------

dea.order <- function(X = NULL, Y = NULL, OX = NULL, OY = NULL,
                      RTS = "crs", ORIENTATION = "out",
                      SURROGATE = "ROC", Weak = TRUE,
                      ALPHA = 2, AUTO_CRS_BOD = TRUE, Gamma = NULL) {

  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("Package 'lpSolve' is required. Install it with install.packages('lpSolve').")
  }

  # ----- Local helpers --------------------------------------------------------
  .base_surrogate <- function(s, method = "ROC", alpha = 2) {
    method <- toupper(method)
    r <- seq_len(s)
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
    } else {
      stop("SURROGATE must be one of 'ROC', 'RS', 'RR', 'RE', or 'EW'.")
    }
    as.numeric(w / sum(w))
  }

  .rank_to_surrogate_no_ties <- function(rank_matrix, method, alpha) {
    rank_matrix <- as.matrix(rank_matrix)
    n <- nrow(rank_matrix)
    s <- ncol(rank_matrix)

    if (!is.numeric(rank_matrix) || any(!is.finite(rank_matrix))) stop("Ranks must be finite numeric values.")
    if (all(rank_matrix == 0)) {
      return(list(active = FALSE,
                  positions = matrix(0, n, s),
                  tau = matrix(0, n, s)))
    }
    if (any(!is.finite(rank_matrix)) || any(rank_matrix <= 0)) {
      stop("A rank matrix must be either all zeros (no restrictions) or contain positive finite ranks.")
    }

    pos <- matrix(NA_real_, n, s)
    tau <- matrix(NA_real_, n, s)
    base <- .base_surrogate(s, method, alpha)

    for (o in seq_len(n)) {
      rr <- as.numeric(rank_matrix[o, ])
      if (anyDuplicated(rr)) {
        stop(paste0("Tied ranks detected in row ", o,
                    ". Use dea.order.ties() for tied ordinal priorities."))
      }
      # Convert arbitrary ordinal labels to positions 1,...,s.
      ord <- order(rr, seq_along(rr))
      pp <- integer(s)
      pp[ord] <- seq_len(s)
      pos[o, ] <- pp
      tau[o, ] <- base[pp]
    }

    list(active = TRUE, positions = pos, tau = tau)
  }

  .append <- function(A, dir, rhs, row, sense, b) {
    A[[length(A) + 1L]] <- row
    dir <- c(dir, sense)
    rhs <- c(rhs, b)
    list(A = A, dir = dir, rhs = rhs)
  }

  # ----- Validation / default BoD setup --------------------------------------
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
  ORIENTATION <- tolower(ORIENTATION)
  if (!ORIENTATION %in% c("out", "in")) stop("ORIENTATION must be 'out' or 'in'.")

  if (is.null(X) && is.null(Y)) stop("At least one of X or Y must be supplied.")

  auto_unit_input <- FALSE
  auto_unit_output <- FALSE

  if (is.null(X)) {
    Y <- as.matrix(Y)
    if (nrow(Y) < 1L || ncol(Y) < 1L) stop("Y must be a non-empty matrix.")
    X <- matrix(1, nrow = nrow(Y), ncol = 1L,
                dimnames = list(rownames(Y), "UnitInput"))
    auto_unit_input <- TRUE
  }
  if (is.null(Y)) {
    X <- as.matrix(X)
    if (nrow(X) < 1L || ncol(X) < 1L) stop("X must be a non-empty matrix.")
    Y <- matrix(1, nrow = nrow(X), ncol = 1L,
                dimnames = list(rownames(X), "UnitOutput"))
    auto_unit_output <- TRUE
  }

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  if (nrow(X) < 1L || ncol(X) < 1L || nrow(Y) < 1L || ncol(Y) < 1L) {
    stop("X and Y must be non-empty matrices.")
  }
  if (!is.logical(AUTO_CRS_BOD) || length(AUTO_CRS_BOD) != 1L || is.na(AUTO_CRS_BOD)) {
    stop("AUTO_CRS_BOD must be TRUE or FALSE.")
  }
  storage.mode(X) <- "double"
  storage.mode(Y) <- "double"

  if (nrow(X) != nrow(Y)) stop("X and Y must have the same number of rows (DMUs).")
  if (any(!is.finite(X)) || any(!is.finite(Y))) stop("X and Y must contain only finite values.")
  if (any(X < 0) || any(Y < 0)) stop("This multiplier formulation requires nonnegative X and Y.")

  n <- nrow(X)
  m <- ncol(X)
  s <- ncol(Y)

  if (is.null(OX)) OX <- matrix(0, nrow = n, ncol = m)
  if (is.null(OY)) OY <- matrix(0, nrow = n, ncol = s)
  OX <- as.matrix(OX)
  OY <- as.matrix(OY)
  if (!all(dim(OX) == c(n, m))) stop("OX must have the same dimensions as X.")
  if (!all(dim(OY) == c(n, s))) stop("OY must have the same dimensions as Y.")

  # A VRS radial multiplier model with an identical unit input is degenerate for
  # standard output-maximizing BoD. The usual BoD model is CRS.
  RTS_requested <- RTS
  if ((auto_unit_input || auto_unit_output) && RTS == "vrs" && AUTO_CRS_BOD) {
    warning("An automatic unit input or output can make radial VRS BoD degenerate; RTS has been set to 'crs'.")
    RTS <- "crs"
  }

  rx <- .rank_to_surrogate_no_ties(OX, SURROGATE, ALPHA)
  ry <- .rank_to_surrogate_no_ties(OY, SURROGATE, ALPHA)

  # ----- Storage --------------------------------------------------------------
  efficiency <- rep(NA_real_, n)
  objective_raw <- rep(NA_real_, n)
  solver_status <- rep(NA_integer_, n)
  weights_input <- matrix(NA_real_, n, m, dimnames = list(rownames(X), colnames(X)))
  weights_output <- matrix(NA_real_, n, s, dimnames = list(rownames(Y), colnames(Y)))
  virtual_input <- matrix(NA_real_, n, m, dimnames = list(rownames(X), colnames(X)))
  virtual_output <- matrix(NA_real_, n, s, dimnames = list(rownames(Y), colnames(Y)))
  potential_improvement <- matrix(NA_real_, n, max(m, s))

  # ----- DEA loop -------------------------------------------------------------
  for (o in seq_len(n)) {

    # Variables are [v_1,...,v_m, u_1,...,u_s] under CRS and
    # [v, u, omega_plus, omega_minus] under VRS, where
    # omega = omega_plus - omega_minus is free.
    if (RTS == "crs") {
      nv <- m + s
      op <- om <- NA_integer_
    } else {
      nv <- m + s + 2L
      op <- m + s + 1L
      om <- m + s + 2L
    }

    A <- list()
    cdir <- character()
    crhs <- numeric()

    add_constraint <- function(row, sense, b) {
      tmp <- .append(A, cdir, crhs, row, sense, b)
      A <<- tmp$A
      cdir <<- tmp$dir
      crhs <<- tmp$rhs
    }

    if (ORIENTATION == "out") {
      # Input-normalized, output-maximizing multiplier / BoD formulation:
      # max u'y_o - omega, s.t. u'y_j - v'x_j - omega <= 0,
      # v'x_o = 1. For CRS omega = 0.
      obj <- numeric(nv)
      obj[(m + 1L):(m + s)] <- Y[o, ]
      if (RTS == "vrs") {
        obj[op] <- -1
        obj[om] <-  1
      }
      direction <- "max"

      for (j in seq_len(n)) {
        row <- numeric(nv)
        row[seq_len(m)] <- -X[j, ]
        row[m + seq_len(s)] <- Y[j, ]
        if (RTS == "vrs") {
          row[op] <- -1
          row[om] <-  1
        }
        add_constraint(row, "<=", 0)
      }
      row <- numeric(nv)
      row[seq_len(m)] <- X[o, ]
      add_constraint(row, "=", 1)

    } else {
      # Output-oriented BCC/CCR dual form:
      # min v'x_o + omega, s.t. u'y_j - v'x_j - omega <= 0,
      # u'y_o = 1. Returned efficiency is 1 / objective.
      obj <- numeric(nv)
      obj[seq_len(m)] <- X[o, ]
      if (RTS == "vrs") {
        obj[op] <-  1
        obj[om] <- -1
      }
      direction <- "min"

      for (j in seq_len(n)) {
        row <- numeric(nv)
        row[seq_len(m)] <- -X[j, ]
        row[m + seq_len(s)] <- Y[j, ]
        if (RTS == "vrs") {
          row[op] <- -1
          row[om] <-  1
        }
        add_constraint(row, "<=", 0)
      }
      row <- numeric(nv)
      row[m + seq_len(s)] <- Y[o, ]
      add_constraint(row, "=", 1)
    }

    # ----- Ordinal restrictions: inputs --------------------------------------
    if (rx$active) {
      pos <- rx$positions[o, ]
      tau <- rx$tau[o, ]
      order_idx <- order(pos)

      # Contribution coefficient vector for each input: x_io * v_i.
      input_contrib_row <- function(i) {
        rr <- numeric(nv)
        rr[i] <- X[o, i]
        rr
      }
      input_total_row <- numeric(nv)
      input_total_row[seq_len(m)] <- X[o, ]

      if (Weak) {
        if (m >= 2L) {
          for (h in seq_len(m - 1L)) {
            a <- order_idx[h]
            b <- order_idx[h + 1L]
            add_constraint(input_contrib_row(a) - input_contrib_row(b) -
                           Gamma * (tau[a] - tau[b]) * input_total_row, ">=", 0)
          }
        }
        least <- order_idx[m]
        add_constraint(input_contrib_row(least) - tau[least] * input_total_row, ">=", 0)
      } else {
        # Exact surrogate virtual shares. The final equality is redundant because
        # sum(tau)=1, so s-1 independent equalities are sufficient.
        if (m >= 2L) {
          for (i in seq_len(m - 1L)) {
            add_constraint(input_contrib_row(i) - tau[i] * input_total_row, "=", 0)
          }
        }
      }
    }

    # ----- Ordinal restrictions: outputs -------------------------------------
    if (ry$active) {
      pos <- ry$positions[o, ]
      tau <- ry$tau[o, ]
      order_idx <- order(pos)

      output_contrib_row <- function(r) {
        rr <- numeric(nv)
        rr[m + r] <- Y[o, r]
        rr
      }
      output_total_row <- numeric(nv)
      output_total_row[m + seq_len(s)] <- Y[o, ]

      if (Weak) {
        if (s >= 2L) {
          for (h in seq_len(s - 1L)) {
            a <- order_idx[h]
            b <- order_idx[h + 1L]
            add_constraint(output_contrib_row(a) - output_contrib_row(b) -
                           Gamma * (tau[a] - tau[b]) * output_total_row, ">=", 0)
          }
        }
        least <- order_idx[s]
        add_constraint(output_contrib_row(least) - tau[least] * output_total_row, ">=", 0)
      } else {
        if (s >= 2L) {
          for (r in seq_len(s - 1L)) {
            add_constraint(output_contrib_row(r) - tau[r] * output_total_row, "=", 0)
          }
        }
      }
    }

    const.mat <- do.call(rbind, A)
    fit <- lpSolve::lp(direction = direction,
                       objective.in = obj,
                       const.mat = const.mat,
                       const.dir = cdir,
                       const.rhs = crhs)

    solver_status[o] <- fit$status
    if (fit$status == 0L) {
      objective_raw[o] <- fit$objval
      weights_input[o, ] <- fit$solution[seq_len(m)]
      weights_output[o, ] <- fit$solution[m + seq_len(s)]

      if (ORIENTATION == "out") {
        efficiency[o] <- fit$objval
      } else {
        efficiency[o] <- if (fit$objval > 0) 1 / fit$objval else NA_real_
      }

      vin <- X[o, ] * weights_input[o, ]
      vout <- Y[o, ] * weights_output[o, ]
      virtual_input[o, ] <- if (sum(vin) > 0) vin / sum(vin) else NA_real_
      virtual_output[o, ] <- if (sum(vout) > 0) vout / sum(vout) else NA_real_

      # Only the CRS output-maximizing BoD score has the direct radial expansion
      # interpretation below. Other cases are intentionally left NA.
      if (ORIENTATION == "out" && RTS == "crs" && is.finite(efficiency[o]) && efficiency[o] > 0) {
        potential_improvement[o, seq_len(s)] <- (1 / efficiency[o] - 1) * Y[o, ]
      }
    }
  }

  if (any(solver_status != 0L)) {
    warning(sum(solver_status != 0L), " LP(s) did not return an optimal solution; inspect solver_status.")
  }

  list(
    weights_input = weights_input,
    weights_output = weights_output,
    efficiency = efficiency,
    objective_raw = objective_raw,
    potential_improvement = potential_improvement,
    virtual_input = virtual_input,
    virtual_output = virtual_output,
    surrogate_input = rx$tau,
    surrogate_output = ry$tau,
    ordinal_positions_input = rx$positions,
    ordinal_positions_output = ry$positions,
    solver_status = solver_status,
    preference = if (Gamma == 0) "weak" else if (Gamma == 1) "strict" else "moderate",
    Gamma = Gamma,
    Weak = Weak,
    RTS = RTS,
    RTS_requested = RTS_requested,
    ORIENTATION = ORIENTATION,
    SURROGATE = toupper(SURROGATE),
    ALPHA = ALPHA,
    auto_unit_input = auto_unit_input,
    auto_unit_output = auto_unit_output
  )
}
