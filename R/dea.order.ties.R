# -----------------------------------------------------------------------------
# dea.order.ties.R
# Corrected ordinal DEA / BoD multiplier model WITH tied priorities.
#
# Ties are handled by pooled positional surrogate mass. If t indicators occupy
# ordinal positions a,...,b, each tied indicator receives
#       mean(base_surrogate[a:b]).
# Example: ROC with four indicators and a two-way tie for first:
# (0.5208, 0.2708, 0.1458, 0.0625) ->
# (0.3958, 0.3958, 0.1458, 0.0625).
# This preserves indifference, ordering, and sum(tau)=1.
#
# Weak = TRUE  : weak/benevolent ordering + least-priority surrogate floor.
# Weak = FALSE : strong/strict exact surrogate virtual shares.
# -----------------------------------------------------------------------------

dea.order.ties <- function(X = NULL, Y = NULL, OX = NULL, OY = NULL,
                           RTS = "crs", ORIENTATION = "out",
                           SURROGATE = "ROC", Weak = TRUE,
                           ALPHA = 2, TIES = "pooled",
                           AUTO_CRS_BOD = TRUE, Gamma = NULL) {

  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("Package 'lpSolve' is required. Install it with install.packages('lpSolve').")
  }

  if (tolower(TIES) != "pooled") {
    stop("The corrected tie treatment is TIES='pooled'. Legacy best/worst/average/sequential rank substitution does not preserve the surrogate unit-sum structure.")
  }

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

  .tie_rank_info <- function(rank_matrix, method, alpha) {
    rank_matrix <- as.matrix(rank_matrix)
    n <- nrow(rank_matrix)
    s <- ncol(rank_matrix)

    if (!is.numeric(rank_matrix) || any(!is.finite(rank_matrix))) stop("Ranks must be finite numeric values.")
    if (all(rank_matrix == 0)) {
      return(list(active = FALSE,
                  positions = matrix(0, n, s),
                  tau = matrix(0, n, s),
                  blocks = vector("list", n)))
    }
    if (any(!is.finite(rank_matrix)) || any(rank_matrix <= 0)) {
      stop("A rank matrix must be either all zeros (no restrictions) or contain positive finite ranks.")
    }

    base <- .base_surrogate(s, method, alpha)
    pos <- matrix(NA_real_, n, s)
    tau <- matrix(NA_real_, n, s)
    blocks_all <- vector("list", n)

    for (o in seq_len(n)) {
      rr <- as.numeric(rank_matrix[o, ])
      levels <- sort(unique(rr))
      blocks <- lapply(levels, function(z) which(rr == z))
      blocks_all[[o]] <- blocks
      cursor <- 1L

      for (bl in blocks) {
        t <- length(bl)
        occupied <- cursor:(cursor + t - 1L)
        pooled <- sum(base[occupied]) / t
        tau[o, bl] <- pooled
        # All members of the tie block share the same ordinal block position.
        pos[o, bl] <- mean(occupied)
        cursor <- cursor + t
      }
      if (abs(sum(tau[o, ]) - 1) > 1e-10) stop("Internal tie-surrogate error: shares do not sum to one.")
    }

    list(active = TRUE, positions = pos, tau = tau, blocks = blocks_all)
  }

  .append <- function(A, dir, rhs, row, sense, b) {
    A[[length(A) + 1L]] <- row
    dir <- c(dir, sense)
    rhs <- c(rhs, b)
    list(A = A, dir = dir, rhs = rhs)
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
  ORIENTATION <- tolower(ORIENTATION)
  if (!ORIENTATION %in% c("out", "in")) stop("ORIENTATION must be 'out' or 'in'.")
  if (is.null(X) && is.null(Y)) stop("At least one of X or Y must be supplied.")

  auto_unit_input <- FALSE
  auto_unit_output <- FALSE
  if (is.null(X)) {
    Y <- as.matrix(Y)
    X <- matrix(1, nrow = nrow(Y), ncol = 1L,
                dimnames = list(rownames(Y), "UnitInput"))
    auto_unit_input <- TRUE
  }
  if (is.null(Y)) {
    X <- as.matrix(X)
    Y <- matrix(1, nrow = nrow(X), ncol = 1L,
                dimnames = list(rownames(X), "UnitOutput"))
    auto_unit_output <- TRUE
  }

  X <- as.matrix(X); Y <- as.matrix(Y)
  if (nrow(X) < 1L || ncol(X) < 1L || nrow(Y) < 1L || ncol(Y) < 1L) {
    stop("X and Y must be non-empty matrices.")
  }
  if (!is.logical(AUTO_CRS_BOD) || length(AUTO_CRS_BOD) != 1L || is.na(AUTO_CRS_BOD)) {
    stop("AUTO_CRS_BOD must be TRUE or FALSE.")
  }
  storage.mode(X) <- "double"; storage.mode(Y) <- "double"
  if (nrow(X) != nrow(Y)) stop("X and Y must have the same number of rows (DMUs).")
  if (any(!is.finite(X)) || any(!is.finite(Y))) stop("X and Y must contain only finite values.")
  if (any(X < 0) || any(Y < 0)) stop("This multiplier formulation requires nonnegative X and Y.")

  n <- nrow(X); m <- ncol(X); s <- ncol(Y)
  if (is.null(OX)) OX <- matrix(0, n, m)
  if (is.null(OY)) OY <- matrix(0, n, s)
  OX <- as.matrix(OX); OY <- as.matrix(OY)
  if (!all(dim(OX) == c(n, m))) stop("OX must have the same dimensions as X.")
  if (!all(dim(OY) == c(n, s))) stop("OY must have the same dimensions as Y.")

  RTS_requested <- RTS
  if ((auto_unit_input || auto_unit_output) && RTS == "vrs" && AUTO_CRS_BOD) {
    warning("An automatic unit input or output can make radial VRS BoD degenerate; RTS has been set to 'crs'.")
    RTS <- "crs"
  }

  rx <- .tie_rank_info(OX, SURROGATE, ALPHA)
  ry <- .tie_rank_info(OY, SURROGATE, ALPHA)

  efficiency <- objective_raw <- rep(NA_real_, n)
  solver_status <- rep(NA_integer_, n)
  weights_input <- matrix(NA_real_, n, m, dimnames = list(rownames(X), colnames(X)))
  weights_output <- matrix(NA_real_, n, s, dimnames = list(rownames(Y), colnames(Y)))
  virtual_input <- matrix(NA_real_, n, m, dimnames = list(rownames(X), colnames(X)))
  virtual_output <- matrix(NA_real_, n, s, dimnames = list(rownames(Y), colnames(Y)))
  potential_improvement <- matrix(NA_real_, n, max(m, s))

  for (o in seq_len(n)) {
    if (RTS == "crs") {
      nv <- m + s; op <- om <- NA_integer_
    } else {
      nv <- m + s + 2L; op <- m + s + 1L; om <- m + s + 2L
    }

    A <- list(); cdir <- character(); crhs <- numeric()
    add_constraint <- function(row, sense, b) {
      tmp <- .append(A, cdir, crhs, row, sense, b)
      A <<- tmp$A; cdir <<- tmp$dir; crhs <<- tmp$rhs
    }

    if (ORIENTATION == "out") {
      obj <- numeric(nv); obj[m + seq_len(s)] <- Y[o, ]
      if (RTS == "vrs") { obj[op] <- -1; obj[om] <- 1 }
      direction <- "max"
      for (j in seq_len(n)) {
        row <- numeric(nv); row[seq_len(m)] <- -X[j, ]; row[m + seq_len(s)] <- Y[j, ]
        if (RTS == "vrs") { row[op] <- -1; row[om] <- 1 }
        add_constraint(row, "<=", 0)
      }
      row <- numeric(nv); row[seq_len(m)] <- X[o, ]; add_constraint(row, "=", 1)
    } else {
      obj <- numeric(nv); obj[seq_len(m)] <- X[o, ]
      if (RTS == "vrs") { obj[op] <- 1; obj[om] <- -1 }
      direction <- "min"
      for (j in seq_len(n)) {
        row <- numeric(nv); row[seq_len(m)] <- -X[j, ]; row[m + seq_len(s)] <- Y[j, ]
        if (RTS == "vrs") { row[op] <- -1; row[om] <- 1 }
        add_constraint(row, "<=", 0)
      }
      row <- numeric(nv); row[m + seq_len(s)] <- Y[o, ]; add_constraint(row, "=", 1)
    }

    # ----- Input ordinal restrictions ----------------------------------------
    if (rx$active) {
      tau <- rx$tau[o, ]; blocks <- rx$blocks[[o]]
      contrib <- lapply(seq_len(m), function(i) {
        z <- numeric(nv); z[i] <- X[o, i]; z
      })
      total <- numeric(nv); total[seq_len(m)] <- X[o, ]

      if (Weak) {
        # Indifference within a tied block.
        for (bl in blocks) if (length(bl) > 1L) {
          for (k in bl[-1L]) add_constraint(contrib[[bl[1L]]] - contrib[[k]], "=", 0)
        }
        # Weak ordering between adjacent priority blocks.
        if (length(blocks) >= 2L) {
          for (h in seq_len(length(blocks) - 1L)) {
            a <- blocks[[h]][1L]; b <- blocks[[h + 1L]][1L]
            add_constraint(contrib[[a]] - contrib[[b]] -
                           Gamma * (tau[a] - tau[b]) * total, ">=", 0)
          }
        }
        # Surrogate lower bound only for the least-priority block.
        for (a in blocks[[length(blocks)]]) {
          add_constraint(contrib[[a]] - tau[a] * total, ">=", 0)
        }
      } else {
        # Exact surrogate shares. One equation is redundant because sum(tau)=1.
        if (m >= 2L) {
          for (i in seq_len(m - 1L)) add_constraint(contrib[[i]] - tau[i] * total, "=", 0)
        }
      }
    }

    # ----- Output ordinal restrictions ---------------------------------------
    if (ry$active) {
      tau <- ry$tau[o, ]; blocks <- ry$blocks[[o]]
      contrib <- lapply(seq_len(s), function(r) {
        z <- numeric(nv); z[m + r] <- Y[o, r]; z
      })
      total <- numeric(nv); total[m + seq_len(s)] <- Y[o, ]

      if (Weak) {
        for (bl in blocks) if (length(bl) > 1L) {
          for (k in bl[-1L]) add_constraint(contrib[[bl[1L]]] - contrib[[k]], "=", 0)
        }
        if (length(blocks) >= 2L) {
          for (h in seq_len(length(blocks) - 1L)) {
            a <- blocks[[h]][1L]; b <- blocks[[h + 1L]][1L]
            add_constraint(contrib[[a]] - contrib[[b]] -
                           Gamma * (tau[a] - tau[b]) * total, ">=", 0)
          }
        }
        for (a in blocks[[length(blocks)]]) {
          add_constraint(contrib[[a]] - tau[a] * total, ">=", 0)
        }
      } else {
        if (s >= 2L) {
          for (r in seq_len(s - 1L)) add_constraint(contrib[[r]] - tau[r] * total, "=", 0)
        }
      }
    }

    fit <- lpSolve::lp(direction = direction,
                       objective.in = obj,
                       const.mat = do.call(rbind, A),
                       const.dir = cdir,
                       const.rhs = crhs)
    solver_status[o] <- fit$status

    if (fit$status == 0L) {
      objective_raw[o] <- fit$objval
      weights_input[o, ] <- fit$solution[seq_len(m)]
      weights_output[o, ] <- fit$solution[m + seq_len(s)]
      efficiency[o] <- if (ORIENTATION == "out") fit$objval else if (fit$objval > 0) 1 / fit$objval else NA_real_

      vin <- X[o, ] * weights_input[o, ]
      vout <- Y[o, ] * weights_output[o, ]
      virtual_input[o, ] <- if (sum(vin) > 0) vin / sum(vin) else NA_real_
      virtual_output[o, ] <- if (sum(vout) > 0) vout / sum(vout) else NA_real_

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
    tie_blocks_input = rx$blocks,
    tie_blocks_output = ry$blocks,
    solver_status = solver_status,
    preference = if (Gamma == 0) "weak" else if (Gamma == 1) "strict" else "moderate",
    Gamma = Gamma,
    Weak = Weak,
    tie_method = "pooled positional surrogate mass",
    RTS = RTS,
    RTS_requested = RTS_requested,
    ORIENTATION = ORIENTATION,
    SURROGATE = toupper(SURROGATE),
    ALPHA = ALPHA,
    auto_unit_input = auto_unit_input,
    auto_unit_output = auto_unit_output
  )
}
