\name{dea.order}
\alias{dea.order}
\title{Ordinal Data Envelopment Analysis and Benefit-of-the-Doubt without Ties}
\description{Estimates radial DEA multiplier scores or unit-input/unit-output BoD scores for rankings without ties.}
\usage{
dea.order(X = NULL, Y = NULL, OX = NULL, OY = NULL,
                      RTS = "crs", ORIENTATION = "out",
                      SURROGATE = "ROC", Weak = TRUE,
                      ALPHA = 2, AUTO_CRS_BOD = TRUE, Gamma = NULL)
}
\arguments{
  \item{X}{Numeric n by m input matrix or data frame (DMUs in rows). Nonnegative, finite values. NULL creates a single column of ones if Y is given.}
  \item{Y}{Numeric n by s desirable-output matrix or data frame. NULL creates a single column of ones if X is given. At least one of X or Y is required.}
  \item{OX}{Numeric ordinal rank matrix with the same dimensions and order as X. Smaller positive labels mean higher importance. NULL or an entire all-zero matrix disables input restrictions; partial zero ranks are invalid.}
  \item{OY}{Numeric ordinal rank matrix with the same dimensions and order as Y. Smaller positive labels mean higher importance. NULL or an entire all-zero matrix disables output restrictions; partial zero ranks are invalid.}
  \item{RTS}{Returns to scale: "crs" or "vrs" (case-insensitive). Radial models default to CRS; directional models default to VRS. VRS uses a free intercept, represented as the difference of two nonnegative variables.}
  \item{ORIENTATION}{"out" (default) maximizes weighted output with weighted input normalized to one; "in" minimizes weighted input with weighted output normalized to one and returns the reciprocal objective. These are historical interface labels, not the conventional names of radial DEA orientations: the first is the input-normalized multiplier form, the second the output-normalized form.}
  \item{SURROGATE}{Surrogate method: "ROC" (rank-order centroid, default), "RS" (rank sum), "RR" (rank reciprocal), "RE" (inverse-power rank exponent), or "EW" (equal weights). Methods are case-insensitive.}
  \item{Weak}{Logical scalar. TRUE (default) gives weak/benevolent preferences; FALSE imposes strong/strict exact surrogate shares. With TRUE, Gamma may specify the moderate family.}
  \item{ALPHA}{Positive finite numeric scalar; exponent for RE, default 2. RE uses r^(-ALPHA) normalized across the rank positions. Other schemes do not use the exponent.}
  \item{AUTO_CRS_BOD}{Logical scalar (default TRUE). If a unit input or output was automatically created and VRS requested, warn and use CRS. FALSE retains the requested VRS model; identical constant inputs or outputs can cause degeneracy.}
  \item{Gamma}{NULL (default) selects 0 when Weak=TRUE and 1 when Weak=FALSE. Supply a value in [0,1] with Weak=TRUE for the nested preference continuum. 0 is weak, 1 is strict, and intermediate values are moderate. Weak=FALSE requires Gamma=NULL or Gamma=1.}
}
\details{
Let x_j and y_j denote the input/output vectors of DMU j and v and u the corresponding nonnegative multipliers. The common supporting inequalities are u'y_j - v'x_j - omega <= 0 for all j. Under CRS omega=0; under VRS omega is free. ORIENTATION="out" maximizes u'y_o - omega subject to v'x_o=1. ORIENTATION="in" minimizes v'x_o + omega subject to u'y_o=1; the returned score is the reciprocal of that objective. These labels are retained for compatibility with the original routines. The CRS scores lie between zero and one for feasible nonnegative data with valid normalization.

For ordinary BoD supply Y only: X is a single unit input and the CRS constraints become u'y_j <= 1. For an all-undesirable radial assessment, supply X only and use ORIENTATION="in" (or the other CRS normalization, which gives the same efficiency under homogeneous restrictions). In conventional DEA, input and output shares are separately normalized; this does not impose a joint ranking across inputs and outputs. Use the directional functions for a joint ranking of good and bad indicators.

Ordinal share constraints are linearized by multiplying by the corresponding contribution total. For example, a moderate output restriction is u_a y_ao - u_b y_bo >= Gamma (tau_a-tau_b) sum_r u_r y_ro. A strict constraint is u_r y_ro = tau_r sum_q u_q y_qo. The VRS intercept is not part of the virtual-share denominator.

Ranks are supplied separately from observed performance: a low ordinal number means high priority regardless of whether an indicator is desirable or undesirable. The solver constructs one LP per evaluated DMU, with all supplied rows as reference observations. Rank labels need not be consecutive integers, but must be positive and finite on active rank matrices. Data and rank matrices are matched by position, not by name.

For an ordered set of s criteria the surrogate vector is normalized to sum to one. ROC uses the sum of reciprocal positions from rank r to s divided by s; RS uses 2(s+1-r)/(s(s+1)); RR uses normalized reciprocal ranks; RE uses normalized inverse-power ranks; EW uses 1/s. The numerical intensities are modeling judgments, not elicited cardinal preferences. EW makes even the weak model equal-share when all criteria on that side are ranked. Compare surrogate profiles with stakeholders and conduct sensitivity analysis.

For an active set with positive contribution total, weak preferences require nonincreasing virtual shares and the least-ranked share to be at least its surrogate value. Moderate preferences additionally require the difference between adjacent shares to be at least Gamma times the corresponding difference in surrogate shares. The least-share floor is unchanged. Strong preferences fix each share to its surrogate value. The latter is equivalent to requiring every share to be at least its surrogate value, since both vectors sum to one. The feasible multiplier sets are nested as Gamma increases, conditional on feasibility and well-defined virtual shares. Under a fixed technology, higher-is-better scores cannot increase with stronger restrictions.

Any repeated positive rank within a ranked row is rejected. Use the corresponding ties function for indifference. All-zero rank matrices are the documented exception and mean no restrictions.
}
\value{
A list with the following elements:
  \item{efficiency}{Numeric vector of higher-is-better scores, one per DMU. The score equals the objective for ORIENTATION="out" and its reciprocal for "in"; NA when the solver fails or the reciprocal is undefined.}
  \item{objective_raw}{Raw optimized LP objective, before the reciprocal transformation, where applicable.}
  \item{weights_input, weights_output}{Matrices of optimal nonnegative DEA multipliers, one row per DMU; an optimum need not be unique.}
  \item{virtual_input, virtual_output}{Matrices of normalized virtual contributions for each side separately, equal to the multiplier times the observed value divided by the corresponding total. A zero total is undefined and reported as NA.}
  \item{surrogate_input, surrogate_output}{Rank-aligned surrogate vectors for each side (all zeros on an unrestricted side). An active row sums to one.}
  \item{ordinal_positions_input, ordinal_positions_output}{Rank labels converted to occupied positions; for ties, the average occupied position is returned.}
  \item{potential_improvement}{Matrix of proportional output expansion amounts, (1/efficiency - 1) times Y, only for CRS and ORIENTATION="out" with positive scores. Unused columns and other cases are NA. With weight restrictions this is an algebraic indicator, not a guaranteed feasible physical projection.}
  \item{solver_status}{Integer status returned by lpSolve, 0 for optimal. Nonzero statuses have missing estimates; the function also issues a warning.}
  \item{preference, Weak, Gamma}{Effective preference name, logical flag (TRUE for weak or moderate), and intensity value.}
  \item{RTS, RTS_requested, ORIENTATION, SURROGATE, ALPHA}{Effective/requested returns to scale and other model settings.}
  \item{auto_unit_input, auto_unit_output}{Flags indicating whether a single artificial unit input/output was generated.}
}
\note{For ranked sides, strictly positive observed values are recommended. A zero observation can force zero virtual contributions; if the entire weighted total is zero, the linearized share restrictions are vacuous and the shares are undefined. Inspect the virtual shares and status before interpreting strict preferences. No artificial epsilon is inserted. The constant-input VRS radial model can return uninformative scores of one. Restrictions depend on the evaluated DMU, so the effective feasible multiplier sets differ across DMUs. Do not interpret a decrease in score as a causal policy effect. Multipliers and rankings may not be unique; inspect alternative optima if needed.}
\author{Thyago Nepomuceno, Department of Statistics, Federal University of Pernambuco. Email: nepomuceno@de.ufpe.be.}
\references{
Cherchye, L., Moesen, W., Rogge, N. and Van Puyenbroeck, T. (2007). An introduction to Benefit of the Doubt composite indicators. Social Indicators Research, 82, 111--145. DOI: 10.1007/s11205-006-9029-7.

Barron, F. H. and Barrett, B. E. (1996). Decision quality using ranked attribute weights. Management Science, 42(11), 1515--1523. DOI: 10.1287/mnsc.42.11.1515.

Zanella, A., Camanho, A. S. and Dias, T. G. (2015). The assessment of cities' livability integrating human wellbeing and environmental impact. Annals of Operations Research, 226, 695--726. DOI: 10.1007/s10479-014-1666-7.
}
\seealso{\code{\link{dea.order.ties}}, \code{\link{dir.dea.order}}, \code{\link{dir.dea.order.ties}}}
\examples{
Y <- matrix(c(5, 8, 8, 5, 6, 6), ncol = 2, byrow = TRUE)
OY <- matrix(rep(c(1, 2), 3), ncol = 2, byrow = TRUE)
base <- dea.order(Y = Y)
weak <- dea.order(Y = Y, OY = OY)
moderate <- dea.order(Y = Y, OY = OY, Gamma = 0.5)
strict <- dea.order(Y = Y, OY = OY, Weak = FALSE)
cbind(base$efficiency, weak$efficiency,
      moderate$efficiency, strict$efficiency)
strict$virtual_output
strict$solver_status

# A genuine input-output DEA example with input priorities
X <- matrix(c(2, 3, 3, 2, 4, 4), ncol = 2, byrow = TRUE)
OX <- matrix(rep(c(1, 2), 3), ncol = 2, byrow = TRUE)
fit <- dea.order(X = X, Y = Y, OX = OX, OY = OY,
                 RTS = "crs", Weak = FALSE)
fit$virtual_input

# All-undesirable radial model: one automatic unit output
bad <- dea.order(X = matrix(c(1, 2, 4), ncol = 1),
                 ORIENTATION = "in")
bad$efficiency
}
\keyword{models}
