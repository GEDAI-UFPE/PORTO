\name{dir.dea.order.ties}
\alias{dir.dea.order.ties}
\title{Directional Ordinal Benefit-of-the-Doubt with Ties}
\description{Estimates directional composite indicators with desirable and undesirable dimensions and pooled positional treatment of tied priorities.}
\usage{
dir.dea.order.ties(Y = NULL, B = NULL, OY = NULL, OB = NULL,
                               RTS = "vrs", SURROGATE = "ROC", Weak = TRUE,
                               ALPHA = 2, TIES = "pooled",
                               GY = NULL, GB = NULL, Gamma = NULL)
}
\arguments{
  \item{Y}{Numeric n by s_g matrix of desirable indicators, larger is better. NULL permits an all-undesirable assessment. Nonnegative finite values, DMUs in rows.}
  \item{B}{Numeric n by s_b matrix of undesirable indicators, smaller is better. NULL permits an all-desirable assessment. At least one of Y and B must be supplied; row counts must agree.}
  \item{OY}{Ordinal ranks for Y, aligned by row and column. Smaller labels represent higher priority. When both Y and B exist, OY and OB must jointly encode the full ordering across both groups, not independent rank sequences.}
  \item{OB}{Ordinal ranks for B, aligned by row and column. Supply positive finite ranks for every included indicator when restrictions are active. Both rank matrices NULL or jointly all zero disable restrictions; mixed zero and positive ranks are invalid.}
  \item{RTS}{Returns to scale: "crs" or "vrs" (case-insensitive). Radial models default to CRS; directional models default to VRS. VRS uses a free intercept, represented as the difference of two nonnegative variables.}
  \item{SURROGATE}{Surrogate method: "ROC" (rank-order centroid, default), "RS" (rank sum), "RR" (rank reciprocal), "RE" (inverse-power rank exponent), or "EW" (equal weights). Methods are case-insensitive.}
  \item{Weak}{Logical scalar. TRUE (default) gives weak/benevolent preferences; FALSE imposes strong/strict exact surrogate shares. With TRUE, Gamma may specify the moderate family.}
  \item{ALPHA}{Positive finite numeric scalar; exponent for RE, default 2. RE uses r^(-ALPHA) normalized across the rank positions. Other schemes do not use the exponent.}
  \item{TIES}{Character scalar, "pooled" only. Tied indicators share the total surrogate mass of the occupied positions equally. Legacy best, worst, average and sequential substitutions are not supported.}
  \item{GY}{Nonnegative direction magnitudes for desirable indicators: vector of length ncol(Y) shared by all DMUs, or matrix of the same dimensions as Y. NULL uses each evaluated DMU's observed desirable values.}
  \item{GB}{Nonnegative direction magnitudes for undesirable indicators: vector of length ncol(B), or matrix of the same dimensions as B. NULL uses each evaluated DMU's observed undesirable values. Supply positive magnitudes for contraction, not negative numbers.}
  \item{Gamma}{NULL (default) selects 0 when Weak=TRUE and 1 when Weak=FALSE. Supply a value in [0,1] with Weak=TRUE for the nested preference continuum. 0 is weak, 1 is strict, and intermediate values are moderate. Weak=FALSE requires Gamma=NULL or Gamma=1.}
}
\details{
For desirable values y_j and undesirable values b_j, solve a directional multiplier minimization with nonnegative multipliers u and p. The objective for unit o is -u'y_o + p'b_o + omega. The common supporting inequalities are -u'y_j + p'b_j + omega >= 0 for every j. Normalize u'g_y + p'g_b = 1, with nonnegative direction magnitudes g_y and g_b; the actual direction expands y and contracts b. Under VRS omega is free; under CRS omega=0. The default data direction uses g_y=y_o and g_b=b_o. At least one direction component must be positive for every DMU. Custom magnitudes change the distance scale and thus the numerical score.

Directional shares are u_r g_y,r and p_k g_b,k, whose combined sum is one by normalization. They encode the magnitude of weighted directional improvement, not signed net achievement. Apply the weak, moderate or strict rules to that combined vector. When using the data direction these are the weighted observed values. No reversal of undesirable indicators is performed. The score is 1/(1+beta), where beta is the minimized objective. This reciprocal scaling is PORTO's convention and differs from other DDF-to-index transformations; do not mix scales.

The nonnegative bad-indicator multipliers represent a monotone (free-disposal) supporting technology for the undesirable dimensions. This is not a weak-disposability equality model with unrestricted bad-output multipliers. The VRS model uses a convex reference set. Without the VRS intercept the CRS directional model uses a conical technology and may be infeasible; VRS is generally preferred for the composite-indicator setting. With ordinal restrictions, the distance is a restricted supporting gap, not necessarily a physically attainable improvement in the unmodified envelopment technology.

Ranks are supplied separately from observed performance: a low ordinal number means high priority regardless of whether an indicator is desirable or undesirable. The solver constructs one LP per evaluated DMU, with all supplied rows as reference observations. Rank labels need not be consecutive integers, but must be positive and finite on active rank matrices. Data and rank matrices are matched by position, not by name.

For an ordered set of s criteria the surrogate vector is normalized to sum to one. ROC uses the sum of reciprocal positions from rank r to s divided by s; RS uses 2(s+1-r)/(s(s+1)); RR uses normalized reciprocal ranks; RE uses normalized inverse-power ranks; EW uses 1/s. The numerical intensities are modeling judgments, not elicited cardinal preferences. EW makes even the weak model equal-share when all criteria on that side are ranked. Compare surrogate profiles with stakeholders and conduct sensitivity analysis.

For an active set with positive contribution total, weak preferences require nonincreasing virtual shares and the least-ranked share to be at least its surrogate value. Moderate preferences additionally require the difference between adjacent shares to be at least Gamma times the corresponding difference in surrogate shares. The least-share floor is unchanged. Strong preferences fix each share to its surrogate value. The latter is equivalent to requiring every share to be at least its surrogate value, since both vectors sum to one. The feasible multiplier sets are nested as Gamma increases, conditional on feasibility and well-defined virtual shares. Under a fixed technology, higher-is-better scores cannot increase with stronger restrictions.

Ties represent indifference, not missing information. If a tied block occupies consecutive positions a through b in the sorted list, each member receives the mean of all base surrogate weights at those positions. For four ROC ranks with the top two tied, the exact vector is (19/48, 19/48, 7/48, 1/16). This preserves the unit sum and ordering. Weak and moderate models impose equality within tied blocks and compare consecutive blocks. Strict models fix each member to its pooled share. The legacy worst-position rule is not equivalent.
}
\value{
A list with the following elements:
  \item{efficiency, score}{Identical numeric vectors of higher-is-better composite scores, computed as 1/(1+beta). NA for nonoptimal solves or invalid negative distances.}
  \item{beta}{Optimal directional distance (or restricted directional support gap). Smaller is better; zero indicates a supporting frontier. Tiny negative round-off is truncated to zero.}
  \item{weights_desirable, weights_undesirable}{Matrices of nonnegative multipliers for Y and B, one row per DMU.}
  \item{virtual_share}{Combined directional contribution shares, in column order cbind(Y,B). With data direction these are multiplier times observed indicator. With custom directions these are multiplier times direction magnitude. Each valid row sums to one.}
  \item{surrogate_share}{Combined rank-aligned surrogate shares, in the same column order. All zeros if restrictions are disabled.}
  \item{ordinal_positions}{Converted positional ranks; tied blocks use the average occupied position.}
  \item{solver_status}{Integer lpSolve status per DMU (0 is optimal). Nonzero status yields NA estimates and a warning.}
  \item{preference, Weak, Gamma}{Effective preference name, flag and intensity.}
  \item{RTS, SURROGATE, ALPHA}{Model settings used.}
  \item{direction}{Label indicating data direction or user-supplied directions.}
  \item{tie_blocks}{List(s) of integer column indices in each priority block, ordered from highest to lowest, one list per DMU. On inactive rank sets the entries are empty.}
  \item{tie_method}{Description of the pooled positional mass rule.}
}
\note{A zero direction component cannot receive a positive share. Weak, moderate or strict constraints can therefore be infeasible when the default data direction contains zeros. Choose substantively defensible positive direction magnitudes where appropriate; no epsilon or automatic imputation is inserted. Nonnegative data and finite ranks are required. Directional restrictions are DMU-specific; comparative scores describe conditional support, not causal effects of policies. Directional targets are not returned, and the restricted gap need not correspond to a feasible physical target. Rescale values and directions together if numerical conditioning is poor.}
\author{Thyago Nepomuceno, Department of Statistics, Federal University of Pernambuco. Email: nepomuceno@de.ufpe.be.}
\references{
Cherchye, L., Moesen, W., Rogge, N. and Van Puyenbroeck, T. (2007). An introduction to Benefit of the Doubt composite indicators. Social Indicators Research, 82, 111--145. DOI: 10.1007/s11205-006-9029-7.

Barron, F. H. and Barrett, B. E. (1996). Decision quality using ranked attribute weights. Management Science, 42(11), 1515--1523. DOI: 10.1287/mnsc.42.11.1515.

Zanella, A., Camanho, A. S. and Dias, T. G. (2015). The assessment of cities' livability integrating human wellbeing and environmental impact. Annals of Operations Research, 226, 695--726. DOI: 10.1007/s10479-014-1666-7.
}
\seealso{\code{\link{dea.order}}, \code{\link{dea.order.ties}}, \code{\link{dir.dea.order}}}
\examples{
Y <- matrix(c(5, 7, 6, 8, 5, 7), ncol = 2)
B <- matrix(c(4, 3, 5), ncol = 1)
# First desirable and undesirable indicator tied for top rank
OY <- matrix(rep(c(1, 3), 3), ncol = 2, byrow = TRUE)
OB <- matrix(1, 3, 1)
w <- dir.dea.order.ties(Y = Y, B = B, OY = OY, OB = OB)
s <- dir.dea.order.ties(Y = Y, B = B, OY = OY, OB = OB,
                        Weak = FALSE)
s$surrogate_share
s$virtual_share
rowSums(s$virtual_share)

# Preference sensitivity with a fixed surrogate rule
gammas <- c(0, 0.25, 0.5, 0.75, 1)
scores <- sapply(gammas, function(g) {
  dir.dea.order.ties(Y = Y, B = B, OY = OY, OB = OB,
                     Gamma = g)$score
})
scores
}
\keyword{models}
