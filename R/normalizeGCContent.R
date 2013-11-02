###########################################################################/**
# @RdocDefault normalizeGcContent
#
# @title "Normalizes signals for GC-content effects"
#
# \description{
#  @get "title". Some or all signals are used to estimated the
#  normalization function.  All signals are normalized.
# }
#
# @synopsis
#
# \arguments{
#   \item{y}{A @numeric @vector of length K of signals to be normalized.}
#   \item{gcContent}{A @numeric @vector of GC fractions.}
#   \item{targetFcn}{An optional @function.
#     If @NULL, the data is normalized to have constant GC content effects
#     (all equal to zero on the log-scale).}
##   \item{subsetToFit}{The subset of data points used to fit the
#      normalization function.
#      If @NULL, all data points are considered.}
#   \item{onMissing}{Specifies how data points for which there is no
#      fragment length is normalized.
#      If \code{"ignore"}, the values are not modified.
#      If \code{"median"}, the values are updated to have the same
#      robust average as the other data points.
#   }
#   \item{.isLogged}{A @logical.}
#   \item{...}{Additional arguments passed to @see "stats::lowess".}
#   \item{.returnFit}{A @logical.}
# }
#
# \value{
#   Returns a @numeric @vector of the normalized signals.
# }
#
# @author "HB"
#
# @keyword "nonparametric"
# @keyword "robust"
# @keyword internal
#*/###########################################################################
setMethodS3("normalizeGcContent", "default", function(y, gcContent, targetFcn=NULL, subsetToFit=NULL, onMissing=c("ignore", "median"), .isLogged=TRUE, ..., .returnFit=FALSE) {
  # predict() for 'lowess' is defined in aroma.light
  require("aroma.light") || throw("Package not loaded: aroma.light");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  y <- as.double(y);
  nbrOfDataPoints <- length(y);
  okY <- is.finite(y);
  # Sanity check
  if (!any(okY, na.rm=TRUE)) {
    throw("Cannot fit normalization function to GC content, because there are no (finite) data points in argument 'y'.");
  }

  # Argument 'gcContent':
  if (!is.vector(gcContent)) {
    throw("Argument 'gcContent' must be a vector: ", class(gcContent)[[1]]);
  }
  if (length(gcContent) != nbrOfDataPoints) {
    throw("Number of rows in argument 'gcContent' does not match the length of argument 'y': ", nrow(gcContent), " != ", nbrOfDataPoints);
  }

  # Assert that there are some finite GC contents
  hasGC <- is.finite(gcContent);
  if (!any(hasGC)) {
    throw("Cannot fit normalization function. Argument 'gcContent' contains no finite values.");
  }

  # Argument 'targetFcn':
  if (!is.null(targetFcn)) {
    if (!is.function(targetFcn)) {
      throw("Argument 'targetFcn' is not a function: ", class(targetFcn)[1]);
    }
  }

  # Argument 'subsetToFit':
  if (!is.null(subsetToFit)) {
    subsetToFit <- as.integer(subsetToFit);
    if (length(subsetToFit) > nbrOfDataPoints) {
      throw("The length of argument 'subsetToFit' does not match the number of data points: ", length(subsetToFit), " != ", nbrOfDataPoints);
    }
  }

  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate normalization function and predict the signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit smooth curve
  naValue <- as.double(NA);
  mu <- rep(naValue, times=nbrOfDataPoints);
  if (!is.null(targetFcn)) {
    muT <- rep(naValue, times=nbrOfDataPoints);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (a) Fit normalization function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ok <- okY & hasGC;
  if (!any(ok)) {
    throw("Cannot fit normalization function, because there exists no data points with finite values in both 'y' and 'gcContent'");
  }

  if (!is.null(subsetToFit)) {
    ok[-subsetToFit] <- FALSE;

    # Sanity check
    if (!any(ok)) {
      throw("Cannot fit normalization function, because after subsetting there exists no data points with finite values in both 'y' and 'gcContent'");
    }
  }

  # Fit finite {(gcContent, log2theta)_j} to data points
  suppressWarnings({
    fit <- lowess(gcContent[ok], y[ok], ...);
    class(fit) <- "lowess";
  })

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (b) Calculate correction factor
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate the correction factor for every data point
  mu[ok] <- predict(fit, newdata=gcContent[ok]);

  # Normalize toward a target function?
  if (!is.null(targetFcn)) {
    muT[ok] <- targetFcn(gcContent[ok]);
  }

  rm(ok, hasGC);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate the *average* predicted signal across enzymes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sum on the non-log scale.
  if (.isLogged) {
    mu <- 2^mu;
    if (!is.null(targetFcn))
      muT <- 2^muT;
  }

  # Special case: Units with unknown fragment lengths
  if (onMissing != "ignore") {
    isMissing <- !hasGC;
    if (any(isMissing)) {
      if (onMissing == "median") {
        # Let the predicted value for these units be the robust average
        # of all other units (based on the assumption that the missing
        # GC content are distributed as the known ones).

        # Identify the set to be used to estimate the target average
        ok <- (okY & !isMissing);
        # Sanity check
        if (!any(ok)) {
          throw("Cannot fit normalization function to loci with unknown GC content, because there are no (finite) data points to be fitted.");
        }

        if (!is.null(subsetToFit)) {
          ok[-subsetToFit] <- FALSE;
          # Sanity check
          if (!any(ok)) {
            throw("Cannot fit normalization function to loci with unknown GC content, because after subsetting there are no (finite) data points to be fitted.");
          }
        }

        # Substitute the predicted means with the median of the already
        # predicted set of loci.
        mu[isMissing] <- median(mu[ok], na.rm=TRUE);
        if (!is.null(targetFcn)) {
          muT[isMissing] <- median(muT[ok], na.rm=TRUE);
        }
        rm(ok);
      } # if (onMissing == "median")
    }
    rm(isMissing);
  }

  if (.isLogged) {
    mu <- log2(mu);
    if (!is.null(targetFcn))
      muT <- log2(muT);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate the correction ("normalization") factor
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate correction factors
  if (is.null(targetFcn)) {
    dy <- mu;
  } else {
    dy <- (mu - muT);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Transform signals
  ok <- is.finite(dy) & okY;
  rm(okY);
  y[ok] <- y[ok] - dy[ok];

  if (.returnFit) {
    fit <- list(fit=fit, mu=mu, targetFcn=targetFcn);
    if (!is.null(targetFcn)) fit$muT <- muT;
    attr(y, "modelFit") <- fit;
  }

  y;
}, private=TRUE)


############################################################################
# HISTORY:
# 2013-11-02
# o BUG FIX: normalizeGcContent() assumed that 'aroma.light' was attached.
# 2012-10-16
# o Created from normalizeFragmentLength.R of aroma.light.
############################################################################
