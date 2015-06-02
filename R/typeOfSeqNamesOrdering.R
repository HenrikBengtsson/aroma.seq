setMethodS3("typeOfSequenceOrdering", "character", function(values, what=c("lexicographic", "canonical", "mixeddecimal", "mixedroman"), unique=FALSE, rank=TRUE, as=c("scores", "humanreadable"), locale="C", ...) {
  # Argument 'as':
  as <- match.arg(as)

  # Argument 'what':
  what <- match.arg(what, several.ok=TRUE)
  what0 <- what

  ## Check unique values; faster if lots of duplicates
  uvalues <- unique(values)
  if (unique) {
    values <- uvalues
    map <- NULL
  } else {
    map <- match(values, table=uvalues)
  }

  ## Coerce to integers?
  names_int <- suppressWarnings(as.integer(uvalues))
  all_ints <- !anyNA(names_int)
  if (all_ints) {
    if (is.element("mixeddecimal", what)) what <- c(what, "canonical")
  }

  n <- length(values)
  counts <- list()

  withLocale({
    ## (a) Lexicographic ordering, e.g. 1,10,11,...,2,20,21,22,MT,X,Y?
    if (is.element("lexicographic", what)) {
      o <- order(uvalues)
      if (!unique) o <- o[map]
      delta <- diff(o)
      counts$lexicographic <- sum(delta == 1L, na.rm=TRUE)
    }

    ## (b) Canonical ordering, e.g. 1,2,...,10,11,...,21,22,X,Y,MT?
    ##     (assumes human or mouse)
    if (is.element("canonical", what)) {
      chrs <- names_int
      if (!all_ints) {
        chrMax <- max(c(0L, chrs), na.rm=TRUE)
        chrs[is.element(uvalues, "X")] <- chrMax+1L
        chrs[is.element(uvalues, "Y")] <- chrMax+2L
        chrs[is.element(uvalues, c("M", "MT"))] <- chrMax+3L
      }
      o <- order(chrs)
      if (!unique) o <- o[map]
      delta <- diff(o)
      counts$canonical <- sum(delta == 1L, na.rm=TRUE)
    }

    ## (c) Mixed decimal ordering, e.g. 1,2,...,10,11,...,21,22,MT,X,Y?
    if (is.element("mixeddecimal", what)) {
      if (all_ints) {
        ## Identical to "canonical"
	counts$mixedorder <- counts$canonical
      } else {
        ## To expensive to calculate?
        ## gtools::mixedorder() is quite slow for > 10e3 elements
        if (length(uvalues) < 10e3) {
          o <- mixedorder(uvalues, numeric.type="decimal")
          if (!unique) o <- o[map]
          delta <- diff(o)
          counts$mixedorder <- sum(delta == 1L, na.rm=TRUE)
        } else {
          counts$mixedorder <- NA_integer_
        }
      }
    }

    ## (d) Mixed roman-numeral ordering, e.g. I,II,III,IV,...,IX,X,...
    if (is.element("mixedroman", what)) {
      if (all_ints) {
        ## Identical to "canonical"
	counts$mixedorder <- counts$canonical
      } else {
        ## To expensive to calculate?
        ## gtools::mixedorder() is quite slow for > 10e3 elements
        if (length(uvalues) < 10e3) {
          o <- mixedorder(uvalues, numeric.type="roman", roman.case="both")
          if (!unique) o <- o[map]
          delta <- diff(o)
          counts$mixedorder <- sum(delta == 1L, na.rm=TRUE)
        } else {
          counts$mixedorder <- NA_integer_
        }
      }
    }
}, category="LC_COLLATE", locale=locale)

  counts <- counts[what0]
  counts <- unlist(counts)
  scores <- (counts + 1L) / n

  if (rank) {
    o <- order(scores, decreasing=TRUE)
    scores <- scores[o]
  }

  if (as == "humanreadable") {
    if (n == 0) {
      scores <- "<any order; empty set>"
    } else if (n == 1) {
      scores <- "<any order; a single item>"
    } else {
      scores <- sprintf("%g%% %s", round(100*scores, digits=1L), names(scores))
    }
  }

  scores
})


setMethodS3("typeOfSequenceOrdering", "numeric", function(values, rank=TRUE, as=c("scores", "humanreadable"), ...) {
  # Argument 'as':
  as <- match.arg(as)

  n <- length(values)
  counts <- list()

  ## (a) Ordered by chromosome length?
  delta <- diff(order(values))
  counts$length <- sum(delta == 1L, na.rm=TRUE)

  counts <- unlist(counts)
  scores <- (counts + 1L) / n

  names <- names(values)
  if (!is.null(names)) {
    scoresN <- typeOfSequenceOrdering(names, rank=FALSE, as="scores", ...)
    scores <- c(scoresN, scores)
  }

  if (rank) {
    o <- order(scores, decreasing=TRUE)
    scores <- scores[o]
  }

  if (as == "humanreadable") {
    if (n == 0) {
      scores <- "<any order; empty set>"
    } else if (n == 1) {
      scores <- "<any order; a single item>"
    } else {
      scores <- sprintf("%g%% %s", round(100*scores, digits=1L), names(scores))
    }
  }

  scores
})


############################################################################
# HISTORY:
# 2015-05-06
# o Added typeOfSequenceOrdering().
############################################################################
