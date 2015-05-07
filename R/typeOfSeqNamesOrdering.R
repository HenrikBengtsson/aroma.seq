setMethodS3("typeOfSequenceOrdering", "character", function(values, sort=TRUE, as=c("scores", "humanreadable"), locale="C", ...) {
  # Argument 'as':
  as <- match.arg(as)

  n <- length(values)
  counts <- list()

  withLocale({
    ## (a) Lexicographic ordering, e.g. 1,10,11,...,2,20,21,22,MT,X,Y?
    delta <- diff(order(values))
    counts$lexicograpic <- sum(delta == 1L, na.rm=TRUE)
   ## (b) Mixed sort ordering, e.g. 1,2,...,10,11,...,21,22,MT,X,Y?
    delta <- diff(mixedorder(values))
    counts$mixedsort <- sum(delta == 1L, na.rm=TRUE)
   ## (d) Classical ordering, e.g. 1,2,...,10,11,...,21,22,X,Y,MT?
    ##     (assumes human or mouse)
    chrs <- suppressWarnings(as.integer(values))
    chrs[is.element(values, "X")] <- 23L
    chrs[is.element(values, "Y")] <- 24L
    chrs[is.element(values, c("M", "MT"))] <- 25L
    delta <- diff(order(chrs))
    counts$classical <- sum(delta == 1L, na.rm=TRUE)
  }, category="LC_COLLATE", locale=locale)

  counts <- unlist(counts)
  scores <- counts / (n-1)

  if (sort) {
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


setMethodS3("typeOfSequenceOrdering", "numeric", function(values, sort=TRUE, as=c("scores", "humanreadable"), ...) {
  # Argument 'as':
  as <- match.arg(as)

  n <- length(values)
  counts <- list()

  ## (a) Ordered by chromosome length?
  delta <- diff(order(values))
  counts$length <- sum(delta == 1L, na.rm=TRUE)

  counts <- unlist(counts)
  scores <- counts / (n-1)

  names <- names(values)
  if (!is.null(names)) {
    scoresN <- typeOfSequenceOrdering(names, sort=FALSE, as="scores", ...)
    scores <- c(scoresN, scores)
  }

  if (sort) {
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
