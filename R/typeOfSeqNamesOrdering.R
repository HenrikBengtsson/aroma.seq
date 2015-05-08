setMethodS3("typeOfSequenceOrdering", "character", function(values, unique=FALSE, rank=TRUE, as=c("scores", "humanreadable"), locale="C", ...) {
  # Argument 'as':
  as <- match.arg(as)

  ## Check unique values; faster if lots of duplicates
  uvalues <- unique(values)
  if (unique) {
    values <- uvalues
    map <- NULL
  } else {
    map <- match(values, table=uvalues)
  }
  
  n <- length(values)
  counts <- list()

  withLocale({
    ## (a) Lexicographic ordering, e.g. 1,10,11,...,2,20,21,22,MT,X,Y?
    o <- order(uvalues)
    if (!unique) o <- o[map]
    delta <- diff(o)
    counts$lexicograpic <- sum(delta == 1L, na.rm=TRUE)

    ## (b) Classical ordering, e.g. 1,2,...,10,11,...,21,22,X,Y,MT?
    ##     (assumes human or mouse)
    chrs <- suppressWarnings(as.integer(uvalues))
    chrs[is.element(uvalues, "X")] <- 23L
    chrs[is.element(uvalues, "Y")] <- 24L
    chrs[is.element(uvalues, c("M", "MT"))] <- 25L
    o <- order(chrs)
    if (!unique) o <- o[map]
    delta <- diff(o)
    counts$classical <- sum(delta == 1L, na.rm=TRUE)
    
    ## (c) Mixed sort ordering, e.g. 1,2,...,10,11,...,21,22,MT,X,Y?
    ## Note: Very slow for large number of items
    o <- mixedorder(uvalues)
    if (!unique) o <- o[map]
    delta <- diff(o)
    counts$mixedsort <- sum(delta == 1L, na.rm=TRUE)

  }, category="LC_COLLATE", locale=locale)

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
