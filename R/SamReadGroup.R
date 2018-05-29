###########################################################################/**
# @RdocClass SamReadGroup
#
# @title "The SamReadGroup class"
#
# \description{
#  @classhierarchy
#
#  A SamReadGroup object represents a SAM Read Group.
# }
#
# @synopsis
#
# \arguments{
#  \item{id}{Specifies the \code{ID} read group field.}
#  \item{sample}{Specifies the \code{SM} read group field.}
#  \item{library}{Specifies the \code{LB} read group field.}
#  \item{platform}{Specifies the \code{PL} read group field.}
#  \item{platformUnit}{Specifies the \code{PU} read group field.}
#  \item{sequencingCenter}{Specifies the \code{CN} read group field.}
#  \item{description}{Specifies the \code{DS} read group field.}
#  \item{runDate}{Specifies the \code{DT} read group field.}
#  \item{flowOrder}{Specifies the \code{FO} read group field.}
#  \item{keySequence}{Specifies the \code{KS} read group field.}
#  \item{program}{Specifies the \code{PG} read group field.}
#  \item{predictedInsertSize}{Specifies the \code{PI} read group field.}
#  \item{...}{Additional named arguments, including two-letter read
#    group keys for the above, e.g. 'SM'.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# \references{
#  [1] The SAM Format Specification Working Group,
#      \emph{The SAM Format Specification}, Sept 7, 2011.\cr
#  [2] \emph{Collected FAQs about BAM files: What's the meaning
#            of the standard read group fields?}, GATK Forums,
#      Mar 2013. \url{http://goo.gl/PWv8MX} \cr
# }
#*/###########################################################################
setConstructorS3("SamReadGroup", function(id=NULL, sample=NULL, library=NULL, platform=NULL, platformUnit=NULL, sequencingCenter=NULL, description=NULL, runDate=NULL, flowOrder=NULL, keySequence=NULL, program=NULL, predictedInsertSize=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'id':
  id <- Arguments$getCharacter(id)

  # Argument 'sample':
  sample <- Arguments$getCharacter(sample)

  # Argument 'library':
  library <- Arguments$getCharacter(library)

  # Argument 'platform':
  platform <- Arguments$getCharacter(platform)

  # Argument 'platformUnit':
  platformUnit <- Arguments$getCharacter(platformUnit)

  # Argument 'sequencingCenter':
  sequencingCenter <- Arguments$getCharacter(sequencingCenter)

  # Argument 'description':
  description <- Arguments$getCharacter(description)

  # Argument 'runDate':
  # Should be ISO 8601 date (e.g. 2012-09-28) or date/time format
  # (e.g. 2012-09-28T02:13Z).
  if (inherits(runDate, "Date")) {
    runDate <- Arguments$getCharacter(runDate)
  } else if (inherits(runDate, "POSIXct") || inherits(runDate, "POSIXlt")) {
    runDate <- Arguments$getCharacter(runDate)
  } else {
    runDate <- Arguments$getCharacter(runDate)
  }


  extend(BasicObject(), "SamReadGroup",
    ID=id,
    SM=sample,
    LB=library,
    PL=platform,
    PU=platformUnit,
    CN=sequencingCenter,
    DS=description,
    DT=runDate,
    FO=flowOrder,
    KS=keySequence,
    PG=program,
    PI=predictedInsertSize,
    ...
  )
})


setMethodS3("as.character", "SamReadGroup", function(x, ...) {
  this <- x
  s <- sprintf("%s:", class(this)[1])
  s <- c(s, asString(this, fmtstr=" %s:%s", ...))
  GenericSummary(s)
}, protected=TRUE)


setMethodS3("byScanBamHeader", "SamReadGroup", function(static, header, ...) {
  text <- header$text
  keys <- names(text)
  rgList <- text[(keys == "@RG")]
  rgList <- lapply(rgList, FUN=function(params) {
    pattern <- "([^:]*):(.*)"
    keys <- gsub(pattern, "\\1", params)
    values <- gsub(pattern, "\\2", params)
    names(values) <- keys
    values
  })

  resList <- vector("list", length=length(rgList))
  for (kk in seq_along(rgList)) {
    rg <- rgList[[kk]]
    res <- newInstance(static)
    for (key in names(rg)) {
      res[[key]] <- rg[key]
    }
    resList[[kk]] <- res
  } # for (kk ...)

  resList
}, static=TRUE) # byScanBamHeader()



setMethodS3("merge", "SamReadGroup", function(x, y, ...) {
  # To please R CMD check
  this <- x
  other <- y

  # Argument 'other':
  other <- Arguments$getInstanceOf(other, "SamReadGroup")

  res <- this
  rgList <- asSamList(other)
  for (key in names(rgList)) {
    res[[key]] <- rgList[[key]]
  }

  res
})


setMethodS3("hasID", "SamReadGroup", function(this, ...) {
  !is.null(this$ID)
})


setMethodS3("hasSample", "SamReadGroup", function(this, ...) {
  !is.null(this$SM)
})

setMethodS3("hasLibrary", "SamReadGroup", function(this, ...) {
  !is.null(this$LB)
})

setMethodS3("hasPlatform", "SamReadGroup", function(this, ...) {
  !is.null(this$PL)
})

setMethodS3("hasPlatformUnit", "SamReadGroup", function(this, ...) {
  !is.null(this$PU)
})

setMethodS3("hasSequencingCenter", "SamReadGroup", function(this, ...) {
  !is.null(this$CN)
})

setMethodS3("hasDescription", "SamReadGroup", function(this, ...) {
  !is.null(this$DS)
})

setMethodS3("hasRunDate", "SamReadGroup", function(this, ...) {
  !is.null(this$DT)
})


setMethodS3("isEmpty", "SamReadGroup", function(this, ...) {
  (length(asSamList(this, ...)) == 0L)
})

setMethodS3("asSamList", "SamReadGroup", function(this, drop=TRUE, ...) {
  res <- attributes(this)
  keep <- (nchar(names(res)) == 2L)
  res <- res[keep]
  if (drop) {
    keep <- !sapply(res, FUN=is.null)
    res <- res[keep]
  }
  res <- lapply(res, FUN=unname)
  res
})


setMethodS3("asString", "SamReadGroup", function(this, fmtstr="@RG\t%s:%s", collapse=NULL, ...) {
  res <- asSamList(this, ...)
  keys <- names(res)
  res <- unlist(res, use.names=TRUE)
  res <- sprintf(fmtstr, keys, res)
  names(res) <- keys
  if (!is.null(collapse)) {
    res <- paste(res, collapse=collapse)
  }
  res
})




setMethodS3("validate", "SamReadGroup", function(this, ...) {
  rgList <- asSamList(this, ...)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Platform
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  platform <- rgList[["PL"]]
  if (!is.null(platform)) {
    knownPlatforms <- c("CAPILLARY", "HELICOS", "ILLUMINA", "IONTORRENT",
                        "LS454", "PACBIO", "SOLID")
    if (!is.element(toupper(platform), knownPlatforms)) {
      throw("Unknown 'PL' (platform) value: ", platform)
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Platform unit
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(platform)) {
    if (toupper(platform) == "ILLUMINA") {
      platformUnit <- rgList[["PL"]]
      if (!is.null(platformUnit)) {
        # convention: Illumina flowcell barcode suffixed with a period and
        # the lane number (and further suffixed with period followed by
        # sample member name for pooled runs) [From NHI/SRA below]
        # TODO ...
        # Example(?): <flowcell barcode>.<lane nbr>[.<sample name>]
      }
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Run date
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  runDate <- rgList[["DT"]]
  # Should be ISO 8601 date (e.g. 2012-09-28) or date/time format
  # (e.g. 2012-09-28T02:13Z).
  # TODO...

  invisible(TRUE)
})
