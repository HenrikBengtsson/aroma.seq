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
#  \item{sample}{Specifies the \code{SM} read group.}
#  \item{library}{Specifies the \code{LB} read group.}
#  \item{platform}{Specifies the \code{PL} read group.}
#  \item{platformUnit}{Specifies the \code{PU} read group.}
#  \item{sequencingCenter}{Specifies the \code{CN} read group.}
#  \item{description}{Specifies the \code{DS} read group.}
#  \item{runDate}{Specifies the \code{DT} read group.}
#  \item{flowOrder}{Specifies the \code{FO} read group.}
#  \item{keySequence}{Specifies the \code{KS} read group.}
#  \item{program}{Specifies the \code{PG} read group.}
#  \item{predictedInsertSize}{Specifies the \code{PI} read group.}
#  \item{...}{Additional named arguments, including two-letter read 
#    group keys for the above, e.g. 'SM'.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
# 
# \references{
#  [1] The SAM Format Specification Working Group,
#      \emph{The SAM Format Specification}, Sept 7, 2011.\cr
# }
#*/###########################################################################
setConstructorS3("SamReadGroup", function(sample=NULL, library=NULL, platform=NULL, platformUnit=NULL, sequencingCenter=NULL, description=NULL, runDate=NULL, flowOrder=NULL, keySequence=NULL, program=NULL, predictedInsertSize=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'sample':
  sample <- Arguments$getCharacter(sample);

  # Argument 'library':
  library <- Arguments$getCharacter(library);

  # Argument 'platform':
  platform <- Arguments$getCharacter(platform);

  # Argument 'platformUnit':
  platformUnit <- Arguments$getCharacter(platformUnit);

  # Argument 'sequencingCenter':
  sequencingCenter <- Arguments$getCharacter(sequencingCenter);

  # Argument 'description':
  description <- Arguments$getCharacter(description);

  # Argument 'runDate':
  if (inherits(runDate, "Date")) {
    runDate <- Arguments$getCharacter(runDate);
  } else if (inherits(runDate, "POSIXct") || inherits(runDate, "POSIXlt")) {
    runDate <- Arguments$getCharacter(runDate);
  } else {
    runDate <- Arguments$getCharacter(runDate);
  }


  extend(Object(), "SamReadGroup",
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
  );
})


setMethodS3("as.character", "SamReadGroup", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, asString(this, fmtstr=" %s:%s", ...));
  class(s) <- c("GenericSummary", class(s));
  s;
})


setMethodS3("byScanBamHeader", "SamReadGroup", function(static, header, ...) {
  text <- header$text;
  keys <- names(text);
  rgList <- text[(keys == "@RG")];
  rgList <- lapply(rgList, FUN=function(params) {
    pattern <- "([^:]*):(.*)";
    keys <- gsub(pattern, "\\1", params);
    values <- gsub(pattern, "\\2", params);
    names(values) <- keys;
    values;
  });

  resList <- vector("list", length=length(rgList));
  for (kk in seq(along=rgList)) {
    rg <- rgList[[kk]];
    res <- newInstance(static);
    for (key in names(rg)) {
      res[[key]] <- rg[key];
    }
    resList[[kk]] <- res;
  } # for (kk ...)

  resList;
}, static=TRUE) # byScanBamHeader()



setMethodS3("assignBy", "SamReadGroup", function(this, other, ...) {
  # Argument 'other':
  other <- Arguments$getInstanceOf(other, "SamReadGroup");

  rgList <- asSamList(other);
  for (key in names(rgList)) {
    this[[key]] <- rgList[[key]];
  }

  invisible(this);
})


setMethodS3("hasSample", "SamReadGroup", function(this, ...) {
  !is.null(this$SM);
})

setMethodS3("hasLibrary", "SamReadGroup", function(this, ...) {
  !is.null(this$LB);
})

setMethodS3("hasPlatform", "SamReadGroup", function(this, ...) {
  !is.null(this$PL);
})

setMethodS3("hasPlatformUnit", "SamReadGroup", function(this, ...) {
  !is.null(this$PU);
})

setMethodS3("hasSequencingCenter", "SamReadGroup", function(this, ...) {
  !is.null(this$CN);
})

setMethodS3("hasDescription", "SamReadGroup", function(this, ...) {
  !is.null(this$DS);
})

setMethodS3("hasRunDate", "SamReadGroup", function(this, ...) {
  !is.null(this$DT);
})


setMethodS3("asSamList", "SamReadGroup", function(this, drop=TRUE, ...) {
  env <- getEnvironment(this);
  res <- as.list(env);
  keep <- (nchar(names(res)) == 2L);
  res <- res[keep];
  if (drop) {
    keep <- !sapply(res, FUN=is.null);
    res <- res[keep];
  }
  res <- lapply(res, FUN=unname);
  res;
})


setMethodS3("asString", "SamReadGroup", function(this, fmtstr="@RG\t%s:%s", collapse=NULL, ...) {
  res <- asSamList(this, ...);
  keys <- names(res);
  res <- unlist(res, use.names=TRUE);
  res <- sprintf(fmtstr, keys, res);
  names(res) <- keys;
  if (!is.null(collapse)) {
    res <- paste(res, collapse=collapse);
  }
  res;
})




setMethodS3("validate", "SamReadGroup", function(this, ...) {
  rgList <- asSamList(this, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate platform
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  platform <- rgList[["PL"]];
  if (!is.null(platform)) {
    knownPlatforms <- c("CAPILLARY", "HELICOS", "ILLUMINA", "IONTORRENT",
                        "LS454", "PACBIO", "SOLID");
    if (!is.element(toupper(platform), knownPlatforms)) {
      throw("Unknown 'PL' (platform) value: ", platform);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate run date
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  runDate <- rgList[["DT"]];
  # TODO...

  invisible(TRUE);
})


############################################################################
# HISTORY:
# 2012-09-28
# o Created.
############################################################################
