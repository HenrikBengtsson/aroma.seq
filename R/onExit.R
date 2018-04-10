###########################################################################/**
# @RdocFunction onExit
#
# @title "Records expressions to be executed when the current function exits"
#
# \description{
#  @get "title" (regardless of cause).
#  This methtod differ from @see "base::on.exit" in that it adds the
#  option to \emph{prepend} an expression to already existing expressed,
#  as an alternative to the default \emph{append}.
#
#  \emph{WARNING: This function is still under development [1,2]. /HB 2013-11-08}
# }
#
# @synopsis
#
# \arguments{
#   \item{expr}{An @expression to be recorded.}
#   \item{where}{A @character string specifying whether the @expression
#    should be prepended (\code{"first"}), appended (\code{"last"}), or
#    replace already recorded @expressions (\code{"replace"}).}
#   \item{...}{Not used.}
#   \item{exits}{(optional) Already recorded.}
#   \item{envir}{The @environment in which @expressions are recorded.}
# }
#
# @author "HB"
#
# \references{
#  [1] R-devel thread 'WISHLIST: on.exit(..., add=TRUE, where="first")
#      to address common use cases', 2013-11-03.
#      \url{https://stat.ethz.ch/pipermail/r-devel/2013-November/067876.html}
#      \cr
#  [2] R-devel thread 'on.exit() & sys.on.exit(): Calling them via eval() does not	work as hoped', 2013-11-03.
#      \url{https://stat.ethz.ch/pipermail/r-devel/2013-November/067866.html}
#      \cr
# }
#
# @keyword internal
#*/###########################################################################
onExit <- function(expr=NULL, where=c("first", "last", "replace"), ..., exits, envir=parent.frame()) {
  # Argument 'expr':
  expr <- substitute(expr);

  # Argument 'exits':
  hasExits <- !missing(exits);
  if (!hasExits) {
    if (exists(".on.exit.expressions", envir=envir, inherits=FALSE)) {
      exits <- get(".on.exit.expressions", envir=envir, inherits=FALSE);
    } else {
      exits <- NULL;
    }
  }
  if (!is.null(exits)) {
    stop_if_not(is.call(exits));
  }

  # Argument 'where':
  where <- match.arg(where);
  if (is.null(exits)) where <- "replace";

  if (where == "replace") {
    exits <- expr;
  } else if (where == "last") {
    exits <- substitute({ old; new }, list(old=exits, new=expr));
  } else if (where == "first") {
    exits <- substitute({ new; old }, list(old=exits, new=expr));
  }

  if (!hasExits) {
    assign(".on.exit.expressions", exits, envir=envir, inherits=FALSE);
  }

  invisible(exits);
} # onExit()
