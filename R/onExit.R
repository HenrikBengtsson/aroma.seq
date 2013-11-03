onExit <- function(expr=NULL, where=c("first", "last", "replace"), ..., exits, envir=parent.frame()) {
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
    stopifnot(is.call(exits));
  }

  # Argument 'expr':
  expr <- substitute(expr);

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


# HISTORY:
# 2013-11-02:
# o Added onExit(), which requires on.exit(eval(onExit())) at the beginning.
#   See today's email to R-devel for why.
# o Created.
