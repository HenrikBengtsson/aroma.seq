#
# Command line call:
#  bowtie2-build path/to/lambda_virus.fa lambda_virus
#
# R system call:
#  bowtie2Build("path/to/lambda_virus.fa");
#  bowtie2Build("path/to/lambda_virus.fa", outPath="path/to/", name="lambda_virus");
#  => calls =>
#  systemBowtie2Build("path/to/lambda_virus.fa", "lambda_virus");
#  => calls =>
#  system("bowtie2-build path/to/lambda_virus.fa lambda_virus");
#

## [
## Q's on the plate:
## - Arg list for bowtie2Build.R: inPathname can be pathname or list?
## - How to do file checking / error handling in bowtie2Build.R?  Does Arguments$getReadablePathname do a stop()?
## - Indent convention = 2 spaces?
## - use dirname instead of getParent?
## - see [] comments below
## ]

setMethodS3("bowtie2Build", "default", function(inPathname, outPath=NULL, outName=NULL, ..., overwrite=FALSE, verbose=FALSE) {
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Validate arguments
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Argument 'inPathname':
    ## inPathname <- Arguments$getReadablePathname(inPathname);

    ## Check if inPathname is readable
    if (is.list(inPathname)) {
        tmp <- sapply(inPathname, function(x)
                  {
                      tryCatch(Arguments$getReadablePathname(x), error=function(e) {e}, warning=function(w) {})
                      ## [ Q: Does this check if x is a *file*, not a dir? ]
                  })
        if (sum(tmp)==0) {throw(paste(inPathname, "contains unreadable file(s)"))}
    } else {
        ## Return error msg but leave out warnings for now
        tryCatch(inPathname <- Arguments$getReadablePathname(inPathname), error=function(e) {e}, warning=function(w) {})
    }

    ## Check outPath
    if (is.null(outPath)) {
        ifelse (is.list(inPathname),
                outPath <- dirname(inPathname[[1]]),  ## [ hack ]
                outPath <- dirname(inPathname)  ## e.g. /path/to
                )
    }
    outPath <- Arguments$getWritablePath(outPath);

    inFiles <- paste(inPathname, collapse=",")

    ## Argument 'outName':
    if (is.null(outName)) {
        ## inFilename <- basename(inPathname);   ## e.g. lambda_virus.fa
        ifelse(is.list(inPathname),
               inFilename <- basename(inPathname[[1]]),
               inFilename <- basename(inPathname) ## e.g. lambda_virus.fa
               )
        ## outName <- gsub("[.][^.]$", ".bam", inFilename);  ## [ Uh, .bam?  Where did this come from? ]
        outName <- gsub("\\.fa$", "", inFilename);  ## e.g. lambda_virus
        ## [ NEED TO TEST: Have to make sure this does not overwrite the input files! ]
    }
    outName <- Arguments$getCharacter(outName);
    ## [ Should this be Arguments$getWriteablePathname() ? ]


    ## Additional arguments
    args <- list(...);

    ## Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
        pushState(verbose);
        on.exit(popState(verbose));   ## [ << what is this ]
    }

    verbose && enter(verbose, "Running bowtie2-build");
    verbose && cat(verbose, "Input pathname: ", inPathname);  ## [ FIX: inPathname can be a list ]
    verbose && cat(verbose, "Arguments:");
    verbose && str(verbose, args);


    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Locate external software
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Locating external software");
    command <- "bowtie2-build";
    verbose && cat(verbose, "Command: ", command)

    bin <- Sys.which(command);
    verbose && cat(verbose, "Located pathname: ", bin)

    ## Assert existence
    if (identical(bin, "") || !isFile(bin)) {
        throw("Failed to located external software (via the system PATH): ", command);
    }
    verbose && exit(verbose);


    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Call external software
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Calling external software");

    ## Argument to bowtie2-index (maybe)
    outPrefix <- file.path(outPath, outName); ## e.g. path/to/lambda_virus

    cmdSwitches <- args$cmdSwitches
    ## binWithArgs <- paste(bin, cmdSwitches, inPathname, outName)
    binWithArgs <- paste(bin, cmdSwitches, inFiles, outName)

    ## cmd <- sprintf("%s", shQuote(binWithArgs));
    cmd <- sprintf("%s", binWithArgs)  ## In emacs/ess, shQuote() version with single quotes added does not run

    verbose && cat(verbose, "System call: ", cmd)



    res <- system(cmd, ...);
    verbose && exit(verbose);
    verbose && exit(verbose);
    invisible(res);
}) ## bowtie2-build()


############################################################################
# HISTORY:
# 2012-07-18
# o Created bowtie2-build() stub.
############################################################################

