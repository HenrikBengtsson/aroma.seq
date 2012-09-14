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

setMethodS3("bowtie2Build", "default", function(inPathnames, outPath=NULL, outName=NULL,
                                                cmdSwitchList,
                                                ..., overwrite=FALSE, verbose=FALSE, command="bowtie2-build") {

    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Validate arguments
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    inPathnames <- sapply(inPathnames, Arguments$getReadablePathname(x))
    inPathsStr <- paste(inPathnames, collapse=",")   ## [ CHECK THIS ]

    ## Check outPath
    if (is.null(outPath)) {
        outPath <- dirname(inPathnames[1])
    }
    outPath <- Arguments$getWritablePath(outPath);
    if (is.null(outPath)) {
        outName <- gsub("\\.fa$", "", basename(inPathnames[1]))   ## e.g. lambda_virus
    }
    outPrefix <- file.path(outPath, outName); ## e.g. path/to/lambda_virus
    outPrefix <- Arguments$getWriteablePathname(outPrefix, mustNotExist=T)

    ## Additional arguments
    args <- list(...);

    ## Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
        pushState(verbose);
        on.exit(popState(verbose));   ## [ << what is this ]
    }

    verbose && enter(verbose, paste("Running", command));
    verbose && cat(verbose, "Input pathname: ", inPathname);  ## [ FIX: inPathname can be a list ]
    verbose && cat(verbose, "Arguments:");
    verbose && str(verbose, args);

    ## Locate external software
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Locating external software");
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

    systemBowtie2Build(inPathname, outPrefix, cmdSwitchList)

    verbose && exit(verbose);
    verbose && exit(verbose);
    invisible(res);
}) ## bowtie2-build()


############################################################################
# HISTORY:
# 2012-07-18
# o Created bowtie2-build() stub. (HB)
# 2012-09-14
# o First real draft of "level 2" code (TAT)
############################################################################

