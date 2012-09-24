## source("systemBowtie2Build.R")

bDebugSystemBowtie2Build <- TRUE
if (bDebugSystemBowtie2Build)
{
    library(R.oo)
    library(R.filesets)
    library(R.utils)
}

## Usage:  systemBowtie2Build referenceIn, bt2IndexBase, optionsList, <unsupported option>=1, ...)

setMethodS3("systemBowtie2Build", "default",
            function(

            ## Required args for bowtie2 index builder
            referenceIn,             ## comma-separated string = list of files with ref sequences
            bt2IndexBase,            ## write .bt2 data to files with this dir/basename
            ## *** Bowtie 2 indexes work only with v2 (not v1). ***

            ## Supported options ( bowtiew2 options w/o values are boolean; nonlogicals are integers, so far )
            optionsList=list(
            f = FALSE,               ## reference files are Fasta (default)
            c = FALSE,               ## reference sequences given on cmd line (as <seq_in>)
            a = FALSE,               ## disable automatic -p/--bmax/--dcv memory-fitting
            noauto = FALSE,          ## disable automatic -p/--bmax/--dcv memory-fitting
            p = FALSE,               ## use packed strings internally; slower, uses less mem
            packed = FALSE,          ## use packed strings internally; slower, uses less mem
            bmax = NULL,             ## <int> = max bucket sz for blockwise suffix-array builder
            bmaxdivn = NULL,         ## <int> = max bucket sz as divisor of ref len (default: 4)
            dcv = NULL,              ## <int> = diff-cover period for blockwise (default: 1024)
            nodc = FALSE,            ## disable diff-cover (algorithm becomes quadratic)
            r = FALSE,               ## don't build .3/.4.bt2 (packed reference) portion
            noref = FALSE,           ## don't build .3/.4.bt2 (packed reference) portion
            `3` = FALSE,             ## just build .3/.4.bt2 (packed reference) portion
            justref = FALSE,         ## just build .3/.4.bt2 (packed reference) portion
            o = NULL,                ## <int>: SA is sampled every 2^offRate BWT chars (default: 5)
            offrate = NULL,          ## <int>: SA is sampled every 2^offRate BWT chars (default: 5)
            t = NULL,                ## <int>: # of chars consumed in initial lookup (default: 10)
            ftabchars = NULL,        ## <int>: # of chars consumed in initial lookup (default: 10)
            seed = NULL,             ## <int>: ## seed for random number generator
            q = FALSE,               ## verbose output (for debugging)
            quiet = FALSE,           ## verbose output (for debugging)
            h = FALSE,               ## print detailed description of tool and its options
            help = FALSE,            ## print detailed description of tool and its options
            usage = FALSE,           ## print this usage message
            version = FALSE          ## print version information and quit
            ),
            ...,
            bin="bowtie2-build",     ## full pathname to bowtie2-build executable
            )
        {
            ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            ## Create string of supported options, and unsupported options; paste together
            cmdOptions <- .optionsList2String(optionsList)
            cmdOptionsAdd <- .optionsList2String(list(...))
            cmdOptions <- paste(cmdOptions, cmdOptionsAdd)

            ## Set up command line
            binWithArgs <- sprintf("%s", bin)
            binWithArgs <- sprintf("%s %s", binWithArgs, cmdOptions)
            binWithArgs <- sprintf("%s %s", binWithArgs, referenceIn)
            binWithArgs <- sprintf("%s %s", binWithArgs, bt2IndexBase)
            binWithArgs <- gsub("  +", " ", binWithArgs)  ## prettify

            ## cmd <- sprintf("%s", binWithArgs)  ## In emacs/ess, shQuote() version with single quotes added does not run
            cmd <- binWithArgs

            verbose <- Arguments$getVerbose(verbose);
            if (verbose) {
                pushState(verbose);
                on.exit(popState(verbose));
            }

            verbose && enter(verbose, "Running bowtie2-build");
            verbose && cat(verbose, "System call: ", cmd)

            ##
            ## ADDITIONAL TESTING FOR BOWTIE2-BUILD IN HERE?
            ##

            if (bDebugSystemBowtie2Build)
            {
                res <- NULL
            } else {
                res <- system(cmd);
            }

            if (verbose) {
                verbose && exit(verbose);
            }

            return(res)
        })

############################################################################
# HISTORY:
# 2012-08-22
# o TT:  First implementation of low-level system wrapper, including all bowtie2-build options; not tested
# 2012-08-21
# o TT:  Implemented working version (turns out this was closer in intent to bowtie2Build.R
# 2012-08-20
# o HB:  Created systemBowtie2Build stub
############################################################################


.optionsList2String <- function(optionsList=list())
{
    ## Return if nothing to process
    if (length(optionsList) < 1) { return(NULL) }

    ## Split into logical (boolean) and non-logical (name-value pair) options
    tmpLogical <- unlist(optionsList[sapply(optionsList, is.logical)])
    tmpNonLogical <- unlist(optionsList[!sapply(optionsList, is.logical)])  ## NULLs are auto-removed

    ## For logical options: 1) get the names of the TRUE ones, 2) add a hyphen (or two)
    tmpLogical <- tmpLogical[tmpLogical==TRUE]
    nms <- names(tmpLogical)
    optionsLogical <-
        sapply(nms, function(nm)
           {
               str <- ifelse(nchar(nm)==1,
                             paste("-", nm, sep=""),
                             paste("--", nm, sep=""))
               return(str)
           })

    ## For "non-logical" options:
    ##  1) create two vectors: names and values; 2) paste together with appropriate hypenation
    nms <- names(tmpNonLogical)
    vals <- as.vector(tmpNonLogical)  ## (Could do some sanity checking of the vals here)
    hyphenVec <-
        sapply(nms, function(nm)
           {
               ifelse(nchar(nm)==1, "-", "--")
           })
    optionsNonLogical <- paste(paste(hyphenVec, nms, sep=""), vals, sep="=")
    ## - This works appropriately even when tmpNonLogical is NULL

    ## Combine vectors into a string and return
    cmdOptions <- paste(c(optionsLogical, optionsNonLogical), collapse=" ")
    return(cmdOptions)
}

