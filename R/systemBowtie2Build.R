## source("systemBowtie2Build.R")

bDebugSystemBowtie2Build <- TRUE
if (bDebugSystemBowtie2Build)
{
    library(R.oo)
    library(R.filesets)
    library(R.utils)
}

setMethodS3("systemBowtie2Build", "default",
            function(

            ## Required args
            bin="bowtie2-build",     ## full pathname to bowtie2-build executable
            referenceIn,             ## comma-separated string = list of files with ref sequences
            bt2IndexBase,            ## write .bt2 data to files with this dir/basename
            ## *** Bowtie 2 indexes work only with v2 (not v1).  Likewise for v1 indexes. ***

            ## Options:
            f = NULL,               ## reference files are Fasta (default)
            c = NULL,               ## reference sequences given on cmd line (as <seq_in>)
            a = NULL,               ## disable automatic -p/--bmax/--dcv memory-fitting
            noauto = NULL,          ## disable automatic -p/--bmax/--dcv memory-fitting
            p = NULL,               ## use packed strings internally; slower, uses less mem
            packed = NULL,          ## use packed strings internally; slower, uses less mem
            bmax = NULL,            ## <int> = max bucket sz for blockwise suffix-array builder
            bmaxdivn = NULL,        ## <int> = max bucket sz as divisor of ref len (default: 4)
            dcv = NULL,             ## <int> = diff-cover period for blockwise (default: 1024)
            nodc = NULL,            ## disable diff-cover (algorithm becomes quadratic)
            r = NULL,               ## don't build .3/.4.bt2 (packed reference) portion
            noref = NULL,           ## don't build .3/.4.bt2 (packed reference) portion
            three = NULL,           ## just build .3/.4.bt2 (packed reference) portion  ## IS THIS ALLOWED?
            ## - The actual bowtie2-build switch is '-3'
            justref = NULL,         ## just build .3/.4.bt2 (packed reference) portion
            o = NULL,               ## <int>: SA is sampled every 2^offRate BWT chars (default: 5)
            offrate = NULL,         ## <int>: SA is sampled every 2^offRate BWT chars (default: 5)
            t = NULL,               ## <int>: # of chars consumed in initial lookup (default: 10)
            ftabchars = NULL,       ## <int>: # of chars consumed in initial lookup (default: 10)
            seed = NULL,            ## <int>: ## seed for random number generator
            q = NULL,               ## verbose output (for debugging)
            quiet = NULL,           ## verbose output (for debugging)
            h = NULL,               ## print detailed description of tool and its options
            help = NULL,            ## print detailed description of tool and its options
            usage = NULL,           ## print this usage message
            version = NULL,         ## print version information and quit
            ...,
            overwrite=FALSE,
            verbose=FALSE
            )
        {
            ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            argsList <- list(...)  ## This automatically leaves out the NULL args
            if (length(argsList) > 0)
            {
                switchVals <- unlist(argsList)
                switchNames <- names(argsList)
                hyphenVec <- sapply(switchNames, function(x) {ifelse(nchar(x) == 1, "-", "--")})
                cmdSwitches <- paste(paste(hyphenVec, switchNames, sep=""), switchVals, collapse=" ")
            } else {
                cmdSwitches <- ""
            }

            binWithArgs <- sprintf("%s", bin)
            binWithArgs <- sprintf("%s %s", binWithArgs, cmdSwitches)
            binWithArgs <- sprintf("%s %s", binWithArgs, referenceIn)
            binWithArgs <- sprintf("%s %s", binWithArgs, bt2IndexBase)

            ## cmd <- sprintf("%s", binWithArgs)  ## In emacs/ess, shQuote() version with single quotes added does not run
            cmd <- binWithArgs

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

            return(res)
        })

############################################################################
# HISTORY:
# 2012-08-22
# o TT:  First implementation of low-level system wrapper, including all bowtie2-build switches; not tested
# 2012-08-21
# o TT:  Implemented working version (turns out this was closer in intent to bowtie2Build.R
# 2012-08-20
# o HB:  Created systemBowtie2Build stub
############################################################################
