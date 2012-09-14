## source("systemBowtie2Build.R")

bDebugSystemBowtie2Build <- TRUE
if (bDebugSystemBowtie2Build)
{
    library(R.oo)
    library(R.filesets)
    library(R.utils)
}

## Usage:  systemBowtie2Build referenceIn="" bt2IndexBase="" f=1 c=1 ...

setMethodS3("systemBowtie2Build", "default",
            function(

            ## Required args for bowtie2 index builder
            referenceIn,             ## comma-separated string = list of files with ref sequences
            bt2IndexBase,            ## write .bt2 data to files with this dir/basename
            ## *** Bowtie 2 indexes work only with v2 (not v1).  Likewise for v1 indexes. ***

            ## Unsupported options [ flags should be binary (boolean)-valued; check the integer-valued ones? ]
            cmdSwitchList=list(
            f = FALSE,               ## reference files are Fasta (default)
            c = FALSE,               ## reference sequences given on cmd line (as <seq_in>)
            a = FALSE,               ## disable automatic -p/--bmax/--dcv memory-fitting
            noauto = FALSE,          ## disable automatic -p/--bmax/--dcv memory-fitting
            p = FALSE,               ## use packed strings internally; slower, uses less mem
            packed = FALSE,          ## use packed strings internally; slower, uses less mem
            bmax = NULL,            ## <int> = max bucket sz for blockwise suffix-array builder
            bmaxdivn = NULL,        ## <int> = max bucket sz as divisor of ref len (default: 4)
            dcv = NULL,             ## <int> = diff-cover period for blockwise (default: 1024)
            nodc = FALSE,            ## disable diff-cover (algorithm becomes quadratic)
            r = FALSE,               ## don't build .3/.4.bt2 (packed reference) portion
            noref = FALSE,           ## don't build .3/.4.bt2 (packed reference) portion
##            three = FALSE,           ## just build .3/.4.bt2 (packed reference) portion  ## IS THIS ALLOWED?
            ## - The actual bowtie2-build switch is '-3'
            `3` = FALSE,
            justref = FALSE,         ## just build .3/.4.bt2 (packed reference) portion
            o = NULL,               ## <int>: SA is sampled every 2^offRate BWT chars (default: 5)
            offrate = NULL,         ## <int>: SA is sampled every 2^offRate BWT chars (default: 5)
            t = NULL,               ## <int>: # of chars consumed in initial lookup (default: 10)
            ftabchars = NULL,       ## <int>: # of chars consumed in initial lookup (default: 10)
            seed = NULL,            ## <int>: ## seed for random number generator
            q = FALSE,               ## verbose output (for debugging)
            quiet = FALSE,           ## verbose output (for debugging)
            h = FALSE,               ## print detailed description of tool and its options
            help = FALSE,            ## print detailed description of tool and its options
            usage = FALSE,           ## print this usage message
            version = FALSE         ## print version information and quit
            ),
            ## overwrite=FALSE, ## Put these args one level up, in bowtie2Build
            ## verbose=FALSE
            bin="bowtie2-build",     ## full pathname to bowtie2-build executable
            ...
            )
        {
            ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            argsList <- as.list(sys.call()) ## This automatically leaves out the NULL args
            argsList <- argsList[-1]  ## remove the function name

            ## ############################################
            ## ## Dummy test code
            ## if (F) {
            ## cmdSwitchList$q <- TRUE
            ## cmdSwitchList$version <- TRUE
            ## cmdSwitchList$t <- 5
            ## }
            ##
            ## cmdSwitchList <- cmdSwitchList[!sapply(cmdSwitchList, is.null)]
            ## cmdSwitchLogical <- cmdSwitchList[sapply(cmdSwitchList, function(x) {identical(x, TRUE)})]
            ## ## cmdSwitchVec <- sapply(cmdSwitchList, function(x) {x})
            ##
            ##
            ## if (F) { ## DEPRECATED
            ## optionVals <- unlist(cmdSwitchList)  ## NULL's automatically removed; [ THIS IS BAD BECAUSE VECTOR HAS TO BE SET TO A SINGLE TYPE! ]
            ## optionNames <- names(optionVals)
            ## optionsLogical <- optionNames[sapply(optionVals, function(x) {identical(x, TRUE)})]
            ## optionsLogical <- sapply(optionNames[sapply(optionVals, function(x) {identical(x, TRUE)})])
            ## if (length(optionsLogical) > 0)
            ## {
            ## hyphenVec <- sapply(optionsLogical, function(x)
            ## {
            ## ifelse(nchar(x) == 1, "-", "--")
            ## })
            ## optionsLogicalStr <- paste(paste(hyphenVec, optionsLogical, sep=""))
            ## } else {
            ## optionsLogicalStr <- NULL
            ## }
            ##
            ##
            ## optionsNonLogic <- optionNames[!is.logical(optionVals)]  ## something like this...
            ## ##strList <- lapply(1:length(optionVals), function(x)
            ## ##              {
            ## optionsStr <- paste(paste("-", namesSet))
            ## }
            ## ############################################


            if (length(argsList) > 0)  ## [ always true as it stands ]
            {
                switchVals <- unlist(argsList)  ## not used
                switchNames <- names(argsList)
                indsDrop <- which(switchNames=="bt2IndexBase"|switchNames=="referenceIn")  ## [ must be a shorter way to do this ]
                ## [ Not needed if move these to one level up in the code ]
                ## indsDrop <- c(indsDrop, which(switchNames=="overwrite"|switchNames=="verbose"))  ## non-bowtie2 args
                if (length(indsDrop) > 0)
                {
                    switchNames <- switchNames[-indsDrop]
                    switchVals <- switchVals[-indsDrop]
                }
                switchNames <- sub("^three$", "3", switchNames)  ## fix the one numerical option that can't be passed as is
                hyphenVec <- sapply(switchNames, function(x) {ifelse(nchar(x) == 1, "-", "--")})
                cmdSwitches <- paste(paste(hyphenVec, switchNames, sep=""), switchVals, collapse=" ")

                ## [ This is not correct yet - some switches take on (integer) values ]

            } else {
                cmdSwitches <- ""
            }

            binWithArgs <- sprintf("%s", bin)
            binWithArgs <- sprintf("%s %s", binWithArgs, cmdSwitches)
            binWithArgs <- sprintf("%s %s", binWithArgs, referenceIn)
            binWithArgs <- sprintf("%s %s", binWithArgs, bt2IndexBase)

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
# o TT:  First implementation of low-level system wrapper, including all bowtie2-build switches; not tested
# 2012-08-21
# o TT:  Implemented working version (turns out this was closer in intent to bowtie2Build.R
# 2012-08-20
# o HB:  Created systemBowtie2Build stub
############################################################################
