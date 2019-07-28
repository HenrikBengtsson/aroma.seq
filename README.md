# aroma.seq: High-Throughput Sequence Analysis using the Aroma Framework


## Installation
R package aroma.seq is only available via [GitHub](https://github.com/HenrikBengtsson/aroma.seq) and can be installed in R as:
```r
remotes::install_github("HenrikBengtsson/aroma.seq")
```

### Pre-release version

To install the pre-release version that is available in Git branch `develop` on GitHub, use:
```r
remotes::install_github("HenrikBengtsson/aroma.seq@develop")
```
This will install the package from source.  



## Contributions

This Git repository uses the [Git Flow](http://nvie.com/posts/a-successful-git-branching-model/) branching model (the [`git flow`](https://github.com/petervanderdoes/gitflow-avh) extension is useful for this).  The [`develop`](https://github.com/HenrikBengtsson/aroma.seq/tree/develop) branch contains the latest contributions and other code that will appear in the next release, and the [`master`](https://github.com/HenrikBengtsson/aroma.seq) branch contains the code of the latest release.

Contributing to this package is easy.  Just send a [pull request](https://help.github.com/articles/using-pull-requests/).  When you send your PR, make sure `develop` is the destination branch on the [aroma.seq repository](https://github.com/HenrikBengtsson/aroma.seq).  Your PR should pass `R CMD check --as-cran`, which will also be checked by <a href="https://travis-ci.org/HenrikBengtsson/aroma.seq">Travis CI</a> and <a href="https://ci.appveyor.com/project/HenrikBengtsson/aroma-seq">AppVeyor CI</a> when the PR is submitted.


## Software status

| Resource:     | GitHub        | Travis CI       | AppVeyor         |
| ------------- | ------------------- | --------------- | ---------------- |
| _Platforms:_  | _Multiple_          | _Linux & macOS_ | _Windows_        |
| R CMD check   |  | <a href="https://travis-ci.org/HenrikBengtsson/aroma.seq"><img src="https://travis-ci.org/HenrikBengtsson/aroma.seq.svg" alt="Build status"></a>   | <a href="https://ci.appveyor.com/project/HenrikBengtsson/aroma-seq"><img src="https://ci.appveyor.com/api/projects/status/github/HenrikBengtsson/aroma.seq?svg=true" alt="Build status"></a> |
| Test coverage |                     | <a href="https://codecov.io/gh/HenrikBengtsson/aroma.seq"><img src="https://codecov.io/gh/HenrikBengtsson/aroma.seq/branch/develop/graph/badge.svg" alt="Coverage Status"/></a>     |                  |
