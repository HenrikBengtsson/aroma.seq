setMethodS3("bowtie2IndexPrefix", "default", function(pathnameFA, method=NULL, subdir="bowtie2", tags="*", ...) {
  createIndexPrefix(pathnameFA, subdir=subdir, tags=tags, asteriskTags=method, ...);
}) # bowtie2IndexPrefix()
