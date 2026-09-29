suppressPackageStartupMessages(library(GenomicSEM))
fn <- GenomicSEM::ldsc
src <- deparse(fn)
# print lines mentioning liability / prev / dnorm / qnorm
idx <- grep("liab|prev|dnorm|qnorm|K\\b|samp|pop_prev|sample_prev|pop.prev", src, ignore.case=TRUE)
cat("=== matching lines ===\n")
for (i in idx) cat(i, ":", src[i], "\n")
