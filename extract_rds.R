library(data.table)

x <- readRDS("D:/IA_AAA_analysis/results/lava/lava_partial.rds")
out.dir <- "D:/IA_AAA_analysis/results/lava"

u <- rbindlist(x$univ, fill = TRUE)
b <- rbindlist(x$bivar, fill = TRUE)
p <- rbindlist(x$pcor, fill = TRUE)

cat("univ rows:", nrow(u), "\n")
cat("bivar rows:", nrow(b), "\n")
cat("pcor rows:", nrow(p), "\n")

write.table(u, file.path(out.dir, "lava_univariate.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(b, file.path(out.dir, "lava_bivar_raw.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(p, file.path(out.dir, "lava_pcor_conditional.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("TSV files written.\n")
cat("\n--- univ head ---\n")
print(head(u))
cat("\n--- bivar head ---\n")
print(head(b))
cat("\n--- pcor head ---\n")
print(head(p))
