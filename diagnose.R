suppressMessages(library(GenomicSEM))
suppressMessages(library(lavaan))

LDSCoutput <- readRDS("D:/IA_AAA_analysis/results/genomic_sem/LDSCoutput_4traits.rds")
cat("LDSCoutput loaded OK\n")

model3 <- '
  IA ~ SBP + SMK
  AAA ~ SBP + SMK
  IA ~~ AAA
'
fit <- usermodel(covstruc = LDSCoutput, model = model3)
cat("Model fitted\n")

cat("\n=== names(fit) ===\n")
print(names(fit))

cat("\n=== names(fit$results) ===\n")
print(names(fit$results))

cat("\n=== str(fit$results) ===\n")
str(fit$results)

cat("\n=== fit$results printed ===\n")
print(fit$results)
