for (f in c('IA_2020','AAA_2023','SBP','SmokingInit')) {
  p <- file.path('D:/IA_AAA_analysis/data/qc', paste0(f,'.sumstats.gz'))
  d <- read.table(gzfile(p), header=TRUE, nrows=200000)
  cat(f, 'median N =', median(d$N, na.rm=TRUE), ' mean N =', round(mean(d$N, na.rm=TRUE)), '\n')
}
