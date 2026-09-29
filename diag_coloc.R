#!/usr/bin/env Rscript
# Quick diagnostic: inspect coloc.susie return structure on chr9
library(data.table)
library(susieR)
library(coloc)

WORK <- "D:/IA_AAA_analysis/results/susie/work"
OUT  <- "D:/IA_AAA_analysis/results/susie"

lchr <- "9"; lstart <- 21600000; lend <- 22600000

bim <- fread(sprintf("%s/region_chr%s_%d_%d.bim", WORK, lchr, lstart, lend),
             col.names=c("chr","rsid","cM","pos","refA1","refA2"))
vars <- scan(sprintf("%s/ld_chr%s_%d_%d.unphased.vcor1.vars", WORK, lchr, lstart, lend), what=character(), quiet=TRUE)
Rfull <- matrix(scan(sprintf("%s/ld_chr%s_%d_%d.unphased.vcor1", WORK, lchr, lstart, lend), what=numeric(), quiet=TRUE), nrow=length(vars), byrow=TRUE)
dimnames(Rfull) <- list(vars, vars)

ia  <- fread(sprintf("%s/region_chr%s_%d_%d_IA.tsv", OUT, lchr, lstart, lend))
aaa <- fread(sprintf("%s/region_chr%s_%d_%d_AAA.tsv", OUT, lchr, lstart, lend))
ia  <- ia[!is.na(beta)&!is.na(se)&se>0&!is.na(maf)&maf>0&maf<0.5]; ia[,pos:=as.integer(pos)]
aaa <- aaa[!is.na(beta)&!is.na(se)&se>0&!is.na(maf)&maf>0&maf<0.5]; aaa[,pos:=as.integer(pos)]
common_pos <- intersect(ia$pos, aaa$pos)
ia <- ia[match(common_pos,pos)]; aaa <- aaa[match(common_pos,pos)]
setkey(bim,pos); m <- bim[J(common_pos)]
keep <- !is.na(m$rsid); ia<-ia[keep]; aaa<-aaa[keep]; m<-m[keep]
al <- function(g1,g2,r1,r2){a1=toupper(g1);a2=toupper(g2);q1=toupper(r1);q2=toupper(r2);s=rep(NA_real_,length(a1));s[a1==q1&a2==q2]=1;s[a1==q2&a2==q1]=-1;s}
s_ia<-al(ia$a1,ia$a2,m$refA1,m$refA2); s_aaa<-al(aaa$a1,aaa$a2,m$refA1,m$refA2)
keep<-!is.na(s_ia)&!is.na(s_aaa); ia<-ia[keep];aaa<-aaa[keep];m<-m[keep];s_ia<-s_ia[keep];s_aaa<-s_aaa[keep]
z_ia<-s_ia*ia$beta/ia$se; z_aaa<-s_aaa*aaa$beta/aaa$se
snps<-m$rsid; R<-Rfull[snps,snps,drop=FALSE]; R[is.na(R)]<-0; diag(R)<-1
ia_N<-median(ia$n); aaa_N<-median(aaa$n)

cat("Fitting susie_rss with estimate_residual_variance=FALSE...\n")
fit_ia <- susie_rss(z=z_ia, R=R, L=5, estimate_residual_variance=FALSE, check_prior=FALSE, coverage=0.95)
fit_aaa <- susie_rss(z=z_aaa, R=R, L=5, estimate_residual_variance=FALSE, check_prior=FALSE, coverage=0.95)
cat("IA CS:", length(fit_ia$sets$cs), "| AAA CS:", length(fit_aaa$sets$cs), "\n")

ds_ia <- list(beta=s_ia*ia$beta, varbeta=ia$se^2, snp=snps, position=m$pos, type="cc", N=ia_N, MAF=ia$maf, s=0.33, LD=R)
ds_aaa <- list(beta=s_aaa*aaa$beta, varbeta=aaa$se^2, snp=snps, position=m$pos, type="cc", N=aaa_N, MAF=aaa$maf, s=0.07, LD=R)

cat("Running coloc.susie...\n")
cr <- coloc.susie(ds_ia, ds_aaa, susie.args=list(L=5, estimate_residual_variance=FALSE))
cat("=== names(cr) ===\n"); print(names(cr))
cat("=== class(cr) ===\n"); print(class(cr))
cat("=== summary(cr) ===\n"); print(summary(cr))
cat("=== colnames(cr$summary) ===\n"); print(colnames(cr$summary))
