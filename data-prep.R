options(repos=c(CRAN="https://cloud.r-project.org"))
need <- c("wooldridge","plm","jsonlite","lmtest","sandwich")
for(p in need) if(!requireNamespace(p,quietly=TRUE)) install.packages(p)
suppressMessages({library(wooldridge);library(plm);library(jsonlite);library(lmtest);library(sandwich)})

# Resolve every output path relative to this script's location so the script
# is portable across machines and CI runners.  When run with `Rscript`,
# commandArgs(trailingOnly=FALSE) yields a "--file=…" entry whose value is the
# path of the script being executed; we strip that and use its directory as
# the project root.
script_path <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
PROJ_ROOT   <- if (length(script_path) && nzchar(script_path)) dirname(normalizePath(script_path)) else normalizePath(getwd())
OUT_JSON    <- file.path(PROJ_ROOT, "data", "wagepan.json")
dir.create(dirname(OUT_JSON), showWarnings = FALSE, recursive = TRUE)

data(wagepan)
d <- wagepan

# derive compact categoricals -------------------------------------------------
indvars <- c("agric","min","construc","trad","tra","fin","bus","per","ent","manuf","pro","pub")
d$industry <- apply(d[,indvars],1,function(r){ i<-which(r==1); if(length(i)) indvars[i[1]] else "other" })
d$region <- with(d, ifelse(nrtheast==1,"northeast",ifelse(nrthcen==1,"northcentral",ifelse(south==1,"south","west"))))
occvars <- paste0("occ",1:9)
d$occ <- apply(d[,occvars],1,function(r){ i<-which(r==1); if(length(i)) i[1] else 0 })

keep <- c("nr","year","lwage","educ","exper","expersq","union","married","black","hisp",
          "hours","poorhlth","industry","region","occ")
out <- d[,keep]
out$lwage   <- round(out$lwage,4)
out$wage    <- round(exp(out$lwage),2)
out$hours   <- as.integer(out$hours)
out$exper   <- as.integer(out$exper)
out$expersq <- as.integer(out$expersq)
write_json(out, OUT_JSON, dataframe="rows", digits=4, auto_unbox=TRUE)
cat("WROTE JSON rows:", nrow(out), " file size:", file.size(OUT_JSON), "bytes\n\n")

# canonical estimates ---------------------------------------------------------
pd <- pdata.frame(d, index=c("nr","year"))
cl <- function(m) coeftest(m, vcov=vcovHC(m, type="HC1", cluster="group"))

show <- function(focal, ctrls=TRUE){
  rhs <- focal
  if(ctrls) rhs <- paste(focal,"+ exper + expersq")
  f <- as.formula(paste("lwage ~", rhs))
  cat("====== focal:", focal, if(ctrls)"(+exper+expersq)"else"(bivariate)","======\n")
  for(mod in c("pooling","between","random","within")){
    m <- plm(f, data=pd, model=mod)
    co <- coef(m); se <- sqrt(diag(vcov(m)))
    if(focal %in% names(co)) cat(sprintf("  %-9s  b(%s)=% .4f  se=%.4f\n", mod, focal, co[focal], se[focal]))
    else cat(sprintf("  %-9s  b(%s)= dropped (collinear/time-invariant)\n", mod, focal))
  }
  # FD
  m <- plm(f, data=pd, model="fd"); co<-coef(m)
  cat(sprintf("  %-9s  b(%s)=%s\n","fd",focal, if(focal%in%names(co)) sprintf("% .4f",co[focal]) else "dropped"))
  cat("\n")
}
for(v in c("union","married","educ")) show(v, TRUE)

# Hausman FE vs RE for union spec
fe <- plm(lwage~union+exper+expersq, pd, model="within")
re <- plm(lwage~union+exper+expersq, pd, model="random")
cat("Hausman (union spec): "); h<-phtest(fe,re); cat(sprintf("chisq=%.3f df=%d p=%.4g\n", h$statistic, h$parameter, h$p.value))
