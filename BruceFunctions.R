#### Variable computing functions for data.table ####

# Examples:
#
library(data.table)
library(dplyr)
library(stringr)
# d=data.table(x1=1:5, x4=c(2,2,5,4,5), x3=c(3,2,NA,NA,5), x2=c(4,4,NA,2,5), x5=c(5,4,1,4,5)); d
# d[,":="(n.na=COUNT(d, "x", 1:5, value=NA),
#         n.2=COUNT(d, "x", 1:5, value=2),
#         sum=SUM(d, "x", 1:5),
#         mean1=MEAN(d, "x", 1:5),
#         mean2=MEAN(d, vars=c("x1", "x4")),
#         mean3=MEAN(d, varrange="x1:x2", rev="x2", likert=1:5),
#         cons1=CONSEC(d, "x", 1:5),
#         cons2=CONSEC(d, varrange="x1:x5"))]; d


RECODE=car::recode

COUNT=function(data, var=NULL, items=NULL,
               vars=NULL,
               varrange=NULL,
               value=NA) {
  Count=function(...) sum(c(...), na.rm=TRUE)
  if(!is.null(varrange)) {
    dn=names(data)
    if(length(varrange)==1) varrange=strsplit(varrange, ":")[[1]]
    varMin=varrange[1]
    varMax=varrange[length(varrange)]
    vars=dn[which(dn==varMin):which(dn==varMax)]
  } else {
    if(is.null(vars)) vars=paste0(var, items)
  }
  vars=paste(deparse(substitute(data)), vars, sep="$")
  if(is.na(value)) {
    varlist=paste0("is.na(", vars, ")")
  } else {
    varlist=paste0(vars, "==", value)
  }
  eval(parse(text=paste0("mapply(Count, ", paste(varlist, collapse=", "), ")")))
}


SUM=function(data, var=NULL, items=NULL,
             vars=NULL,
             varrange=NULL,
             rev=NULL, likert=NULL,
             na.rm=TRUE) {
  Sum=function(...) sum(..., na.rm=na.rm)
  if(!is.null(varrange)) {
    dn=names(data)
    if(length(varrange)==1) varrange=strsplit(varrange, ":")[[1]]
    varMin=varrange[1]
    varMax=varrange[length(varrange)]
    vars=dn[which(dn==varMin):which(dn==varMax)]
  } else {
    if(is.null(vars)) vars=paste0(var, items)
  }
  if(is.character(rev)) rev=which(vars %in% rev)
  vars=paste(deparse(substitute(data)), vars, sep="$")
  pre=rep("", length(vars))
  pre[rev]=ifelse(is.null(likert), "", paste0(min(likert)+max(likert), "-"))
  varlist=paste0(pre, vars)
  eval(parse(text=paste0("mapply(Sum, ", paste(varlist, collapse=", "), ")")))
}


MEAN=function(data, var=NULL, items=NULL,
              vars=NULL,
              varrange=NULL,
              rev=NULL, likert=NULL,
              na.rm=TRUE) {
  Mean=function(...) mean(c(...), na.rm=na.rm)
  if(!is.null(varrange)) {
    dn=names(data)
    if(length(varrange)==1) varrange=strsplit(varrange, ":")[[1]]
    varMin=varrange[1]
    varMax=varrange[length(varrange)]
    vars=dn[which(dn==varMin):which(dn==varMax)]
  } else {
    if(is.null(vars)) vars=paste0(var, items)
  }
  if(is.character(rev)) rev=which(vars %in% rev)
  vars=paste(deparse(substitute(data)), vars, sep="$")
  pre=rep("", length(vars))
  pre[rev]=ifelse(is.null(likert), "", paste0(min(likert)+max(likert), "-"))
  varlist=paste0(pre, vars)
  eval(parse(text=paste0("mapply(Mean, ", paste(varlist, collapse=", "), ")")))
}


STD=function(data, var=NULL, items=NULL,
             vars=NULL,
             varrange=NULL,
             rev=NULL, likert=NULL,
             na.rm=TRUE) {
  Std=function(...) sd(c(...), na.rm=na.rm)
  if(!is.null(varrange)) {
    dn=names(data)
    if(length(varrange)==1) varrange=strsplit(varrange, ":")[[1]]
    varMin=varrange[1]
    varMax=varrange[length(varrange)]
    vars=dn[which(dn==varMin):which(dn==varMax)]
  } else {
    if(is.null(vars)) vars=paste0(var, items)
  }
  if(is.character(rev)) rev=which(vars %in% rev)
  vars=paste(deparse(substitute(data)), vars, sep="$")
  pre=rep("", length(vars))
  pre[rev]=ifelse(is.null(likert), "", paste0(min(likert)+max(likert), "-"))
  varlist=paste0(pre, vars)
  eval(parse(text=paste0("mapply(Std, ", paste(varlist, collapse=", "), ")")))
}


CONSEC=function(data, var=NULL, items=NULL,
                vars=NULL,
                varrange=NULL,
                values=0:9) {
  Conseq=function(string, number=values) {
    # Consecutive Identical Digits
    library(stringr)
    pattern=paste(paste0(number, "{2,}"), collapse="|")
    ifelse(grepl(pattern, string), max(nchar(str_extract_all(string=string, pattern=pattern, simplify=TRUE))), 0)
  }
  if(!is.null(varrange)) {
    dn=names(data)
    if(length(varrange)==1) varrange=strsplit(varrange, ":")[[1]]
    varMin=varrange[1]
    varMax=varrange[length(varrange)]
    vars=dn[which(dn==varMin):which(dn==varMax)]
  } else {
    if(is.null(vars)) vars=paste0(var, items)
  }
  vars=paste(deparse(substitute(data)), vars, sep="$")
  varlist=vars
  eval(parse(text=paste0("mapply(Conseq, paste0(", paste(varlist, collapse=", "), "))")))
}



#### Basic analysis functions ####


p.z=function(z) pnorm(abs(z), lower.tail=FALSE)*2

p.t=function(t, df) pt(abs(t), df, lower.tail=FALSE)*2

p.f=function(f, df1, df2) pf(f, df1, df2, lower.tail=FALSE)

p.r=function(r, n) p.t(r/sqrt((1-r^2)/(n-2)), n-2)

p.chi2=function(chi2, df) pchisq(chi2, df, lower.tail=FALSE)


CI=function(var) {
  output=sprintf("Mean = %.2f, 95%% CI = [%.2f, %.2f]",
                 mean(var), quantile(var, 0.025), quantile(var, 0.975))
  cat(output, "\n")
}


freq=function(var, label=NULL, sort=FALSE, plot=FALSE, digits=1) {
  tableVar=table(var)
  N=sum(tableVar)
  if(is.null(label)) label=names(tableVar)
  if(plot) hist(var, xlab=deparse(substitute(var)),
                main=paste("Histogram of", deparse(substitute(var))))
  output=cbind(matrix(tableVar,
                      dimnames=list(label, "N")),
               matrix(round(tableVar/N*100, digits),
                      dimnames=list(label, "%")))
  if(sort==FALSE) print(output)
  if(sort=="-") print(output[order(output[,"N"], decreasing=T),])
  if(sort=="+") print(output[order(output[,"N"], decreasing=F),])
  cat("\nTotal N:", N, "\n")
  invisible(output)
}


corr=function(data, method="pearson", short=FALSE, plot=FALSE, pcorr=FALSE) {
  library(psych)
  cr=corr.test(data, adjust="none", method=method)
  print(cr, digits=4, short=short)
  if(plot) {
    cor.plot(cr$r, adjust="none", numbers=TRUE, diag=FALSE, pval=cr$p, stars=TRUE)
  }
  if(pcorr) {
    library(corpcor)
    print(cor2pcor(cor(data)), digits=4)
  }
}


Alpha=function(data, var, items, rev=NULL) {
  library(jmv)
  vars=paste0(var, items)
  if(!is.null(rev)) rev=paste0(var, rev)
  reliability(data, vars=eval(vars), revItems=eval(rev),
              meanScale=TRUE, sdScale=TRUE,
              itemRestCor=TRUE, alphaItems=TRUE)
}
