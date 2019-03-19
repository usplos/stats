HLM_summary=function(model,
                     vartypes=c("intercept",
                                "L1fixed",
                                "L1random-[grouptag]-[L1varname]",
                                "L2-[grouptag]",
                                "cross-[grouptag]-[L1varname]"),
                     digits=0,
                     nsmall=3,
                     nsmall.p=5,
                     check.vartypes=FALSE) {
  library(lme4)
  library(lmerTest)
  p.trans=function(p) {
    ifelse(p<2e-16, "< 2e-16",
           ifelse(p<10^-nsmall.p, format(p, digits=2, scientific=T),
                  format(p, digits=0, nsmall=nsmall.p, scientific=F)))
  }
  sig.trans=function(p) {
    ifelse(p<.0001, "****",
           ifelse(p<.001, "***",
                  ifelse(p<.01, "**",
                         ifelse(p<.05, "*",
                                ifelse(p<.10, ".", "")))))
  }

  sumModel=summary(model, cor=F)
  paras=sumModel[["devcomp"]][["dims"]][["p"]]
  varTemplate=c("intercept", "L1fixed", "L1random-[grouptag]-[L1varname]", "L2-[grouptag]", "cross-[grouptag]-[L1varname]")
  if(identical(vartypes, varTemplate) | length(vartypes)!=paras) {
    print(sumModel)
    cat("\nPlease re-define 'vartypes'! We have", paras, "parameters (including the intercept)!\n")
    cat("You should define 'vartypes' by a vector 'c()' with the same variable order as in the 'summary' output. You may include these terms if necessary:\n")
    cat(" - 'intercept'\n - 'L1fixed'\n - 'L1random-[grouptag]-[L1varname]'\n - 'L2-[grouptag]'\n - 'cross-[grouptag]-[L1varname]' (i.e., cross-level interaction)\n")
    cat("Note. '[grouptag]' should be replaced by the clustering/grouping variable name in your data.")
    cat("'[L1varname]' should be replaced by the level-1 variable name in your data.\n")
  } else {
    ## print: Model formula ##
    cat(sumModel[["methTitle"]], "\n")
    print(sumModel[["call"]])
    cat("\n")

    # print: Information Criteria ##
    # AIC = -2LL + 2p  [p = number of parameters]
    # BIC = -2LL + p*ln(N)  [N = number of cases]
    logLik=sumModel[["logLik"]]
    print(logLik)
    cat("-2 Log Likelihood (-2LL, deviance):  ", -2*logLik, "\n")
    cat("Akaike's Information Criterion (AIC):", AIC(logLik), "\n")
    cat("Schwarz's Bayesian Criterion (BIC):  ", BIC(logLik), "\n")

    ## print: random effects ##
    cat("\nRandom effects:\n")
    RE=sumModel[["varcor"]]
    res=sumModel[["sigma"]]^2
    print(RE, comp="Variance")

    ## print: Omega^2 ##
    omg2 = 1-var(residuals(model))/(var(model.response(model.frame(model))))
    cat("Omega^2 =", omg2, "\n")

    ## print: sample size ##
    .prt.grps(ngrps=ngrps(model), nobs=nobs(model))
    # cat("\nSample size: N =", sumModel[["devcomp"]][["dims"]][["N"]])
    # cat("\n       valid n =", sumModel[["devcomp"]][["dims"]][["n"]])
    # cat("\nGroups:\n")
    # print(sumModel[["ngrps"]])

    ## print: fixed effects ##
    cat("\nFixed effects:\n")
    FE=sumModel[["coefficients"]]
    # dimnames(FE)[[2]][1]="Gamma"
    # dimnames(FE)[[2]][2]="S.E."
    # dimnames(FE)[[2]][3]="t"
    # compute df, p, r (t-to-r transformation)
    df.l1=sumModel[["devcomp"]][["dims"]][["nmp"]] # N - all parameters
    df.l2=sumModel[["ngrps"]]
    Sq=sum(grepl("L2", vartypes)) # number of level-2 predictors
    q=df.l2
    for(grouptag in names(df.l2)) {
      q[grouptag]=sum(grepl(paste0("L2-",grouptag), vartypes))
    }
    ts=as.vector(FE[,"t value"])
    dfs=ps=rs=R2s=sig=c()
    for(i in 1:paras) {
      if(vartypes[i]=="intercept") {
        # df=min(df.l2)-Sq-1
        df=NA
      } else if(vartypes[i]=="L1fixed") {
        df=df.l1
      } else {
        vartemp=strsplit(vartypes[i], "-")[[1]]
        vartype=vartemp[1]
        grouptag=vartemp[2]
        if(vartype=="L2") {
          df=df.l2[grouptag]-q[grouptag]-1
        } else {
          # vartype=="L1random" | vartype=="cross"
          l1var=vartemp[3]
          qc=sum(grepl(paste0("cross-",grouptag,"-",l1var), vartypes))
          df=df.l2[grouptag]-qc-1
        }
      }
      t=ts[i]
      p=pt(abs(t), df, lower.tail=FALSE)*2
      dfs[i]=df
      ps[i]=p
      rs[i]=sign(t)*sqrt(t^2/(t^2+df))
      R2s[i]=sqrt(t^2/(t^2+df))^2
      sig[i]=sig.trans(p)
    }
    # combine columns
    varnames=dimnames(FE)[[1]]
    FE=as.data.frame(FE)
    row.names(FE)=varnames
    names(FE)=c("Gamma", "S.E.", "df.apprx", "t", "p.apprx") # abbreviate("approx", 5)
    FE=cbind(FE[c(1,2,4,3,5)], df.HLM=dfs, p.HLM=ps, r.prtl=rs, R2.prtl=R2s, sig=sig)
    if(check.vartypes) FE=cbind(FE, vartypes)
    # digits transform
    FE$p.apprx=mapply("p.trans", FE$p.apprx)
    FE$p.HLM=mapply("p.trans", FE$p.HLM)
    FE$t=mapply("format", FE$t, digits=0, nsmall=2, scientific=F)
    FE$df.apprx=mapply("format", FE$df.apprx, digits=0, nsmall=1, scientific=F)
    FE$df.HLM=mapply("format", FE$df.HLM, digits=0, nsmall=0, scientific=F)
    # print
    print(format(FE, digits=digits, nsmall=nsmall, scientific=F), quote=F)

    ## print: notes ##
    cat("---\nSignif. codes: 0 '****' .0001 '***' .001 '**' .01 '*' .05 '.' .10 ' ' 1\n")
    cat("Note #1: df.apprx is estimated by Satterthwaite's (1946) approximation.\n")
    cat("Note #2: r.partial is calculated by t-to-r transformation.\n")
  }
}
