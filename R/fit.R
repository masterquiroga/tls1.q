# PQ-TLS Conway-Maxwell-Poisson fitter
#' 
#' Fits PQ-TLS stage clock times as COM-Poisson.
#' 
#' @param file file to read data
#' @param name name of the data
#' @param estimate_auth logical; if TRUE, authentication penalty times are estimated
#'
#' @references 
#' 
#' @references
#' Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
#' Count Data. Annals of Applied Statistics, 4(2), 943-961.
#' 
#' 
#' 
#' @author Víctor G. G. Quiroga
#' @keywords COM-Poisson distribution, PQ-TLS
#' @name fit
NULL

library("tidyverse")
library("fitdistrplus")
library("COMPoissonReg")
#library("compoisson")

#source("./R/com.R")

optim_methods <- c(
  "L-BFGS-B", 
  "SANN", # Simulated-annealing
  "Nelder-Mead", # Nelder-Mead
  "CG" # conjugate gradients
)

distros <- c(
  'pois',
  'cmp',
  'tweedie',
  'bnbinom'
)

search <- function(x, params) {
    options(warn = -1)
    result = optim(c(params$lambda, params$nu), function(p) {
        return(-compoisson::com.loglikelihood(x, p[1], p[2]))
    }, method = "L-BFGS-B", lower = c(1e-10, 1e-10))
    options(warn = 0)
    lambda = result$par[1]
    nu = result$par[2]
    fits = list(lambda = lambda, nu = nu, z = compoisson::com.compute.z(lambda, 
        nu), fitted.values = sum(x[, 2]) * compoisson::dcom(x[, 1], lambda, 
        nu), log.likelihood = compoisson::com.loglikelihood(x, lambda, nu))
    return(fits)
}

fit <- function(file, name, estimate_auth = FALSE) {
  dir <- getwd()
  data <- read_csv(file, 
    col_names = c(
      "auth",
      "exchange",
      "cypher",
      "total"
    ))
  if (estimate_auth) {
    data <- data %>% mutate(auth = total - exchange - cypher)
  }

  fits <- c("total") %>% map(function(variable) {
    cat(date(),"Fitting ", variable, "\n")

    
    


    cat(date(),"  Finding starting Poisson models\n")
    freqs <- plyr::count(data[[variable]])
    (params <- list(
      lambda = (freqs[,1] %*% freqs[,2])/sum(freqs[,2]),
      nu = 1
    ))
    pmle <- data[[variable]] %>% fitdist(
          distr = 'pois', 
          method = 'mle', 
          start = list(
            lambda = params$lambda
          ),
          lower = c(1e-10),
          control = list(maxit = 100000))
    mle <- pmle
    params <- mle$estimate %>% as.list

    F <- ecdf(data[[variable]])
    probs <- c(boxplot.stats(data[[variable]])$conf %>% F) %>% sort
    pqme <- data[[variable]] %>% fitdist(
          distr = 'pois', 
          method = 'qme', 
          start = list(
            lambda = params$lambda
          ),
          probs = probs[2],
          control = list(maxit = 100000))
    qme <- pqme
    params <- qme$estimate %>% as.list
    #mme <- NA

    cat(date(),"  Finding optimal starting parameters via Maximum Likelihood\n")
    (params <- compoisson::com.fit(freqs))
    setwd("./data/output/")
    sink(paste0(name,"_",variable,"_params.txt"))
    print(params)
    sink()
    setwd(dir)

    cat(date(),"  Fitting CMP distribution via Maximum Likelihood\n")
    for (optim.method in optim_methods) {
      tryCatch({
        cat(date(),"    Attempting fit with ", optim.method, " optimization algorithm\n")
        (mle <- data[[variable]] %>% fitdist(
          distr = 'cmp', 
          method = 'mle', 
          start = list(
            lambda = params$lambda,
            nu = params$nu
          ),
          optim.method = optim.method,
          lower = c(1e-10, 1e-10),
          control = list(maxit = 100000))) %>% print
        setwd("./data/output/")
        sink(paste0(name,"_",variable,"_mle.txt"))
        mle %>% summary %>% print
        sink()
        setwd(dir)
        params <- mle$estimate %>% as.list
        break
      }, error = function(err) {
        cat(date(),"    Failed to converge...\n")
        warning(err)
      })
    }

    cat(date(),"  Fitting CMP distribution via Quantile Matching\n")
    
    for (optim.method in optim_methods) {
      tryCatch({
        cat(date(),"    Attempting fit with ", optim.method, " optimization algorithm\n")
        (qme <- data[[variable]] %>% fitdist(
          distr = 'cmp', 
          method = 'qme', 
          start = list(
            lambda = params$lambda,
            nu = params$nu
          ),
          probs = probs,
          optim.method = optim.method,
          lower = c(1e-10, 1e-10),
          control = list(maxit = 1000000))) %>% print
        setwd("./data/output/")
        sink(paste0(name,"_",variable,"_qme.txt"))
        qme %>% summary %>% print
        sink()
        setwd(dir)
        params <- qme$estimate %>% as.list
        break
      }, error = function(err) {
        cat(date(),"    Failed to converge...\n")
        warning(err)
      })
    }

    # cat(date(),"  Fitting CMP distribution via Matching Moments\n")
    # for (optim.method in optim_methods) {
    #   tryCatch({
    #     cat(date(),"    Attempting fit with ", optim.method, " optimization algorithm\n")
    #     (mme <- data[[variable]] %>% fitdist(
    #       distr = 'cmp', 
    #       method = 'mme', 
    #       start = list(
    #         lambda = params$lambda,
    #         nu = params$nu
    #       ),
    #       order = 1:2,
    #       memp = "memp",
    #       optim.method = optim.method,
    #       control = list(maxit = 10000)))
    #     break
    #   }, error = function(err) {
    #     cat(date(),"    Failed to converge...\n")
    #     warning(err)
    #     next
    #   })
    # }

    models <- list(pmle, pqme, mle, qme)

    try({
      legends <- c(
        "Poisson Maximum Likelihood", 
        "Poisson Quantile Matching", 
        "COM-Poisson Maximum Likelihood", 
        "COM-Poisson Quantile Matching"
      )
      setwd("./data/output/")
      cat(date(),"  Selecting best model\n")
      png(
        filename = paste0(name,"_",variable,"_gof.png"),
        width = 1600, height = 800)
      par(mfrow = c(1,2))
      denscomp(
        models, 
        legendtext = legends, 
        discrete = TRUE,
        fitlty = 1)
      cdfcomp(
        models, 
        legendtext = legends, 
        discrete = TRUE,
        fitlty = 1)
      dev.off()

      sink(paste0(name,"_",variable,"_gof.txt"))
      (gof <- gofstat(models, fitnames = legends)) %>% print
      sink()

      print(gof)
      cat(date(),"  Best model for ",variable," found in ", names(which(gof$bic == min(gof$bic))), "\n")

      sink(paste0(name,"_",variable,"_gof.txt"))
      (gof <- gofstat(models, fitnames = legends)) %>% print
      sink()

    })

    setwd(dir)
    return(list(
        pmle,
        pmle,
        mle,
        qme))
  })
  return(fits)
}


x86 <- list(
  no = fit("./data/input/new/x86/No-Dilitium-r2", "x86_no_auth", TRUE),
  single = fit("./data/input/new/x86/Single-Dilitium-r2", "x86_single_auth", FALSE),
  double = fit("./data/input/new/x86/Double-Dilitium-r2", "x86_double_auth", FALSE)
)


arm <- list(
  no = fit("./data/input/new/arm/No-Dilitium", "arm_no_auth", TRUE),
  single = fit("./data/input/new/arm/Single-Dilitium", "arm_single_auth", FALSE),
  double = fit("./data/input/new/arm/Double-Dilitium", "arm_double_auth", FALSE)
)
