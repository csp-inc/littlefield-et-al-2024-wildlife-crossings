
kfoldRSF<-function (mod, k = 5, nrepet = 10, nbins = 10, jitter = TRUE, 
                    random = TRUE, method = method, x = m, form_ls = ls, reproducible = TRUE){ 
  require(glmmTMB)  
  require(tidyverse)
  
  # Function for updating model fit
  modUpdate <- function(mod, data){
    mod_up <- update(mod, used ~ ., data = data)
    mod_up$parameters$theta[1] <- log(1e3)
    ntheta = length(mod_up$parameters$theta)
    mod_up$mapArg <- list(theta=factor(c(NA, rep(1,ntheta-1))))
    mod_up_fit <- glmmTMB:::fitTMB(mod_up)
    return(mod_up_fit)
  }

  dt <- model.frame(mod)
  kfold <- rd <- numeric(length = nrepet)
  resp <- as.character(attr(terms(mod), "variables"))[attr(terms(mod), 
                                                           "response") + 1]
  for (i in 1:nrepet) {
    print(paste("loop position:", i, "of", nrepet))
    dt$sets <- "train"
    if (reproducible) 
      set.seed(i)
    dt$sets[sample(which(dt[, resp] == 1), sum(dt[, resp] == 
                                                 1)/k)] <- "test"
    reg <- modUpdate(mod, data = subset(dt, sets == "train"))
    cof <- coef(reg)$cond$nid
    cof[is.na(cof)] <- 0
    predall <- exp(as.numeric(model.matrix(terms(reg), dt) %*% fixef(reg)$cond))
    if (jitter) {
      if (reproducible) 
        set.seed(i)
      predall <- jitter(predall)
    }
    quant <- quantile(predall[dt[, resp] == 0], probs = seq(from = 0, 
                                                            to = 1, length.out = nbins + 1))
    quant[1] <- -Inf
    quant[length(quant)] <- Inf
    int <- factor(findInterval(predall[dt$sets == "test"], 
                               quant), levels = 1:nbins)
    kfold[i] <- cor(1:nbins, table(int), method = "spearman")
    if (random) {
      if (reproducible) 
        set.seed(i)
      dt$sets[sample(which(dt[, resp] == 0), sum(dt[, resp] == 
                                                   1)/k)] <- "rd"
      int <- factor(findInterval(predall[dt$sets == "rd"], 
                                 quant), levels = 1:nbins)
      rd[i] <- cor(1:nbins, table(int), method = "spearman")
    }
  }
  if (random) 
    return(data.frame(kfold = c(kfold, rd), type = rep(c("obs", 
                                                         "rand"), each = nrepet)))
  else return(kfold)
  }


kfoldSSF <- function(mod, strata_name = "step_id", k = 5, nrepet = 100, jitter = FALSE,
                     reproducible = FALSE, details = FALSE){
    require(glmmTMB)
    require(tidyverse)
    
    # Get model frame
    dt <- model.frame(mod)
    names(dt)[1] <- "case"
    # Remove case-control sets with no case...
    full <- dt %>% group_by(step_id) %>% summarize(case_sum = sum(case))
    keeps <- full$step_id[full$case_sum > 0]
    dt <- dt %>% dplyr::filter(step_id %in% keeps)
    
    # Function for updating model fit
    modUpdate <- function(mod, data){
      mod_up <- update(mod, case ~ ., data = data)
      mod_up$parameters$theta[1] <- log(1e3)
      ntheta = length(mod_up$parameters$theta)
      mod_up$mapArg <- list(theta=factor(c(NA, rep(1,ntheta-1))))
      mod_up_fit <- glmmTMB:::fitTMB(mod_up)
      return(mod_up_fit)
    }
    
    kfold <- rd <- warn <- numeric(length = nrepet)
    if (details)
      dbg <- list()
    
    # Loop through folds and refit model
    print(paste("start:", Sys.time()))
    for (i in 1:nrepet) {
      print(paste("loop position:", i, "of", nrepet))
      
      flag = 1
      while(flag == 1){
        dt$sets <- "train"
        
        # Take a random subset of strata (1/k of total) and set as test set. There rest are train set
        dt$sets[dt[, strata_name] %in% sample(unique(dt[, strata_name]),
                                              length(unique(dt[, strata_name]))/k)] <- "test"
        
        # Refit model with training set only
        reg = NULL
        tryCatch({
          reg <- modUpdate(mod, data = filter(dt, sets == 'train'))
        }, error=function(e){
          cat("ERROR :",conditionMessage(e), "\n")}
        )
        
        if(!is.null(reg)) flag = flag +1
        
        print(paste("while loop flag = ", flag))
      }
      
      # Predict on test set
      dtest <- droplevels(subset(dt, sets == "test"))
      coefs <- fixef(reg)$cond
      if("dev:slope" %in% names(coefs)) dtest$`dev:slope` <- dtest$dev * dtest$slope
      covs <- dtest %>% select(names(coefs)) %>% as.matrix
      dtest$predall <- exp(covs %*% coefs)
      
      if (jitter) {
        if (reproducible)
          set.seed(i)
        dtest$predall <- jitter(dtest$predall)
      }
      samplepred <- function(df) {
        nrand <- sum(df$case == 0)
        obs <- rank(df$predall)[df$case == 1]
        if (length(obs)==0) obs <- NA
        if (reproducible)
          set.seed(i)
        rand <- sample(rank(df$predall[df$case == 0]),
                       1)
        return(data.frame(obs = obs, rand = rand, nrand = nrand))
      }
      ranks <- do.call(rbind, by(dtest, dtest[, strata_name], samplepred))
      nrand <- unique(ranks$nrand)
      if (length(nrand) != 1) {
        nrand <- max(nrand)
        warn[i] <- 1
      }
      kfold[i] <- cor(1:(nrand+1), table(factor(ranks$obs,
                                                levels = 1:(nrand+1))), method = "spearman")
      rd[i] <- cor(1:(nrand), table(factor(ranks$rand, levels = 1:(nrand))),
                   method = "spearman")
      
      # report out
      print(paste("kfold =", kfold[i], " ... ", "rand =", rd[i]))
      
      if (details)
        dbg[[i]] <- ranks
    }
    res <- data.frame(kfold = c(kfold, rd), type = rep(c("obs",
                                                         "rand"), each = nrepet))
    if (details)
      attr(res, "details") <- dbg
    if (sum(warn) > 0)
      warning(paste("The number of controls was not equal among stratas for",
                    sum(warn), "repetitions. Correlations might be biased.",
                    "Use 'details = TRUE' to get more details."))
    print(paste("end:", Sys.time()))
    return(res)
  }
