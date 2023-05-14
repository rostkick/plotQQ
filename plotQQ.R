library(ggplot2)


lambdaGC <- function (data, plot = FALSE, hline=-log10(0.05/2500),
                      proportion = 1, method = "regression", filter = TRUE, df = 1, exp=NA, ...){
  #` Reworked GenABEL::estlambda(). Only plot function was changing.
  pvals <- copy(data)
  len_pvals <- length(pvals)
  data <- data[which(!is.na(data))]
  if (proportion > 1 || proportion <= 0) 
    stop("proportion argument should be greater then zero and less than or equal to one")
  ntp <- round(proportion * length(data))
  if (ntp < 1) 
    stop("no valid measurements")
  if (ntp == 1) {
    warning(paste("One measurement, lambda = 1 returned"))
    return(list(estimate = 1, se = 999.99))
  }
  if (ntp < 10)
    warning(paste("number of points is too small:", ntp))
  if (min(data) < 0) 
    stop("data argument has values <0")
  if (max(data) <= 1) {
    data <- qchisq(data, 1, lower.tail = FALSE)
  }
  if (filter) {
    data[which(abs(data) < 1e-08)] <- NA
  }
  data <- sort(data)
  ppoi <- ppoints(data)
  ppoi <- sort(qchisq(ppoi, df = df, lower.tail = FALSE))
  data <- data[1:ntp]
  ppoi <- ppoi[1:ntp]
  out <- list()
  
  if(is.na(exp)) exp=-log10(1:len_pvals/len_pvals)
  
  dd <- data_frame(exp=exp, obs=-log10(sort(pvals)))
  if (method == "regression") {
    s <- summary(lm(data ~ 0 + ppoi))$coeff
    out$estimate <- s[1, 1]
    out$se <- s[1, 2]
  }
  else if (method == "median") {
    out$estimate <- median(data, na.rm = TRUE)/qchisq(0.5, df)
    out$se <- NA
  }
  else {
    stop("'method' should be either 'regression' or 'median'!")
  }

  clower = -log10(qbeta(
    p = (1 - 0.95) / 2,
    shape1 = 1:len_pvals,
    shape2 = len_pvals:1))
  cupper   = -log10(qbeta(
    p = (1 + 0.95) / 2,
    shape1 = 1:len_pvals,
    shape2 = len_pvals:1))
  dd$clower <- clower
  dd$cupper <- cupper
  if (plot) {
    p <- ggplot(dd, aes(x=exp, y=obs))+
      geom_point(shape=19, size=3, alpha=0.8, color='#3C5488FF')+
      geom_hline(yintercept=hline, linetype=5, color='#DC0000FF')+
      geom_abline(intercept = 0, slope = 1, alpha = 0.8, color='black')+
      geom_line(data = dd, aes(exp, cupper), linetype = 2)+
      geom_line(data = dd, aes(exp, clower), linetype = 2)+
      scale_x_continuous(expression(paste("Expected -log"[10], plain(P))))+
      scale_y_continuous(expression(paste("Observed -log"[10], plain(P))))+
      theme_classic()
    show(p)
  }
  print(out)
  return(p)
}

