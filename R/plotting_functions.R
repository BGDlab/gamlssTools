
################################################

# misc. plotting functions

################################################

#' Worm Plot (Taki's Version)
#' 
#' Output more robust worm plots for assessing gamlss fit in ggplot
#' 
#' This function fixes an inconsistency in the original package - that if do a worm plot with a formula 
#' it gives the wrong answer compared to if you provide a vector of values to plot against. This is 
#' rectified by adding a check that the partition of the covariate is actually providing a partition.
#' The percentage of points that fall outside the pointwise confidence intervals are shown in blue
#' (5% indicates good fit). Written and contributed by the illustrious Taki Shinohara :)
#' 
#' @param object gamlss model object
#' @param xvar vector containing values of predictor used for plotting model (requires `resid`)
#' @param resid vector containing the residuals of the same model object (requires `xvar`)
#' @param n.inter (optional) number of subsets with ~equal number of points to plot across the range of `xvar`. Defaults to 4.
#' @param xlim.worm control plot range(s)
#' @param ylim.worm control plot range(s)
#' 
#' @returns list with ggplot objects of worm plot(s) ($plot) and df of pts outside dotted CI ($outliers)
#' 
#' @examples
#' iris_model <- gamlss(formula = Sepal.Width ~ Sepal.Length + Species, sigma.formula = ~ Sepal.Length, data=iris, family=BCCG)
#' 
#' #get one worm plot by just passing the model object
#' wp.taki(iris_model)
#'
#' #or get subplots across the range of a covariate
#' wp.taki(xvar=iris$Sepal.Length, resid=iris_model$residuals)
#' wp.taki(xvar=iris$Sepal.Length, resid=resid(iris_model))
#' wp.taki(xvar=iris$Sepal.Length, resid=resid(iris_model), n.inter=10)
#' 
#' @export
wp.taki<-function (object = NULL, xvar = NULL, resid = NULL, n.inter = 4, 
                   xlim.worm = 3.5, ylim.worm = 12 * sqrt(n.inter/length(resid)),
                   #show.given = TRUE, line = TRUE, ylim.all = 12 * sqrt(1/length(resid)), xlim.all = 4,
                   #pch = 21, bg = "wheat", col = "red", bar.bg = c(num = "light blue"),
                   #cex = 1, cex.lab = 1, xcut.points = NULL,
                   ...) 
{
  
  ### ggplotification from 4-5/23
  
  ## Functions for worm plotting:
  # This function estimates the worm plot function from some residuals
  get.wp.x<-function(y.in) {
    qq <- as.data.frame(qqnorm(y.in, plot = FALSE))
    return(qq$x)
  }
  
  # This function estimates the worm plot function from some residuals
  get.wp.y<-function(y.in) {
    qq <- as.data.frame(qqnorm(y.in, plot = FALSE))
    return(qq$y - qq$x)
  }
  
  # find the limits based on the normal distribution
  get.lims<-function(resid) {
    lims.df <- data.frame(zval=seq(-xlim.worm, xlim.worm, 0.1))
    lims.df$p <- pnorm(lims.df$zval)
    lims.df$se <- (1/dnorm(lims.df$zval)) * (sqrt(lims.df$p * (1 - lims.df$p)/length(resid)))
    lims.df$low <- qnorm((1 - 0.95)/2) * lims.df$se; lims.df$high <- -lims.df$low
    return(lims.df)
  }
  
  
  ## Helpful function for dealing with formulas (from original wp())
  deparen <- function(expr) {
    while (is.language(expr) && !is.name(expr) && deparse(expr[[1L]])[1L] == 
           "(") expr <- expr[[2L]]
    expr
  }
  
  
  ## Make sure there's either an object or a set of residuals
  if (is.null(object) && is.null(resid)) 
    stop(paste("A fitted object with resid() method or the argument resid should be used ", 
               "\n", ""))
  
  ## If there's an object, capture its residuals
  resid <- if (is.null(object)) { resid } else { resid(object) }
  
  ## Check if there's a dataframe in the object, parse it
  DataExist <- FALSE
  if (!is.null(object) && any(grepl("data", names(object$call))) && 
      !(object$call["data"] == "sys.parent()()")) {
    DaTa <- eval(object$call[["data"]])
    DataExist <- TRUE
  }
  
  ## If no covariate values are provided, just look at normality of the residuals
  if (is.null(xvar)) {
    
    # find the limits based on the normal distribution
    lims.df <- data.frame(zval=seq(-xlim.worm, xlim.worm, 0.1))
    lims.df$p <- pnorm(lims.df$zval)
    lims.df$se <- (1/dnorm(lims.df$zval)) * (sqrt(lims.df$p * (1 - lims.df$p)/length(resid)))
    lims.df$low <- qnorm((1 - 0.95)/2) * lims.df$se; lims.df$high <- -lims.df$low
    
    #Construct the worm plot dataframe
    wp.df <- data.frame(y = resid %>% get.wp.y, x = resid %>% get.wp.x)
    
    #Count outer points
   wp.dt <- wp.df %>% arrange(x) %>% as.data.table()
   lims.dt <- as.data.table(lims.df)
   combo.dt <- lims.dt[wp.dt, on = .(zval == x), roll=TRUE]
   if (sum(is.na(combo.dt))>0){
     warning("missing some CI values, try increasing xlim.worm")
   }
   n_outer <- combo.dt %>%
     mutate(outer = ifelse((y < low | y > high), 1, 0)) %>%
     summarise(n = n(),
               n_outer = sum(outer)) %>%
     mutate(pcnt = n_outer/n,
            x = xlim.worm * 0.95,
            y = ylim.worm * 0.95)
    
    #Return the plot
    p <- ggplot(wp.df,aes(x=x,y=y)) + geom_smooth(method=lm,formula=y~poly(x,3)) + 
      geom_point() + theme_classic() +
      xlab("Unit Normal Quantile") + ylab("Deviation") + 
      {if (is.finite(xlim.worm)) xlim(c(-xlim.worm, xlim.worm))} + 
      #{if (is.finite(ylim.worm)) { ylim(c(-ylim.worm, ylim.worm)) } else { ylim(c(-1,1))} } +
      {if (is.finite(ylim.worm)) { coord_cartesian(ylim = c(-ylim.worm, ylim.worm)) } else { coord_cartesian(ylim = c(-1,1))} } + 
      geom_line(data = lims.df, aes(x=zval,y=low),linetype = "dashed") +
      geom_line(data = lims.df, aes(x=zval,y=high),linetype = "dashed")  + 
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_text(data=n_outer, 
                mapping = aes(x = x, 
                              y = y, 
                              label = scales::percent(pcnt)),
                inherit.aes = FALSE,
                #label.size = 0.15,
                color="blue")
    
  } else {
    
    ## If a covariate is provided as a vector ...
    if (!is(xvar, "formula")) {
      if (length(resid) != length(xvar)) stop("Error - incorrect length of predictor vector...")
      
      if (is.factor(xvar)) {
        wp.df<-data.frame(z=xvar)
      } else { 
        wp.df<-data.frame(z=cut(xvar,breaks=quantile(xvar,probs=seq(0,1,length.out=n.inter+1)),include.lowest=TRUE))
      }
      
    }
    
    ## If a covariate is provided as a formula ... 
    if (is(xvar, "formula")) {
      
      if (DataExist) {
        
        # extract the facet variable as the first variable in the provided formula
        ## note - if we want to do this by multiple variables, this will require a bit more coding
        xvar.vec<-eval(deparen(deparen(xvar)[[2L]]), envir = as.environment(DaTa))
        if (is.factor(xvar.vec)) {
          wp.df<-data.frame(z=xvar.vec)
        } else { 
          wp.df<-data.frame(z=cut(xvar.vec,breaks=quantile(xvar.vec,probs=seq(0,1,length.out=n.inter+1)),include.lowest=TRUE))
        }
        
      } else { 
        stop("Dataframe missing, exiting...")
      }
    }
    
    #Construct data frame for main plot
    wp.df$resid<-resid
    wp.df <- wp.df %>% group_by(z) %>% mutate(y = resid %>% get.wp.y, x = resid %>% get.wp.x)
    
    #Construct data frame for confidence limits
    lims.df<-data.frame(z=wp.df$z,resid=wp.df$resid)
    lims.df<-lims.df %>% group_by(z) %>% reframe(lims.df=get.lims(resid)) %>% tidyr::unnest(cols=c(lims.df))
    
    #Count outer points
    wp.dt <- wp.df %>% group_by(z) %>% arrange(x, .by_group = TRUE) %>% as.data.table()
    lims.dt <- as.data.table(lims.df)
    combo.dt <- lims.dt[wp.dt, on = .(z ==z, zval == x), roll=TRUE]
    if (sum(is.na(combo.dt))>0){
      warning("missing some CI values, try increasing xlim.worm")
    }
    n_outer <- combo.dt %>%
      mutate(outer = ifelse((y < low | y > high), 1, 0)) %>%
      group_by(z) %>%
      summarise(n = n(),
                n_outer = sum(outer)) %>%
      mutate(pcnt = n_outer/n)
    
   #Return the plot
    p <- ggplot(wp.df,aes(x=x,y=y)) + geom_smooth(method=lm,formula=y~poly(x,3)) + 
      geom_point() + facet_wrap(~z) + theme_classic() +
      xlab("Unit Normal Quantile") + ylab("Deviation") +
      {if (is.finite(xlim.worm)) xlim(c(-xlim.worm, xlim.worm))} + 
      {if (is.finite(ylim.worm)) { ylim(c(-ylim.worm, ylim.worm)) } else { ylim(c(-1,1))} } +
      geom_line(data = lims.df, aes(x=zval,y=low),linetype = "dashed") +
      geom_line(data = lims.df, aes(x=zval,y=high),linetype = "dashed") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_text(data=n_outer, 
                mapping = aes(x = Inf, y = Inf, 
                label = scales::percent(pcnt)),
                hjust = 1.5,
                vjust = 1.5,
                label.size = 0.15,
                color="blue")
  }
  out <- list()
  out$plot <- p
  out$outliers <- n_outer
  return(out)
}