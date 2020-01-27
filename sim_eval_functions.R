#Functions for evaluation of simulations

#Function to calcluate MSE
mse = function(x,y) mean((x-y)^2)
rmse = function(x,y) sqrt(mean((x-y)^2))
bias = function(x,y) mean(abs(x-y))

mse_by_el = function(x, y) {
  m = matrix(0, NROW(x), NCOL(x))
  for (i in 1:NROW(x)) {
    for (j in 1:NCOL(x)) {
      m[i,j] = as.numeric(mse(x[i,j], y[i,j]))
    }
  }
  return(m)
}

rmse_by_el = function(x,y) {
  m = mse_by_el(x,y)
  return(sqrt(m))
}

abs_pct_error = function(x, y) {
  mean(abs((x-y)/y))
}

#Function to calculate MSE on lower triangular part of two matrices (diag included)
mse_lower = function(x,y) {
  xl = x[lower.tri(x, diag=TRUE)]
  yl = y[lower.tri(y, diag=TRUE)]
  mean((xl-yl)^2)
}

rmse_lower = function(x,y) {
  sqrt(mse_lower(x,y))
  #abs(m/mean(y[lower.tri(y, diag=TRUE)]))
}

#Function to calculate bias on lower triangular part of two matrices (diag included, second value is true parameter value)
bias_lower = function(x,y) {
  xl = x[lower.tri(x, diag=TRUE)]
  yl = y[lower.tri(y, diag=TRUE)]
  mean(abs(xl-yl))
}


#Function to calculate MSE on lower triangular part of two matrices (no diag included)
mse_lower_nodiag = function(x,y) { 
  xl = x[lower.tri(x)]
  yl = y[lower.tri(y)]
  mean((xl-yl)^2)
}

rmse_lower_nodiag = function(x,y) {
  sqrt(mse_lower_nodiag(x,y))
  #abs(m/mean(y[lower.tri(y)]))
}

#Function to calculate bias on lower triangular part of two matrices (no diag included, second value is true parameter value)
bias_lower_nodiag = function(x,y) { 
  xl = x[lower.tri(x)]
  yl = y[lower.tri(y)]
  mean(abs(xl-yl))
}

#Calculate mean absolute bias for selected elements of matrix
bias_select = function(x, y, select) {
  #return(mean(abs(x[select]-y[select])))
  return(mean(abs(x[select])-abs(y[select])))
}

#Function to calculate MSE on only diagonal part
mse_diag = function(x,y) {
  xd = diag(x)
  yd = diag(y)
  mean((xd-yd)^2)
}

rmse_diag = function(x,y) {
  sqrt(mse_diag(x,y))
  #abs(m/mean(diag(y)))
}

#Function to calculate bias on only diagonal part (second value is true parameter value)
bias_diag = function(x,y) {
  xd = diag(x)
  yd = diag(y)
  mean(abs(xd-yd))
}

#Function to look at mean absolute magnitude of off-diagonal elements of precision matrix
mean_abs_lower_nodiag = function(x) {
  lt = lower.tri(x)
  return(mean(abs(x[lt])))
}

#Function to look at mean absolute magnitude of diagonal elements of precision matrix
mean_abs_lower_diag = function(x) {
  return(mean(abs(diag(x))))
}

#Convert precision matrix to (unweighted) adjacency matrix
prec_to_adj = function(x) {
  adj = abs(x) > 0
  diag(adj) = 0
  return(adj)
}

#Convert precision matrix to partial correlation matrix
prec2partial = function(x) {
  x <- -cov2cor(x)
  diag(x) <- 1
  return(x)
}

#Function that checks whether confidence intervals cross zero to create (unweighted) network adjacency matrix
ci_to_adj = function(lower, upper) {
  adj = (lower>0) | (upper<0)
  diag(adj) = 0
  return(adj)
}

num_zeros = function(x, prop=FALSE) {
  lt = lower.tri(x)
  r = ifelse(prop, mean(x[lt]==0), sum(x[lt]==0))
  return(r)
}

#Function to compare network connections based on precision matrices (only looks at edge vs. no edge)
# y contains the true precision matrix
compareNetwork = function(x,y,x.lower=NULL,x.upper=NULL) {
  if (is.null(x.lower) & is.null(x.upper)) {
    x.adj = prec_to_adj(x)
  } else if (!is.null(x.lower) & !is.null(x.upper)) {
    x.adj = ci_to_adj(x.lower, x.upper)
  } else {
    stop("x.lower and x.upper must both be NULL or both be non-NULL")
  }
  
  y.adj = prec_to_adj(y)
  xl = x.adj[lower.tri(x.adj)]
  yl = y.adj[lower.tri(y.adj)]
  acc = mean(xl==yl)
  sens = sum(xl[yl==1])/sum(yl==1)
  spec = sum(xl[yl==0]==0)/sum(yl==0)
  
  rvec = c(acc,sens,spec)
  names(rvec) = c("acc","sens","spec")
  return(rvec)
}

#Comparing real prec matrix to simulated one (to get magnitude of difference)
compareMat_plot = function(x, y, lower.ci, upper.ci, diag=FALSE, ylim=NULL, ...) {
  lt = lower.tri(x, diag=diag)
  
  if(is.null(ylim)) {
      ylim = c(min(lower.ci[lt]), max(upper.ci[lt]))
  }
  plot(x[lt], y[lt], ylim=ylim, ...)
  segments(x[lt], lower.ci[lt], x[lt], upper.ci[lt])
  abline(0,1)
  abline(h=0, col="red", lty=2)
  abline(v=0, col="red", lty=2)
}

compareMat_plot_nonsym = function(x, y, lower.ci, upper.ci, ylim=NULL, ...) {
  if(is.null(ylim)) {
    ylim = c(min(lower.ci), max(upper.ci))
  }
  plot(x, y, ylim=ylim, ...)
  segments(x, lower.ci, x, upper.ci)
  abline(0,1)
  abline(h=0, col="red", lty=2)
  abline(v=0, col="red", lty=2)
}

# See if credible intervals of lower triangular part of estimated precision matrix overlaps with true values
credInt_overlap_lower = function(y, lower.ci, upper.ci, diag=FALSE) {
  lt = lower.tri(y, diag=diag)
  return(mean((y[lt] < upper.ci[lt]) & (y[lt] > lower.ci[lt])))
}

# Function to scale down matrix based on maximum abs off-diagonal entry
scale_by_max <- function(x) {
  return(x/max(abs(x)))
}

# Function to calculated weighted natural connectivity
w_natcon <- function(x) {
  if (any(is.na(x) | is.infinite(x))) {
    e <- NA
  } else {
    e <- eigen(abs(x))$values
  }
  return(log(mean(exp(e))))
}

# Function to scale by maximum value and calculated weighted natural connectivity
scale_natcon <- function(x) {
  s <- scale_by_max(x)
  return(w_natcon(s))
}


