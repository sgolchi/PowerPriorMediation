#' Generalized Mahalanobis Distance (GMD) weights
#'
#' Returns weights computed based on the GMD and truncated at 95% quantile of the current study weight distribution
#'
#' @param cat.v matrix of categorical variables, columns represent variables
#' @param cont.v matrix of continuous variables
#' @param s.ind study index (vector); 1 for current, 0 for historical
#' @return vector of weights
#'
#' @export
# Changing the value of p in line 62 changes the cut-off value for including observations in the power prior
GMDw = function(cat.v, cont.v, s.ind) {
  d = ncol(cont.v)
  if (is.null(dim(cat.v))) z = sapply(cat.v, paste, collapse = '') else z = apply(cat.v, 1, paste, collapse = '')
  z = as.numeric(as.factor(z))
  z1 = z[s.ind==1]
  z2 = z[s.ind==0]
  lev = sort(unique(z1))
  if(is.null(dim(cont.v))) {
    cont.v1 = cont.v[s.ind==1]
    mu = mean(cont.v1[z1==1])
    Sigma = var(cont.v1[z1==1])
  } else {
    cont.v1 = cont.v[s.ind==1,]
    if (!is.null(dim(cont.v1[z1==1,]))) {
      mu = apply(cont.v1[z1==1,], 2, mean)
      Sigma = cov(cont.v1[z1==1,])
    } else {
      mu = cont.v1[z1==1,]
      Sigma =  matrix(NA, ncol(cont.v1), ncol(cont.v1))
    }
  }
  for (i in 2:length(lev)) {
    if(is.null(dim(cont.v1))) {
      mu = c(mu, mean(cont.v1[z1==i]))
      Sigma = c(Sigma, var(cont.v1[z1==i]))
    } else {
      if (!is.null(dim(cont.v1[z1==i,]))) {
        mu = rbind(mu, apply(cont.v1[z1==i,], 2, mean))
        Sigma = Matrix::bdiag(Sigma, cov(cont.v1[z1==i,]))
      } else {
        mu = rbind(mu, cont.v1[z1==i,])
        Sigma = Matrix::bdiag(Sigma, matrix(NA, ncol(cont.v1), ncol(cont.v1)))
      }
    }}
  pc = ph = c()
  J2 = c()
  for (j in lev) {
    pc[j] = mean(z1==j)
    ph[j] = mean(z2==j)
    if(is.null(dim(cont.v))) {
      J2 = c(J2, mean(ph[j], pc[j])*(cont.v[z==j] - mu[j])^2/Sigma[j])
    } else {
      J2 = c(J2, mean(ph[j], pc[j])*apply(cont.v[z==j,], 1, mahalanobis, center = mu[j,], cov = matrix(Sigma[(d*(j-1)+1):(d*j), (d*(j-1)+1):(d*j)], d, d)))
    }
    }
  pc[is.na(pc)] = 0
  ph[is.na(ph)] = 0
  J1 = sum((ph - pc)*log(ph/pc))
  d = J1 + J2
  w2 = 1- (d - min(d))/(max(d) - min(d))
  eps = quantile(w2[s.ind==1], p = .05) 
  w2[w2<eps] =0
  w2[s.ind==1] = 1
  return(w2)
}


