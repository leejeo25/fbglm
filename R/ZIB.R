#' Zero-inflated binomial regression
#'@description
#'Fit zero-inflated binomial regression model via maximum likelihood.
#'
#'@details
#'The model regresses all the parameters-- zero-inflation component \eqn{\pi} (with logit link), and the probability of success \eqn{p} (with log link)-- on covariates.
#'
#'@references
#' Lee, J, and Breece, C. (2025) Fractional binomial regression model for count data with excess zeros.\url{https://arxiv.org/html/2410.08488v1}
#'
#' @param y A response vector.
#' @param x A data frame with covariates.
#' @param w (optional) A vector of weights.  
#'
#' @return A list of AIC, log-likelihood, estimated coefficients, and maximum likelihood estimation results.
#'
#'
#'@importFrom stats C dnbinom dbinom dpois pnorm sd setNames
#'
#' @export
#'
#' @examples
#'
#' ZIB(y=c(0,1,4,2,2,1,3,0,0,6,4), x=data.frame(x1=c(1,3,4,5,4,3,6,7,4,3,1 )))
#'ZIB(y=c(0,1,4,2,2,1,3,0,0,6,4), x=data.frame(x1=c(1,3,4,5,4,3,6,7,4,3,1 )), w=c(1,1,1,1,1,1,1,1,1,1,1))
#'
ZIB<-function(y,x,w=NULL){
  if(!is.null(w)){
  w_0<-which(y==0)
n<-max(y)

data3<-cbind(intercept=rep(1, length(y)),x)
nn<-dim(data3)[2]
nn0<-nn+1
nn10<-2*nn

P<-paste0(colnames( data3),"(P)")
Pi<-paste0(colnames( data3),"(pi)")

likelihood_b<-function( X){
  P<-1/(1+exp(-(rowSums( mapply(`*`,data3,X[1:nn]) ) )))
  pi<-1/(1+exp(-(rowSums( mapply(`*`,data3,X[nn0:nn10]) ) )))
  
  DAT2<-rbind(y,P)[,-w_0]
  
  -sum(w[w_0]*log(pi[w_0]+(1-pi[w_0])*dbinom(0,n,P[w_0] )))-
    sum(w[-w_0]*log((1- pi[-w_0])* apply(DAT2  ,2,
                                         function(x){dbinom(x[1],n,x[2] )} )  ))
}

bbmle::parnames(likelihood_b) <- c(P, Pi)

est<-bbmle::mle2(
  minuslogl = likelihood_b ,
  start = setNames(rep(0, nn10),c(P, Pi ) ),
  vecpar = TRUE)
my_list <- list( "AIC"=AIC(est),"Log-likelihood"= -est@details$value, "coef" = est@coef, "summary" = bbmle::summary(est))
return(my_list) }
  else{ w_0<-which(y==0)
  n<-max(y)
  
  data3<-cbind(intercept=rep(1, length(y)),x)
  nn<-dim(data3)[2]
  nn0<-nn+1
  nn10<-2*nn
  
  P<-paste0(colnames( data3),"(P)")
  Pi<-paste0(colnames( data3),"(pi)")
  
  likelihood_b<-function( X){
    P<-1/(1+exp(-(rowSums( mapply(`*`,data3,X[1:nn]) ) )))
    pi<-1/(1+exp(-(rowSums( mapply(`*`,data3,X[nn0:nn10]) ) )))
    
    DAT2<-rbind(y,P)[,-w_0]
    
    -sum(log(pi[w_0]+(1-pi[w_0])*dbinom(0,n,P[w_0] )))-
      sum(log((1- pi[-w_0])* apply(DAT2  ,2,
                                   function(x){dbinom(x[1],n,x[2] )} )  ))
  }
  
  bbmle::parnames(likelihood_b) <- c(P, Pi)
  
  est<-bbmle::mle2(
    minuslogl = likelihood_b ,
    start = setNames(rep(0, nn10),c(P, Pi ) ),
    vecpar = TRUE)
  my_list <- list( "AIC"=AIC(est),"Log-likelihood"= -est@details$value, "coef" = est@coef, "summary" = bbmle::summary(est))
  return(my_list) 
  }
  
}

ZIB(y=c(0,1,4,2,2,1,3,0,0,6,4), x=data.frame(x1=c(1,3,4,5,4,3,6,7,4,3,1 )))

ZIB(y=c(0,1,4,2,2,1,3,0,0,6,4), x=data.frame(x1=c(1,3,4,5,4,3,6,7,4,3,1 )), w=c(1,1,1,1,1,1,1,1,1,1,1))
