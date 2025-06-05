#' Extended zero-inflated negative binomial regression
#'@description
#'Fit extended zero-inflated negative binomial regression model via maximum likelihood.
#'
#'@details
#'The model regresses all the parameters-- zero-inflation component \eqn{\pi} (with logit link), and both the mean \eqn{\mu}
#'and dispersion parameter \eqn{\theta} (with log link)-- on covariates.
#'
#'@references
#' Lee, J. and Breece, C. (2025) Fractional binomial regression model for count data with excess zeros.\url{https://arxiv.org/html/2410.08488v1}
#'
#' @param y A response vector.
#' @param x A data frame with covariates.
#' @param w A vector of weights (optional argument). 
#' @return A list of log-likelihood, estimated coefficients, and maximum likelihood estimation results.
#'
#'
#'@importFrom stats C dnbinom dpois pnorm sd setNames
#'
#' @export
#'
#' @examples
#' Example 1
#' library(agridat)
#' sample<-sample(270, 30)
#' my_y<-ridout.appleshoots$roots[sample]
#' my_x<-data.frame(pho=ridout.appleshoots$pho[sample])
#' ZINB2(y=my_y, x=my_x  )
#' Example 2
#'ZINB2(y=c(0,1,4,2,2,1,3,0,0,6,4), x=data.frame(x1=c(1,3,4,5,4,3,6,7,4,3,1 )))
#' Example 3
#'ZINB2(y=c(0,1,4,2,2,1,3,0,0,6,4), x=data.frame(x1=c(1,3,4,5,4,3,6,7,4,3,1 )), w=c(1,1,1,1,1,1,1,1,1,1,1))
#'
ZINB2<-function(y,x,w=NULL){
  if(!is.null(w)){ 
  w_0<-which(y==0)
  dnb<-function(k,theta,mu){prob<-theta/(mu+theta);dnbinom(x=k, size=theta, prob=prob, log = FALSE)}
  
  data3<-cbind(intercept=rep(1, length(y)),x)
  nn<-dim(data3)[2]
  nn0<-nn+1
  nn10<-2*nn
  nn1<-2*nn+1
  nn2<-3*nn
  
  Theta<-paste0(colnames( data3),"(theta)")
  Mu<-paste0(colnames( data3),"(mu)")
  Pi<-paste0(colnames( data3),"(pi)")
  
  likelihood_nb<-function( X){
    theta<-exp(rowSums( mapply(`*`,data3,X[1:nn]) ))
    mu<-exp(rowSums( mapply(`*`,data3,X[nn0:nn10])   ))
    pi<-1/(1+exp(-(rowSums( mapply(`*`,data3,X[nn1:nn2]) ) )))
    DAT<-rbind(theta, mu)
    DAT2<-rbind(y,theta, mu)[,-w_0]
    
    -sum(w[w_0]*log(pi[w_0]+(1-pi[w_0])* apply(DAT[, w_0]  ,2,
                                               function(x){dnb(0,x[1],x[2] )} )))-
      sum(w[-w_0]*log((1- pi[-w_0])* apply(DAT2  ,2,
                                           function(x){dnb(x[1],x[2],x[3] )} )  ))
  }
  
  bbmle::parnames(likelihood_nb) <- c(Theta, Mu, Pi)
  
  est<-bbmle::mle2(
    minuslogl = likelihood_nb ,
    start = setNames(rep(0, nn2),c(Theta, Mu, Pi ) ),
    vecpar = TRUE)
  my_list <- list( "AIC"=AIC(est),"Log-likelihood"= -est@details$value, "coef" = est@coef, "summary" = bbmle::summary(est))
  return(my_list)}  
  else{  w_0<-which(y==0)
  dnb<-function(k,theta,mu){prob<-theta/(mu+theta);dnbinom(x=k, size=theta, prob=prob, log = FALSE)}
  
  data3<-cbind(intercept=rep(1, length(y)),x)
  nn<-dim(data3)[2]
  nn0<-nn+1
  nn10<-2*nn
  nn1<-2*nn+1
  nn2<-3*nn
  
  Theta<-paste0(colnames( data3),"(theta)")
  Mu<-paste0(colnames( data3),"(mu)")
  Pi<-paste0(colnames( data3),"(pi)")
  
  likelihood_nb<-function( X){
    theta<-exp(rowSums( mapply(`*`,data3,X[1:nn]) ))
    mu<-exp(rowSums( mapply(`*`,data3,X[nn0:nn10])   ))
    pi<-1/(1+exp(-(rowSums( mapply(`*`,data3,X[nn1:nn2]) ) )))
    DAT<-rbind(theta, mu)
    DAT2<-rbind(y,theta, mu)[,-w_0]
    
    -sum(log(pi[w_0]+(1-pi[w_0])* apply(DAT[, w_0]  ,2,
                                        function(x){dnb(0,x[1],x[2] )} )))-
      sum(log((1- pi[-w_0])* apply(DAT2  ,2,
                                   function(x){dnb(x[1],x[2],x[3] )} )  ))
  }
  
  bbmle::parnames(likelihood_nb) <- c(Theta, Mu, Pi)
  
  est<-bbmle::mle2(
    minuslogl = likelihood_nb ,
    start = setNames(rep(0, nn2),c(Theta, Mu, Pi ) ),
    vecpar = TRUE)
  my_list <- list( "AIC"=AIC(est),"Log-likelihood"= -est@details$value, "coef" = est@coef, "summary" = bbmle::summary(est))
  return(my_list)}
  
}


