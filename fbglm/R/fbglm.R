#' Fractional binomial regression model
#'@description
#'Fit a fractional binomial regression model via maximum likelihood.
#'
#'@details
#'Fractional binomial distribution can be considered as zero-inflated, over-dispersed binomial model, and it has three parameters
#' \eqn{(p,H,c)} in addition to the number of trials \eqn{n}.
#' We use a specific parametrization such that  \eqn{p,H,c \in (0,1)}, and
#' regress these parameters with logit link on the covariates, while letting \eqn{n} as the maximum of the response `y`.
#'
#' @param y A response vector.
#' @param x A data frame with covariates. 
#'
#' @return A list of log-likelihood, estimated coefficients, and maximum likelihood estimation results. 
#'
#'@references
#'Breece, C. and Lee, J. (2024) Fractional binomial regression model for count data with excess zeros.\url{https://arxiv.org/html/2410.08488v1}
#'
#'
#'
#' @export
#' 
#' @examples
#' library(agridat)
#' library(bbmle)
#' data <-ridout.appleshoots
#' my_y<-data$roots
#' my_x<-data.frame(pho=data$pho)
#' fbglm(y=my_y, x=my_x  )
#' 
fbglm<-function( y, x ){
    n<-max(y)
data3<-cbind(intercept=rep(1, length(y)),x) 
nn<-dim(data3)[2]
nn0<-nn+1
nn10<-2*nn
nn1<-2*nn+1
nn2<-3*nn

P<-paste0(colnames( data3),"(p)")
H<-paste0(colnames( data3),"(h)")
C<-paste0(colnames( data3),"(c)")

X<-matrix(nrow=1, ncol=3*dim(data3)[2])

likelihoodfunction2<-function(X){ 
  p<-1/(1+exp(-(rowSums( mapply(`*`,data3,X[1:nn]))
  )) );
  h<-1/(1+exp(-(rowSums( mapply(`*`,data3,X[nn0:nn10] ))
  )) );
  c00<- apply(rbind(p,h) ,2, function(x){min( .5*(-2*x[1]+2^(2*x[2]-2)+(4*x[1]-x[1]*2^(2*x[2])+2^(4*x[2]-4))^(1/2)),1-x[1]  ) })
  c<-c00/(1+exp(-(rowSums( mapply(`*`,data3, X[nn1:nn2] ))
  )) );
  DATA<-rbind(y,p, h,c )
  -sum(apply(DATA,  2 , function(x){log(frbinom::dfrbinom(x[1],n,x[2],x[3],x[4] ))}))
}
parnames(likelihoodfunction2) <- c(P,H,C)

est<-bbmle::mle2(
  minuslogl = likelihoodfunction2 ,
  start = setNames(rep(0, nn2),c(P,H,C) ),
  vecpar = TRUE)
my_list <- list( "Log-likelihood"= -est@details$value, "coef" = est@coef, "summary" = summary(est))
return(my_list)  }


  

