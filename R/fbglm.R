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
#' @param w A vector of weights (optional argument).
#' @return A list of log-likelihood, estimated coefficients, and maximum likelihood estimation results.
#'
#'@references
#'Lee, J and Breece, C. (2025) Fractional binomial regression model for count data with excess zeros.\url{https://arxiv.org/html/2410.08488v1}
#'
#'@importFrom stats C dnbinom dpois pnorm sd setNames
#'
#' @export
#'
#' @examples
#' library(agridat)
#' sample<-sample(270, 30)
#' my_y<-ridout.appleshoots$roots[sample]
#' my_x<-data.frame(pho=ridout.appleshoots$pho[sample])
#' fbglm(y=my_y, x=my_x  )
#' fbglm(y=my_y, x=my_x , w=rep(1, 30) )
#'
fbglm<-function( y, x , w=NULL){
  dfrbinom<-function(x, size, prob, h, c, start = FALSE) {
    if( size==1){dbinom(x, size, prob)} else{
      if (1 < prob || prob < 0) {
        stop("Invalid value for probability parameter.")
      }
      
      if (1 < h || h < 0) {
        stop("Invalid value for h parameter.")
      }
      if (.5*(-2*prob+2^(2*h-2)+(4*prob-prob*2^(2*h)+2^(4*h-4))^(1/2))<0){
        c<-0
      }else if (min(.5*(-2*prob+2^(2*h-2)+(4*prob-prob*2^(2*h)+2^(4*h-4))^(1/2)),1-prob) < c
                || c < 0) {
        stop("Inavlid value for c parameter.")
      }
      max2<-size-1
      
      p.fun<-function(X){
        p<-X[1]
        H<-X[2]
        c<-X[3]
        P<-rep(0,size)
        P[1]<-p+c
        d<-0
        for(i in 2:size){ d<-0; i1<-i-1; for(j in 1:i1){d<-d+(p+c*(i-j)^(2*H-2))*P[j]}
        P[i]<-p+c*i^(2*H-2)-d }; P;   }
      
      p.fun_0<-function(X){
        p<-X[1]
        H<-X[2]
        c<-X[3]
        P<-rep(0,size)
        P[1]<-p
        d<-0
        for(i in 2:size){ d<-0; i1<-i-1; for(j in 1:i1){d<-d+(p+c*(i-j)^(2*H-2))*P[j]}
        P[i]<-p-d }; P;}
      
      theo.p<-theo.p_0<-NULL
      PPm<-matrix(0, ncol=size, nrow=size); PP1<-c()
      
      yes_in<-(x %in% seq(0, size,1))
      final<-rep(0,length(x))
      if(start==FALSE){ theo.p_00<-p.fun_0(c(prob, h, c))
      theo.p<-p.fun(c(prob, h, c))
      theo.p_0<-1-cumsum(theo.p); PPm[1,]<-theo.p_00
      for(k in 1:max2){  k1<-k+1; PPm[2:k1,k1]<-PPm[1:k,1:k]%*%rev(theo.p[1:k])  }
      PP1<-PPm[,size]+c(PPm[1:max2, 1:max2]%*%rev(theo.p_0[1:max2]),0)
      final[yes_in]<-c(1-sum(theo.p_00),PP1)[(x + 1)[yes_in]]
      final
      }
      else{
        theo.p<-p.fun(c(prob, h, c))
        theo.p_0<-1-cumsum(theo.p) ;PPm[1,]<-theo.p
        for(k in 1:max2){  k1<-k+1; PPm[2:k1,k1]<-PPm[1:k,1:k]%*%rev(theo.p[1:k])  }
        PP1<-PPm[,size]+c(PPm[1:max2, 1:max2]%*%rev(theo.p_0[1:max2]),0)
        final[yes_in]<-c(1-sum(theo.p),PP1)[(x + 1)[yes_in]]
        final }
    }}
  
  if(!is.null(w)){  n<-max(y)
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
    DATA<-rbind(y,p, h,c,w )
    -sum(apply(DATA,  2 , function(x){x[5]*log(dfrbinom(x[1],n,x[2],x[3],x[4] ))}))
  }
  bbmle::parnames(likelihoodfunction2) <- c(P,H,C)
  
  est<-bbmle::mle2(
    minuslogl = likelihoodfunction2 ,
    start = setNames(rep(0, nn2),c(P,H,C) ), 
    vecpar = TRUE)
  my_list <- list( "AIC"=AIC(est),"Log-likelihood"= -est@details$value, "coef" = est@coef, "summary" = bbmle::summary(est))
  return(my_list)}

  
  else{
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
    -sum(apply(DATA,  2 , function(x){log(dfrbinom(x[1],n,x[2],x[3],x[4] ))}))
  }
  bbmle::parnames(likelihoodfunction2) <- c(P,H,C)
  
  est<-bbmle::mle2(
    minuslogl = likelihoodfunction2 ,
    start = setNames(rep(0, nn2),c(P,H,C) ), 
    vecpar = TRUE)
  my_list <- list( "AIC"=AIC(est),"Log-likelihood"= -est@details$value, "coef" = est@coef, "summary" = bbmle::summary(est))
  return(my_list)} }
