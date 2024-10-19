#' Vuong closeness test for zero-inflated models 
#'
#'@description 
#'Compare zero-inflated regression models via Vuong closeness test.
#'
#'@details 
#'Perform one-tailed Vuong closeness test with the null hypothesis that the two models are equally close to the true data generating process, against the alternative that one model 1 is closer than model 2. 
#'Choose `model1` and `model2` from zero-inflated negative binomial regression ("ZINB"), extended zero-inflated negative
#'binomial regression ("ZINB2"), zero-inflated Poisson regression ("ZIP"), and fractional binomial regression ("fbglm").
#'For "ZINB2" and "fbglm", see "fbglm::ZINB2" and "fbglm::fbglm" for details. 
#'In "ZIP" and "ZINB", all the covariates are used as regressors in both the count and zero-inflation component.
#'
#' @param y A response vector.
#' @param x A data frame with covariates. 
#'
#' @param model1 A character; one of "ZINB", "ZIP", "ZINB2", and "fbglm".
#' @param model2 A character; one of "ZINB", "ZIP", "ZINB2", and "fbglm".
#'
#'@references
#'Vuong, Quang H. (1989). Likelihood Ratio Tests for Model Selection and non-nested Hypotheses. Econometrica. 57 (2): 307–333.
#'
#'@return One-sided p-value will be returned.
#'
#' @export
#' 
#' @example 
#' library(agridat)
#' library(bbmle)
#' data <-ridout.appleshoots
#' my_y<-data$roots
#' my_x<-data.frame(pho=data$pho)
#' test( y=my_y, x=my_x , "fbglm", "ZINB2" )
#' 
test<-function( y,x,model1 ,model2){
  compute<-function(xx, y,x )
  {if(xx=="ZINB"){    
    y<-my_y;x<-my_x
    data3<-cbind(intercept=rep(1, length(y)),x) 
    nn<-dim(data3)[2]
    nn0<-nn+1
    nn10<-nn+2
    nn1<-2*nn+1
    
    ind<-rep(0, length(y))
    ind[which(y==0)]<-1
    
    dnb<-function(k,theta,mu){prob<-theta/(mu+theta);dnbinom(x=k, size=theta, prob=prob, log = FALSE)}
    likelihood_nb_1<-function(X){
      theta<-X[1]
      mu<-exp(  rowSums( mapply(`*`,data3,X[2:nn0]) )  )
      pi<-1/(1+exp(-(  rowSums( mapply(`*`,data3,X[nn10:nn1]) )  ) ))
      
      DAT2<-rbind(y,theta, mu)
      log(    pi*ind  +(1- pi)* apply(DAT2  ,2, 
                                      function(x){dnbinom(x[1],size=x[2],mu=x[3] )} )  )
    }
    data<-data.frame(y=y, x=x)
    C.ZINB<-pscl::zeroinfl(   y~.|., data=data, dist='negbin')$coef
    input<-c(c(C$theta)[[1]],c(C.ZINB)$count, c(C.ZINB)$zero)
    
    my.list<-list( "likelihood"=likelihood_nb_1(input ) ,
                   "par"= nn1)
    
    return( my.list)
  }
    if(xx=="ZINB2"){  
      
      data3<-cbind(intercept=rep(1, length(y)),x) 
      nn<-dim(data3)[2]
      nn0<-nn+1
      nn10<-2*nn
      nn1<-2*nn+1
      nn2<-3*nn
      
      ind<-rep(0, length(dat$roots))
      ind[which(dat$roots==0)]<-1
      
      dnb<-function(k,theta,mu){prob<-theta/(mu+theta);dnbinom(x=k, size=theta, prob=prob, log = FALSE)}
      
      likelihood_nb_2<-function( X  ){
        theta<-exp(rowSums( mapply(`*`,data3,X[1:nn]) )  )
        mu<-exp(rowSums( mapply(`*`,data3,X[nn0:nn10]) )  )
        pi<-1/(1+exp(-(  rowSums( mapply(`*`,data3,X[nn1:nn2]) ) ) ))
        
        DAT2<-rbind(y,theta, mu)
        log(    pi*ind  +(1- pi)* apply(DAT2  ,2, 
                                        function(x){dnbinom(x[1],size=x[2],mu=x[3] )} )  )
      }
      M.ZINB2<-ZINB2(y,x)
      C.ZINB2<-M.ZINB2$coef
      my.list<-list( "likelihood"=likelihood_nb_2( C.ZINB2),
                     "par"= nn2)
      return(my.list )}
    if(xx=="fbglm"){
      n<-max(y)
      data3<-cbind(intercept=rep(1, length(y)),x) 
      nn<-dim(data3)[2]
      nn0<-nn+1
      nn10<-2*nn
      nn1<-2*nn+1
      nn2<-3*nn
      
      likelihoodfunction3<-function(X){ 
        p<-1/(1+exp(-( rowSums( mapply(`*`,data3,X[1:nn]))  )));
        h<-1/(1+exp(-( rowSums( mapply(`*`,data3,X[nn0:nn10])) )));
        c00<- apply(rbind(p,h) ,2, function(x){min( .5*(-2*x[1]+2^(2*x[2]-2)+(4*x[1]-x[1]*2^(2*x[2])+2^(4*x[2]-4))^(1/2)),1-x[1]  ) })
        c<-c00/(1+exp(-( rowSums( mapply(`*`,data3,X[nn1:nn2])))));
        
        DATA<-rbind(y,p, h,c )
        (apply(DATA,  2 , function(x){log(frbinom::dfrbinom(x[1],n,x[2],x[3],x[4] ))}))
      }
      FB<-fbglm(y,x)
      C.FB<-FB$coef
      my.list<-list( "likelihood"=likelihoodfunction3(C.FB) ,
                     "par"= nn2)
      
      return( my.list  )  }
    if(xx=="ZIP"){ 
      x<-my_x
      y<-my_y
      ind<-rep(0, length(y))
      ind[which(dat$roots==0)]<-1
      data3<-cbind(intercept=rep(1, length(y)),x) 
      nn<-dim(data3)[2]
      nn0<-nn+1
      nn10<-2*nn
      
      likelihood_p<-function(  X ){
        la<-exp(rowSums( mapply(`*`,data3,X[1:nn])) )
        pi<-1/(1+exp(-(rowSums( mapply(`*`,data3,X[nn0:nn10])) ) ))
        
        DAT2<-rbind(y,la)
        log(    pi*ind  +(1- pi)* apply(DAT2  ,2, 
                                        function(x){dpois(x[1],lambda=x[2] )} )  )
      }
      
      data<-data.frame(y=y, x=x)
      C.P<-pscl::zeroinfl(   y~.|., data=data, dist='poisson')$coef
      input2<-c(c(C.P$count), c(C.P$zero))
      
      my.list<-list( "likelihood"=likelihood_p( input2) ,
                     "par"= nn10)
      
      return( my.list    )}}
  
  
  aa<- compute(model1, y,x)
  bb<-compute(model2, y,x)
  a1<-aa$likelihood
  b1<-bb$likelihood
  final<-( sum(a1)-sum(b1)- (aa$par -bb$par)*log(length(y))/2)/(sqrt(length(y))*sd(a1-b1))
  final.p<- 1-pnorm(final,0,1)
  return(final.p)
}

