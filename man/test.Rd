% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test.R
\name{test}
\alias{test}
\title{Vuong closeness test for zero-inflated models}
\usage{
test(y, x, model1, model2)
}
\arguments{
\item{y}{A response vector.}

\item{x}{A data frame with covariates.}

\item{model1}{A character; one of "ZINB", "ZIP", "ZINB2", and "fbglm".}

\item{model2}{A character; one of "ZINB", "ZIP", "ZINB2", and "fbglm".}
}
\value{
One-sided p-value will be returned.
}
\description{
Compare zero-inflated regression models via Vuong closeness test.
}
\details{
Perform one-tailed Vuong closeness test with the null hypothesis that the two models are equally close to the true data generating process, against the alternative that one model 1 is closer than model 2.
Choose \code{model1} and \code{model2} from zero-inflated negative binomial regression ("ZINB"), extended zero-inflated negative
binomial regression ("ZINB2"), zero-inflated Poisson regression ("ZIP"), and fractional binomial regression ("fbglm").
For "ZINB2" and "fbglm", see "fbglm::ZINB2" and "fbglm::fbglm" for details.
In "ZIP" and "ZINB", all the covariates are used as regressors in both the count and zero-inflation component.
}
\examples{
library(agridat)
library(bbmle)
data <-ridout.appleshoots
my_y<-data$roots
my_x<-data.frame(pho=data$pho)
test( y=my_y, x=my_x , "fbglm", "ZINB2" )

}
\references{
Vuong, Quang H. (1989). Likelihood Ratio Tests for Model Selection and non-nested Hypotheses. Econometrica. 57 (2): 307–333.
}
