\name{MD_OLS}
\alias{MD_OLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{R Package for Quadratic Functional of Moderate-Dimensional OLS}
\description{
This package conduct statistical inference for quadratic functional of the moderate-dimensional linear model in Guo and Cheng (2020).
}
\usage{
MD_OLS(type_inference, level_significance, null_value, X, Y, V, W)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type_inference}{the type of inference: 0-confidence interval for quadratic functional; 1-signal detection; 2-FVE; 3-error variance; 4-equality of two parameters; 5-co-heritability; 6-confidence ball; 7-inference for quadratic functional}
  \item{level_significance}{nominal significance level, default value is 0.05}
  \item{null_value}{the null value for hypothesis testing}
  \item{X}{n by p matrix of predictor for model 1}
  \item{Y}{n by 1 vector of response for model 1}
  \item{V}{n by p matrix of predictor for model 2}
  \item{W}{n by 1 vector of response for model 2}
}
\details{"type_inference", "X" and "Y" are required inputs.
}
\value{
\item{eta_hat_proposed}{proposed estimator for eta}
\item{rho_hat_proposed}{proposed estimator for rho}
\item{var_proposed_test}{proposed estimator for the variance part}
\item{var_convnetional_test}{conventional estimator for the variance part}
\item{CI_proposed_method_lower}{lower bound of the confidence interval}
\item{CI_proposed_method_upper}{upper bound of the confidence interval}
\item{test_stat_proposed}{proposed test statistic}
\item{test_stat_conventional}{conventional test statistic}
\item{p_value_proposed}{p-value of the proposed test}
\item{p_value_conventional}{p-value of the conventional test}
\item{norm_beta_estimated}{estimated norm of beta}
\item{norm_gamma_estimated}{estimated norm of gamma}
\item{co_herit_hat_proposed}{proposed estimate of co-heritability}
\item{co_herit_hat_conventional}{conventional estimate of co-heritability}
\item{sigma_star_hat_square}{sigma star hat square}
\item{ID_correctly_covered_by_confidence_ball_two_sided}{indicator that the parameter is covered by the two-sided confidence ball}
\item{ID_correctly_covered_by_confidence_ball_one_sided}{indicator that the parameter is covered by the one-sided confidence ball}
}
\references{
Moderate-Dimensional Inferences on Quadratic Functionals in Ordinary Least Squares
}
\author{
Xiao Guo and Guang Cheng
}
\note{
NA
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
NA
}
\examples{
type<-4
level<-0.05
null_value<-0
n_1 <- 500
p <-100
n_2 <- 600
X<-matrix(rnorm(n_1*p),nrow=n_1,ncol=p)
eps<-matrix(rnorm(n_1),nrow=n_1,ncol=1)
Y <- X[,1] + eps
V<-matrix(rnorm(n_2*p),nrow=n_2,ncol=p)
delta<-matrix(rnorm(n_2),nrow=n_2,ncol=1)
W<-V[,1] + delta
MD_OLS(type, level, null_value, X, Y,V,W)
}
