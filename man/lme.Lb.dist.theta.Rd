% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pass.lme.R
\name{lme.Lb.dist.theta}
\alias{lme.Lb.dist.theta}
\title{Calculate mean and variance  \cr
for linear combination of the  \cr
Best Linear Unbiased Estimator (BLUE) \cr
for Linear Mixed Effect (LME) Model}
\usage{
lme.Lb.dist.theta(B, D, R, X, Z, m = NULL, L)
}
\arguments{
\item{B}{fixed effect beta in px1 matrix}

\item{D}{list of qxq random effect variance matrix;
where the first element corresponding to the highest-level effect,
the last element corresponding to the level 1 effect \cr \cr}

\item{R}{residual variance \cr \cr}

\item{X}{nxp matrix representing the covariates for the fixed effects \cr \cr}

\item{Z}{nxq matrix representing the covariates for each level of random effects \cr \cr}

\item{m}{vector of repeated measures from the highest to lowest level
(level 1 effects are addressed by Z and X and no need to be specified) \cr \cr}

\item{L}{lxp matrix, representing l-linear-combinations of beta interested, \cr
if L is not defined, it will be auto-created to select the last coefficient \cr \cr}
}
\value{
theta:  parameters (mu and Sigma) of the normal distribution for Lb
}
\description{
Consider the following model: \cr \cr
      Y=XB+Zu+e, u~N(0,D), e~N(0,R) \cr
      Yi~N(XBi,Vi), Vi=Zi*D*Zi'+Ri, \cr \cr
          for each independent observation i \cr \cr \cr
estimate of fixed effect coefficient B, denoted by b: \cr \cr
      b=inv(sum(Xi'*inv(Vi)*Xi))*(sum(Xi'*inv(Vi)*Yi)) \cr \cr \cr
variance of b: \cr \cr
      var(b)=Vb/n=inv(sum(Xi'*inv(Vi)*Xi)) \cr \cr
          where Vb=inv(Xi'*inv(Vi)*Xi) \cr \cr \cr
distribution of any linear combinations L of b is given by: \cr \cr
      Lb~N(mu,Sigma/n) \cr \cr
            where mu = LB, Sigma = L*Vb*L' \cr \cr \cr
}
\details{
\href{https://drive.google.com/open?id=1MQvzp9na-Zm9mFT5jC1I_3U0JHEFky_k}{Technical note}
}
\examples{
#Example 1
# calc BLUE for 1-level LME model,
# with covariates X, Z: (1,t), t=1,2,3
# for both fixed and random effects,
# with fixed effect coefficients B: (100,-0.5),
# random effect variance D: (2 1;1 2),
# residual variance R: 0.2
B <- matrix(c(100,-0.5),2,1)
D <- matrix(c(2,1,1,2),2,2)
R <- 0.2
X <- cbind(rep(1,3),1:3)
Z <- X
lme.Lb.dist.theta(B,D,R,X,Z)

#Example 2
# calc BLUE for 3-levels LME model,
# with level 1 same as the above example
# with 3 repeated-measures in level 2
# and 5 repeated-measures in the highest level,
# assuming random effect variance for level 2 = (3 1;1 3),
# and random effect variance for highest level = (5 1;1 5)
D <- list(matrix(c(2,1,1,2),2,2),matrix(c(3,1,1,3),2,2), matrix(c(5,1,1,5),2,2))
lme.Lb.dist.theta(B,D,R,X,Z,m=c(5,3))
}
\seealso{
\code{\link{pass.lme.CLb.test}}
}
