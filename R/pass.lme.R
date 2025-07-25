##################################################
#' @title Power and Sample Size for Linear Mixed Effect Models Package  \cr
#' @author Marco Chak Yan YU \cr
#' Maintainer: Marco Chak Yan YU <marcocyyu@gmail.com>
#' @note License: GPL-3
#' @seealso \code{\link{lme.Lb.dist.theta}}, \code{\link{pass.lme.CLb.test}}
#' @details \href{https://drive.google.com/open?id=1MQvzp9na-Zm9mFT5jC1I_3U0JHEFky_k}{Technical note}
#' @description
#' Consider the following model: \cr \cr
#'       Y=XB+Zu+e, u~N(0,D), e~N(0,R) \cr
#'       Yi~N(XBi,Vi), Vi=Zi*D*Zi'+Ri, \cr \cr
#'           for each independent observation i \cr \cr \cr
#' estimate of fixed effect coefficient B, denoted by b: \cr \cr
#'       b=inv(sum(Xi'*inv(Vi)*Xi))*(sum(Xi'*inv(Vi)*Yi)) \cr \cr \cr
#' variance of b: \cr \cr
#'       var(b)=Vb/n=inv(sum(Xi'*inv(Vi)*Xi)) \cr \cr
#'           where Vb=inv(Xi'*inv(Vi)*Xi) \cr \cr \cr
#' distribution of any linear combinations L of b is given by: \cr \cr
#'       Lb~N(mu,Sigma/n) \cr \cr
#'             where mu = LB, Sigma = L*Vb*L' \cr \cr \cr
#' Interested parameters/linear combinations LB from more than one
#' independent populations can be aggregrate togeter by
#' appending mu vertically and Sigma/n diagonally \cr \cr
#' Consider Lb~N(MU,SIGMA) as the aggregrated estimates\cr
#' Any comparison of interested parameters can be formulated by
#' multiplying a contrast matrix C on LB and set \cr \cr
#' H0: C*LB=d for any vector of value d to be tested \cr \cr
#' We then have \cr \cr
#'       C*Lb~N(C*MU,C*SIGMA*C') \cr \cr
#' and \cr \cr
#'       (C*Lb-d)'*inv(C*SIGMA*C')*(C*Lb-d)~chisq(q,lambda) \cr \cr
#' where degree of freedom q=rank(C*SIGMA*C'), \cr
#' non-centrality parameter lambda=(C*LB-d)'*inv(C*SIGMA*C')*(C*LB-d) \cr \cr \cr
#' Power of the test H0 is given by 1-beta=P(chisq(q,lambda)>qchisq(1-alpha,lambda)) \cr
#' Required sample size for desired power can be obtained by bisection method.
##
.onAttach <- function(libname,pkgname) {
  packageStartupMessage('Title: Power and Sample Size for Linear Mixed Effect Models')
  packageStartupMessage('Author: Marco Chak Yan YU')
  packageStartupMessage('Maintainer: Marco Chak Yan YU <marcocyyu@gmail.com>')
  packageStartupMessage('License: GPL-3')
  packageStartupMessage('see also: lme.Lb.dist.theta, pass.lme.CLb.test')
}


##################################################
#' @title Calculate mean and variance  \cr
#' for linear combination of the  \cr
#' Best Linear Unbiased Estimator (BLUE) \cr
#' for Linear Mixed Effect (LME) Model
#' @description
#' Consider the following model: \cr \cr
#'       Y=XB+Zu+e, u~N(0,D), e~N(0,R) \cr
#'       Yi~N(XBi,Vi), Vi=Zi*D*Zi'+Ri, \cr \cr
#'           for each independent observation i \cr \cr \cr
#' estimate of fixed effect coefficient B, denoted by b: \cr \cr
#'       b=inv(sum(Xi'*inv(Vi)*Xi))*(sum(Xi'*inv(Vi)*Yi)) \cr \cr \cr
#' variance of b: \cr \cr
#'       var(b)=Vb/n=inv(sum(Xi'*inv(Vi)*Xi)) \cr \cr
#'           where Vb=inv(Xi'*inv(Vi)*Xi) \cr \cr \cr
#' distribution of any linear combinations L of b is given by: \cr \cr
#'       Lb~N(mu,Sigma/n) \cr \cr
#'             where mu = LB, Sigma = L*Vb*L' \cr \cr \cr
#' @param B      fixed effect beta in px1 matrix
#' @param D      list of qxq random effect variance matrix;
#'                  where the first element corresponding to the highest-level effect,
#'                  the last element corresponding to the level 1 effect \cr \cr
#' @param R      residual variance \cr \cr
#' @param L      lxp matrix, representing l-linear-combinations of beta interested, \cr
#'                  if L is not defined, it will be auto-created to select the last coefficient \cr \cr
#' @param X      nxp matrix representing the covariates for the fixed effects \cr \cr
#' @param Z      nxq matrix representing the covariates for each level of random effects \cr \cr
#' @param m      vector of repeated measures from the highest to lowest level
#'                  (level 1 effects are addressed by Z and X and no need to be specified) \cr \cr
#' @return
#' theta:  parameters (mu and Sigma) of the normal distribution for Lb
#' @seealso \code{\link{pass.lme.CLb.test}}
#' @details \href{https://drive.google.com/open?id=1MQvzp9na-Zm9mFT5jC1I_3U0JHEFky_k}{Technical note}
#' @examples
#' #Example 1
#' # calc BLUE for 1-level LME model,
#' # with covariates X, Z: (1,t), t=1,2,3
#' # for both fixed and random effects,
#' # with fixed effect coefficients B: (100,-0.5),
#' # random effect variance D: (2 1;1 2),
#' # residual variance R: 0.2
#' B <- matrix(c(100,-0.5),2,1)
#' D <- matrix(c(2,1,1,2),2,2)
#' R <- 0.2
#' X <- cbind(rep(1,3),1:3)
#' Z <- X
#' lme.Lb.dist.theta(B,D,R,X,Z)
#'
#' #Example 2
#' # calc BLUE for 3-levels LME model,
#' # with level 1 same as the above example
#' # with 3 repeated-measures in level 2
#' # and 5 repeated-measures in the highest level,
#' # assuming random effect variance for level 2 = (3 1;1 3),
#' # and random effect variance for highest level = (5 1;1 5)
#' D <- list(matrix(c(2,1,1,2),2,2),matrix(c(3,1,1,3),2,2), matrix(c(5,1,1,5),2,2))
#' lme.Lb.dist.theta(B,D,R,X,Z,m=c(5,3))
##
lme.Lb.dist.theta <- function(B,D,R,X,Z,m=NULL,L) {
  ##########
  # format input
  B <- matrix(c(B),ncol=1)
  if (missing(L)) {
    L <- matrix(c(rep(0,nrow(B)-1),1),nrow=1)
  } else if (is.vector(L)) {
    L <- matrix(c(L),nrow=1)
  } else if (dim(L)[2]!=dim(B)[1]) {
    L <- t(L)
  }
  if (!is.list(D)) D <- list(as.matrix(D))
  ##########
  # input check
  if (dim(B)[1]!=dim(X)[2]) {
    stop('inconsistent matrix dimension between X and B, requires dim(X)[2] == dim(B)[1]')
  } else if (dim(L)[2]!=dim(B)[1]) {
    stop('inconsistent matrix dimension between L and B, requires dim(L)[2] == dim(B)[1]')
  } else if (length(m)!=length(D)-1) {
    stop('inconsistent input length between D and m, requires length(D)-1 == length(m)')
  } else if (any(sapply(
    1:length(D),function(i)
      (dim(D[[i]])[1]!=dim(D[[i]])[2])&
    (dim(D[[i]])[1]!=dim(Z)[2])))) {
    stop('inconsistent matrix dimension between D and Z, requires dim(Z)[2] == dim(D[[i]])[1] == dim(D[[i]])[2] for all i')
  } else if (dim(X)[1]!=dim(Z)[1]) {
    stop('inconsistent matrix dimension between X and Z, requires dim(X)[1] == dim(Z)[1]')
  } else if (length(c(R))!=1) {
    stop('wrong input R, requires a single number')
  } else {
    ##########
    # perform the calculation
    LB <- L%*%B
    # calc Zi*D*Zi' in Vi=Zi*D*Zi'+Ri
    Zi <- Z
    Xi <- X
    Vi <- Zi %*% matrix(D[[length(D)]],dim(D[[length(D)]])) %*% t(Zi)
    if (length(m)>=1) {
      for (i in seq(length(m),1,-1)) {
        Zi <- kronecker(matrix(rep(1,m[i]),ncol=1),Zi)
        Xi <- kronecker(matrix(rep(1,m[i]),ncol=1),Xi)
        Vi <- kronecker(diag(rep(1,m[i])),Vi)
        Vi <- Vi + Zi %*% matrix(D[[i]],dim(D[[i]])) %*% t(Zi)
      }
    }
    # calc Vi=Zi*D*Zi'+Ri
    Vi <- Vi + diag(R,dim(Vi)[1])
    # Vb/n=var(b)=inv(sum(Xi'*inv(Vi)*Xi))
    Vb <- solve(t(Xi) %*% solve(Vi) %*% Xi)
    # Lb~N(LB,L*Vb*L'/n)
    VLb <- L %*% Vb %*% t(L)
    # generate output
    theta <- NULL
    theta$mu <- LB
    theta$Sigma <- VLb
    return(theta)
  }
}


##################################################
#' @title Calculate Power or Sample Size required \cr
#' for Contrasts of linear combinations \cr
#' of fixed effect parameters \cr
#' in Linear Mixed Effect (LME) Model
#' @description
#' Interested parameters/linear combinations LB from more than one
#' independent populations can be aggregrate togeter by
#' appending mu vertically and Sigma/n diagonally \cr \cr
#' Consider Lb~N(MU,SIGMA) as the aggregrated estimates\cr
#' Any comparison of interested parameters can be formulated by
#' multiplying a contrast matrix C on LB and set \cr \cr
#' H0: C*LB=d for any vector of value d to be tested \cr \cr
#' We then have \cr \cr
#'       C*Lb~N(C*MU,C*SIGMA*C') \cr \cr
#' and \cr \cr
#'       (C*Lb-d)'*inv(C*SIGMA*C')*(C*Lb-d)~chisq(q,lambda) \cr \cr
#' where degree of freedom q=rank(C*SIGMA*C'), \cr
#' non-centrality parameter lambda=(C*LB-d)'*inv(C*SIGMA*C')*(C*LB-d) \cr \cr \cr
#' Power of the test H0 is given by 1-beta=P(chisq(q,lambda)>qchisq(1-alpha,lambda)) \cr
#' Required sample size for desired power can be obtained by bisection method.
#' @param thetas  list of theta (LB and VLb), can be different for each group \cr \cr
#' @param C       Contrast of Matrix \cr \cr
#' @param d       Value vector to be tested for all contrast \cr \cr
#' @param alpha   significant level \cr \cr
#' @param power   desired power for sample size calculation \cr \cr
#' @param n       sample size for power calculation / \cr
#'                    or sample size ratio with power for sample size calculation
#'                       (NULL for balanced design) \cr \cr
#' @return
#'   solved.power  given sample size n, this gives the power for testing H0 \cr
#'   solved.n      given the desired power, this gives the sample size for H0 \cr
#' @seealso \code{\link{lme.Lb.dist.theta}}
#' @details \href{https://drive.google.com/open?id=1MQvzp9na-Zm9mFT5jC1I_3U0JHEFky_k}{Technical note}
#' @examples
#' #Example 1 (test fixed effect coefficient 2=0) with power of 80%
#' # for 1-level LME model, with covariates X, Z: (1,t), t=1,2,3
#' # for both fixed and random effects, with fixed effect coefficients B: (100,-0.5),
#' # random effect variance D: (2 1;1 2), residual variance R: 0.2
#' B <- matrix(c(100,-0.5),2,1)
#' D <- matrix(c(2,1,1,2),2,2)
#' R <- 0.2
#' X <- cbind(rep(1,3),1:3)
#' Z <- X
#' theta <- lme.Lb.dist.theta(B,D,R,X,Z)
#' pass.lme.CLb.test(list(theta),alpha=0.05,power=0.8)
#' pass.lme.CLb.test(list(theta),alpha=0.05,n=66)
#'
#' #Example 2 (compare two fixed effect coefficient 2) with power of 80%
#' # Consider above model as a control group model,
#' # with an independent treatment group with model same as the control
#' # except a different fixed effect coefficient 2 for treatment
#' # = fixed effect coefficient 2 for control x 0.7
#' theta2 <- theta
#' theta2$mu <- theta$mu *0.7
#' C <- matrix(c(1,-1),1,2)
#' pass.lme.CLb.test(list(theta,theta2),C,alpha=0.05,power=0.8)
#' pass.lme.CLb.test(list(theta,theta2),C,alpha=0.05,n=1468)
#'
#' #Example 3 (compare two fixed effect coefficient 2) with power of 80%
#' # with sample size ratio, control:treatment = 1:2
#' pass.lme.CLb.test(list(theta,theta2),C,alpha=0.05,power=0.8,n=c(1,2))
#' pass.lme.CLb.test(list(theta,theta2),C,alpha=0.05,n=c(1101,2202))
#'
#' #Example 4 (repeated-measures ANOVA for comparing 3 group means) with power of 80%
#' # for 1-level LME model with mean for group 1, 2 and 3 are 100, 99, 102, respectively,
#' # each subject to be measured 2 times, with within-subject variance = 15, residual variance = 10
#' B <- 100
#' D <- 15
#' R <- 10
#' X <- matrix(1,2,1)
#' Z <- X
#' theta <- lme.Lb.dist.theta(B,D,R,X,Z)
#' theta2 <- theta
#' theta3 <- theta
#' theta2$mu <- 99
#' theta3$mu <- 102
#' C <- rbind(c(1,-1,0),c(1,0,-1))
#' pass.lme.CLb.test(list(theta,theta2,theta3),C,alpha=0.05,power=0.8)
#' pass.lme.CLb.test(list(theta,theta2,theta3),C,alpha=0.05,n=41)
##
pass.lme.CLb.test <- function(thetas,C=NULL,d=NULL,alpha=0.05,power=NULL,n=NULL) {
  if (is.null(n)) {
    n <- rep(1,length(thetas))
  } else if ((length(n)==1)&(length(thetas)>1)) {
    n <- rep(n,length(thetas))
  }
  LB <- NULL
  for (i in 1:length(thetas)) {
    LB <- rbind(LB,thetas[[i]]$mu)
  }
  VLb <- matrix(0,dim(LB)[1],dim(LB)[1])
  j <- 0
  for (i in 1:length(thetas)) {
    VLbi <- thetas[[i]]$Sigma
    VLb[(j+1):(j+dim(VLbi)[1]),(j+1):(j+dim(VLbi)[1])] <- VLbi/n[i]
    j <- j+dim(VLbi)[1]
  }

  if (is.null(C)) C <- diag(dim(LB)[1])
  if (is.null(d)) d <- matrix(0,nrow(C),1)
  VCLb <- C%*%VLb%*%t(C)
  df <- qr(VCLb)$rank
  chisq.ncp <- t(C%*%LB-d) %*% solve(VCLb) %*% (C%*%LB-d)

  if (is.null(power)) {
    solved.power <- 1-pchisq(qchisq(1-alpha,df,0),df,chisq.ncp)
    return(solved.power)
  } else {
    solved.n <- 100
    solved.n.lb <- 0
    solved <- FALSE
    while (!solved) {
      if ((1-pchisq(qchisq(1-alpha,df,0),df,chisq.ncp*solved.n) < power)) {
        solved.n.lb <- solved.n
        solved.n <- solved.n*10
      } else if ((1-pchisq(qchisq(1-alpha,df,0),df,chisq.ncp*solved.n) > power+0.001)) {
        solved.n <- (solved.n + solved.n.lb)/2
      } else {
        solved <- TRUE
      }
    }
    solved.n <- solved.n *n
    return(solved.n)
  }
}

