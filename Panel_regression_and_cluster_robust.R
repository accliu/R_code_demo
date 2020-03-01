# install.packages("plm")
library(plm)
library(sandwich)
library(lmtest)
library(car)
library(numDeriv)

########  The computation of White-Arellano robust s.e.   ########

data("Produc", package = "plm")
plm_fit <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, data = Produc, effect = "individual", model = "pooling", index = c("state","year"))
# plm_fit <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, data = Produc, index = c("state","year"))
lm_fit <- lm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, data = Produc)

##equivalent calculations in plm package
plm::vcovHC(plm_fit)
plm::vcovHC(plm_fit, method = "arellano", type = "HC0", cluster = "group")
vcovG(plm_fit, cluster = "group", inner = "cluster", l = 0)

#########################################################################
##Use one way cluster from "Cluster-robust standard errors using R"
#########################################################################

##Small sample adjustment
clx <- function(fm, dfcw, cluster){
     library(sandwich)
     library(lmtest)
     M <- length(unique(cluster))
     N <- length(cluster)
     dfc <- (M/(M-1))*((N-1)/(N-fm$rank))
     u <- apply(estfun(fm),2,function(x) tapply(x, cluster, sum))
     vcovCL <- dfc*sandwich(fm, meat=crossprod(u)/N)*dfcw
     vcovCL
}

##No adjustment
clx2 <- function(fm, cluster){
  library(sandwich)
  library(lmtest)
  M <- length(unique(cluster))
  N <- length(cluster)
  u <- apply(estfun(fm),2,function(x) tapply(x, cluster, sum))
  vcovCL <- sandwich(fm, meat=crossprod(u)/N)
  vcovCL
}

M <- length(unique(Produc$state))
dfcw <- lm_fit$df/(lm_fit$df - (M -1))
clx(lm_fit, dfcw, Produc$state)
##Without small sample adjustment, the vcov equals the plm package calculation
clx2(lm_fit, Produc$state)
vcovG(plm_fit, cluster = "group", inner = "cluster", l = 0)

#########################################################################
##Use glm and glm robust formula to calculate the same thing. glm cluster robust is in Greene 14.8.4
##The code is from Molly cluster robust slides at Harvard
#########################################################################
glm_fit = glm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, family = gaussian(link = "identity"), data = Produc)

cluster_var <- Produc$state
m <- length(unique(cluster_var))
k <- length(coef(glm_fit))
u <- estfun(glm_fit)
u.clust <- matrix(NA, nrow=m, ncol=k)
for(j in 1:k){
  u.clust[,j] <- tapply(u[,j], cluster_var, sum)
}
##with df adjustment
sandwich(glm_fit,meat. = ((m/(m-1)) * t(u.clust) %*% (u.clust))/nrow(Produc))
##without. This will reproduce the vcov by plm. Note "/nrow(Produc)" is used. 
sandwich(glm_fit,meat. = crossprod(u.clust)/nrow(Produc))
sandwich(glm_fit,meat. = (t(u.clust) %*% (u.clust))/nrow(Produc))
vcovG(plm_fit, cluster = "group", inner = "cluster", l = 0)

##Equivalently, after some degree of freedom adjustment. 
summary(glm_fit) -> sum_glm
adj_fac = sum_glm$df[2]/sum(sum_glm$df[1:2])
vcov(glm_fit) %*% (t(u.clust*adj_fac) %*% (u.clust*adj_fac)) %*% vcov(glm_fit) 
##The vcov function uses already the df.residual to get the average vcov, whereas sandwich package operates on number of data points. 
 
  
##The bread of glm is not vcov extraction. A comparison:
bread(glm_fit);vcov(glm_fit)
bread(glm_fit)/vcov(glm_fit)
##The bread function in sandwich returns the total sum of the values of the score/estimation function. The vcov is already adjusted by residual degree of freedom. 


#########################################################################
##And, it's not appropriate to directly use lm fit to calculate the scores/gradients at data points 
#########################################################################
head(estfun(lm_fit))
head(estfun(glm_fit))
head(estfun(lm_fit))/head(estfun(glm_fit))
cat("They are not the same!","\nThe reason is lm uses different objective function. So, don't use these glm formula to lm object.")
##But, the difference is a scalar. 
bread(lm_fit)/vcov(lm_fit)

##Try to calculate Arellano's cluster estimator using sandwich function
n_individual = length(unique(Produc$state))
n_param = lm_fit$rank
x_u_cluster = apply(estfun(lm_fit),2,function(x,cluster_vars=cluster_var) tapply(x,list(cluster_vars),sum))
sandwich(lm_fit,meat. = t(x_u_cluster)%*%x_u_cluster/length(cluster_var))

##Try to calculate it manually. See "Robust Standard Error Estimators for Panel Models" 3.1
##The calculation of estimation/score function for linear regression. Show the equivalence. 
head(estfun(lm_fit))
head(cbind(ones=1,lm_fit$model[,-1])*lm_fit$residuals)

lm_mod_matrix = as.matrix(cbind(ones=1,lm_fit$model[,-1]))
x_u_cluster2 = apply(lm_mod_matrix*lm_fit$residuals,2,function(x,cluster_vars=cluster_var) tapply(x,list(cluster_vars),sum))
##bread is what's looked for. 
##The cluster estimator.
solve(t(lm_mod_matrix)%*%lm_mod_matrix) %*% (t(x_u_cluster2)%*%x_u_cluster2) %*% solve(t(lm_mod_matrix)%*%lm_mod_matrix)
##Yes. The match. So, the bread from the lm is the inverse of t(X)%*%X

##The glm and lm have different objective functions - one loglikehood and the other least squares. The estimated bread and meat thus will be different. The constant as pi and the sigma parameter in MLE produce a scalar difference. 
##Try to find out how glm score/estimation function is calculated. 
x <- as.matrix(cbind(1,lm_fit$model[,-1]));colnames(x)[1] = "ones"
y <- as.vector(lm_fit$model[,1])
glm_par = c(glm_fit$coefficients,sqrt(summary(glm_fit)$dispersion))

##For single observation. Make sure x, y are matrix and vector. 
# x=1;y=1;glm_par=c(1,1)
llik_gauss_single <- function(pars=glm_par,X,Y) {
  
  n_predictor = length(X)
  Y <- as.vector(Y); X <- as.vector(X)
  xbeta <- t(X)%*%as.vector(pars[1:n_predictor]); Sig <- pars[length(pars)]
  (-(1/2)*log(2*pi)-(1/2)*log(Sig^2)-(1/(2*Sig^2))*(Y-xbeta)^2)
}
# llik_gauss_single(pars=glm_par,X=x,Y=y)

gauss_grad = function(z=c(x[1,],y[1]),parm=glm_par) {
  i=z[1:(length(z)-1)]
  j=z[length(z)]
  numDeriv::grad(llik_gauss_single,method="Richardson",parm,X=i,Y=j)
}

score_glm = t(apply(cbind(x,y), 1, gauss_grad, parm=glm_par))
score_glm[1:6,-6]-head(estfun(glm_fit))
score_glm[1:6,-6]/head(estfun(glm_fit))
##The difference is due to degree of freedom adjustment.
##Truned out the function estfun adjust the score function by multipying n/(n-p). 
(score_glm[1:6,-6]/adj_fac)/head(estfun(glm_fit))


##try a different function. The same result. 
# library(pracma)
# gauss_grad = function(z=c(x[1,],y[1]),parm=glm_par) {
#   i=z[1:(length(z)-1)]
#   j=z[length(z)]
#   pracma::grad(llik_gauss_single,parm,X=i,Y=j)
# }

score_glm = t(apply(cbind(x,y), 1, gauss_grad, parm=glm_par))
score_glm[1:6,-6]/head(estfun(glm_fit))

###############################################################
## Use optim to estimate hessian
###############################################################
# First, put the data into matrices for the MLE procedure
x <- cbind(1,lm_fit$model[,-1]);colnames(x)[1] = "ones"
y <- lm_fit$model[,1]

llik_gauss <- function(pars,X=x,Y=y) {
  
  n_predictor = ncol(X); n_par = n_predictor + 1
  
  Y <- as.vector(Y); X <- as.matrix(X)
  
  xbeta <- X%*%pars[1:n_predictor]; Sig <- pars[n_par]
  
  sum(-(1/2)*log(2*pi)-(1/2)*log(Sig^2)-(1/(2*Sig^2))*(Y-xbeta)^2)
}

lm_mle <- optim(rep(1,ncol(x)+1),llik_gauss, method = "BFGS", control = list(trace=1,maxit=100,fnscale = -1), hessian = TRUE)
##If default values of X,Y are not supplied in objective function
lm_mle <- optim(rep(1,ncol(x)+1),llik_gauss,X=x,Y=y, method = "BFGS", control = list(trace=1,maxit=100,fnscale = -1), hessian = TRUE)

vcov_mle <- -solve(lm_mle$hessian)
se_mle <- sqrt( diag(vcov_mle))
vcov_mle[-6,-6]->vcovMLE
vcovMLE
vcov(glm_fit)
vcov(lm_fit)

##The difference may be due to degree of freedom adjustment? 
vcov_mle[-6,-6]->vcovMLE
vcovMLE

##Try optim procedure parameters in robust error calculation. Don't divide by adj factor if use robust formula below. 
score_glm = t(apply(cbind(x,y), 1, gauss_grad, parm=lm_mle$par))

cluster_var <- Produc$state
m <- length(unique(cluster_var))
k <- length(coef(glm_fit))
u <- score_glm
u.clust <- matrix(NA, nrow=m, ncol=k)
for(j in 1:k){
  u.clust[,j] <- tapply(u[,j], cluster_var, sum)
}

##A comparison
vcovG(plm_fit, cluster = "group", inner = "cluster", l = 0) -> cluster_plm

##Equivalently, after some degree of freedom adjustment.
vcovMLE %*% (t(u.clust) %*% (u.clust)) %*% vcovMLE -> cluster_mle

cluster_plm - cluster_mle
cluster_plm/cluster_mle

#########################################################################
##They are fairly close. Replicated glm cluster robust. 
#########################################################################



#########################################################################
##Try a few test.  
#########################################################################
plm_fit_pool <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, data = Produc, effect = "individual", model = "pooling", index = c("state","year"))
plm_fit_fixed <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, data = Produc, effect = "individual", model = "within", index = c("state","year"))
plm_fit_random <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, data = Produc, effect = "individual", model = "random", index = c("state","year"))


##F test for individual effect restriction. Restrict individual effect equal in the fixed effects model. 
pFtest(plm_fit_fixed,plm_fit_pool)
##test if random effect is present. Supply a formula or any model above. It's the formula that matters. 
plmtest(plm_fit_random,type = "bp")
##Hausman for fix and random
phtest(plm_fit_fixed,plm_fit_random)
##Test serial correlation
pbgtest(plm_fit_pool)
pbgtest(plm_fit_fixed)
pbgtest(plm_fit_random)







