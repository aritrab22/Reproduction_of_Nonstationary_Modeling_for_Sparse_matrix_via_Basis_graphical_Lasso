#######################################################################
BGLBasisSetup <- function(y,locs,basis="LatticeKrig",Phi=NULL,crossvalidation=FALSE,distance.penalty=FALSE,...)
{
  tmp <- match.call(expand.dots=TRUE)
  if(!is.null(Phi)){
    returnPhi = FALSE
    Phi_Phi <- crossprod(Phi)
    Phi_y <- crossprod(Phi, y)
    trS <- norm(y, "F")^2/ifelse(ncol(y)==1,1,ncol(y)-1)
    
  } else if(is.null(Phi) & basis=="LatticeKrig"){
    returnPhi=TRUE
    if(!all(c("NC","nlevel") %in% names(tmp))){
      stop("Need to specify NC and nlevel for LatticeKrig basis")
    }
    sDomain <- apply(locs, 2, "range")
    # a.wght and nu have no effect on the basis structure but must be specified for LKrigSetup to run
    
    LKinfo <- LKrigSetup(sDomain, a.wght=4.05, nu=0.5,...)
    Phi <- as.matrix(LKrig.basis(locs, LKinfo))
    Phi_Phi <- crossprod(Phi)
    Phi_y <- crossprod(Phi, y)
    trS <- norm(y, "F")^2/ifelse(ncol(y)==1,1,ncol(y)-1)
  } else {
    stop("LatticeKrig is the only currently supported basis")
  }
  
  if(crossvalidation){
    mylist <- list(Phi_Phi=Phi_Phi, Phi_y=Phi_y, trS=trS)
  } else {
    Phi_S_Phi <- tcrossprod(Phi_y)/ ifelse(ncol(Phi_y)==1,1,ncol(Phi_y)-1)
    mylist <- list(Phi_Phi=Phi_Phi, Phi_S_Phi=Phi_S_Phi, trS=trS)
  }
  
  if(distance.penalty){
    lookvec <- vector('list',LKinfo$nlevel)
    for(i in 1:LKinfo$nlevel){
      look <- LKrigLatticeCenters(LKinfo,Level=i)
      lookvec[[i]] <- expand.grid(look)
    }
    
    basiscenters <- as.matrix(do.call("rbind",lookvec))
    rm(look,lookvec)
    mylist$basiscenters <- basiscenters
    mylist$basisdistancematrix <- rdist(basiscenters)
  }
  
  if(returnPhi) {
    mylist$Phi <- Phi
  }
  
  return(mylist)
}

########################################################################

basis.setup <- BGLBasisSetup(y=dat, locs=locs, nlevel=1, NC=2)
rowSums(basis.setup$Phi)
names(basis.setup)

########################################################################
nugget_estimate <- function(Phi_Phi,Phi_S_Phi,trS,n)
{
  l <- dim(Phi_Phi)[1]
  
  nugget_function <- function(par)
  {
    diag_parameter <- par[1]
    nug <- par[2]
    Phi_Phi_over_nugget <- Phi_Phi / nug
    Phi_S_Phi_over_nugsq <- Phi_S_Phi / (nug^2)
    cholQ_plus_Phi <- chol(diag_parameter * diag(l) + Phi_Phi_over_nugget)
    
    return(2*sum(log(diag(cholQ_plus_Phi))) - l*log(diag_parameter) - sum(diag(backsolve(cholQ_plus_Phi,backsolve(cholQ_plus_Phi,Phi_S_Phi_over_nugsq,transpose=TRUE),transpose=FALSE))) + n*log(nug) + trS/nug)
  }
  
  out <- optim(par=c(1,0.1),fn=nugget_function,hessian=FALSE,method="L-BFGS-B",lower=c(0.01,0.01))
  my_tau_sq <- out$par[2]
  return(my_tau_sq)
}

#######################################################################

basis.setup <- BGLBasisSetup(y=tmin$data,locs=tmin$lon.lat.proj,basis="LatticeKrig",
                                   crossvalidation=FALSE, NC=30,nlevel=1)
Phi_Phi <- basis.setup$Phi_Phi
Phi_S_Phi <- basis.setup$Phi_S_Phi
trS <- basis.setup$trS
tau_sq <- nugget_estimate(Phi_Phi,Phi_S_Phi,trS,n=dim(tmin$lon.lat)[1])

#######################################################################

penaltymatrixsetup <- function(lambda,zero.diagonal.penalty=TRUE,l,basisdistancematrix=NULL)
{
  # Set up penalty matrix if not prespecified and do not penalize marginal precision parameters
  
  if(length(lambda)==1 & !is.null(basisdistancematrix))
  {
    biglambda <- lambda*basisdistancematrix
  }  else if(length(lambda)==1 & is.null(basisdistancematrix))
  {
    biglambda <- matrix(lambda,nrow=l,ncol=l)
  } else {
    biglambda <- lambda
  }
  
  if(zero.diagonal.penalty==TRUE)
  {
    diag(biglambda) <- 0
  }
  
  return(biglambda)
}

#########################################################

BGL_DC <- function(lambda,Phi_Dinv_Phi,Phi_Dinv_S_Dinv_Phi,guess,outer_tol=NULL,MAX_ITER=NULL,MAX_RUNTIME_SECONDS=NULL,verbose=TRUE)
{
  
  if(is.null(outer_tol))
  {
    outer_tol=0.01
  }
  if(is.null(MAX_ITER))
  {
    MAX_ITER=50
  }
  if(is.null(MAX_RUNTIME_SECONDS))
  {
    MAX_RUNTIME_SECONDS=604800
  }
  
  my_norm = 1
  counter = 1L
  start.time <- Sys.time()
  elapsed <- 0
  while (my_norm >= outer_tol && counter <= MAX_ITER && elapsed <= MAX_RUNTIME_SECONDS)
  {
    
    M = chol2inv(chol(guess+Phi_Dinv_Phi))
    S_star <-  M + M %*% Phi_Dinv_S_Dinv_Phi %*% M
    
    new_guess <- QUIC(S_star, rho= lambda,msg=0)$X
    current.time <- Sys.time()
    elapsed <- difftime(current.time,start.time,units="secs")
    my_norm <- norm(new_guess - guess,type="F")/norm(guess,type="F")
    
    if(verbose) {
      cat('Iteration ', counter, '. Relative error: ', my_norm, '. Time elapsed: ', elapsed,"\n", sep="")
    }
    
    guess <- new_guess
    counter <- counter+1
  }
  
  if(my_norm < outer_tol){
    cat("Convergence. Relative error: ", sprintf("%10f",my_norm), sep = '')
    cat("\n")
  }  else if(elapsed>MAX_RUNTIME_SECONDS) {
    cat("MAX RUNTIME REACHED. Relative error: ", sprintf("%10f",my_norm), sep = '')
    cat("\n")
  } else {
    cat("MAX ITERATIONS REACHED. Relative error: ", sprintf("%10f",my_norm), sep = '')
    cat("\n")
  }
  return(guess)
}

###################################################

basis.setup <- BGLBasisSetup(y=tmin$data,locs=tmin$lon.lat.proj,basis="LatticeKrig", 
                                   crossvalidation=FALSE,NC=20,nlevel=1)
Phi_Phi <- basis.setup$Phi_Phi
Phi_S_Phi <- basis.setup$Phi_S_Phi
tau_sq <- 2
lambda <- matrix(10,nrow=dim(Phi_Phi)[1],ncol=dim(Phi_Phi)[1])
diag(lambda) <- 0
BGLguess <- BGL_DC(lambda=lambda,Phi_Dinv_Phi=Phi_Phi/tau_sq,
                   Phi_Dinv_S_Dinv_Phi=Phi_S_Phi/(tau_sq^2), guess=diag(dim(Phi_Phi)[1]),
                  outer_tol=0.05,MAX_ITER=50,MAX_RUNTIME_SECONDS=86400)

#############################################################

unpenalized_neglikelihood_Q <- function(Q,Phi_Dinv_Phi,Phi_Dinv_S_Dinv_Phi)
{
  cholQ_plus_Phi <- chol(Q + Phi_Dinv_Phi)
  return(2*sum(log(diag(cholQ_plus_Phi))) - 2*sum(log(diag(chol(Q)))) - sum(diag(backsolve(cholQ_plus_Phi,backsolve(cholQ_plus_Phi,Phi_Dinv_S_Dinv_Phi,transpose=TRUE),transpose=FALSE))))
}

BGL_CV <- function(kfolds=5,y,locs,lambdalist,zero.diagonal.penalty=TRUE,basis="LatticeKrig",Phi=NULL,guess=NULL,outer_tol=NULL,MAX_ITER=NULL,MAX_RUNTIME_SECONDS=NULL,tau_sq=NULL,verbose=TRUE,distance.penalty=FALSE,final.guess=TRUE,...)
{
  basis.setup <- BGLBasisSetup(y=y,locs=locs,basis=basis,Phi=Phi,crossvalidation=TRUE,distance.penalty=distance.penalty,...)
  Phi_Phi <- basis.setup$Phi_Phi
  Phi_y <- basis.setup$Phi_y
  Phi_S_Phi <- tcrossprod(Phi_y)/ifelse(ncol(Phi_y)==1,1,ncol(Phi_y)-1)
  trS <- basis.setup$trS
  basiscenters <- basis.setup$basiscenters
  basisdistancematrix <- basis.setup$basisdistancematrix
  
  if(is.null(tau_sq))
  {
    cat('Estimating nugget effect: ')
    tau_sq <- nugget_estimate(Phi_Phi,Phi_S_Phi,trS,n=dim(locs)[1])
    cat(tau_sq, '\n', sep="")
  }
  
  if(is.null(guess))
  {
    guess <- diag(dim(Phi_Phi)[1])
  }
  
  set.seed(1)
  folds <- sample(rep(1:kfolds,length.out=ncol(y)))
  likelihood_matrix <- matrix(0,nrow=kfolds,ncol=length(lambdalist))
  
  if(zero.diagonal.penalty==TRUE)
  {  
    zero.diagonal.penalty.list <- vector("list",length=length(lambdalist))
    for(i in 1:length(lambdalist))
    {
      thislambda <- lambdalist[[i]]
      if(length(c(thislambda))==1)
      {
        thislambda <- matrix(thislambda,nrow=dim(Phi_Phi)[1],ncol=dim(Phi_Phi)[1])
        diag(thislambda) <- 0
      } else {
        diag(thislambda) <- 0
      }
      zero.diagonal.penalty.list[[i]] <- thislambda
    }
  }
  
  cat("Starting CV...")
  
  for(k in 1:kfolds){
    test.indices <- which(folds==k)
    train.indices <- which(folds!=k)
    
    Phi_Stest_Phi <- tcrossprod(Phi_y[, test.indices])/ ifelse(ncol(Phi_y[,test.indices])==1,1,ncol(Phi_y[,test.indices])-1)
    Phi_Strain_Phi <- tcrossprod(Phi_y[, train.indices])/ ifelse(ncol(Phi_y[,train.indices])==1,1,ncol(Phi_y[,train.indices])-1)
    
    for(penaltyindex in 1:length(lambdalist))
    {
      penaltymatrix <- penaltymatrixsetup(lambdalist[[penaltyindex]],zero.diagonal.penalty,l,basisdistancematrix)
      foldguess <- BGL_DC(penaltymatrix,Phi_Dinv_Phi=Phi_Phi/tau_sq,Phi_Dinv_S_Dinv_Phi=Phi_Strain_Phi/(tau_sq^2),guess=guess,outer_tol,MAX_ITER,MAX_RUNTIME_SECONDS)
      likelihood_matrix[k,penaltyindex] <- unpenalized_neglikelihood_Q(foldguess,Phi_Dinv_Phi=Phi_Phi/tau_sq,Phi_Dinv_S_Dinv_Phi=Phi_Stest_Phi/(tau_sq^2))
    }
  }
  
  cat("Starting final estimate with best penalty parameter:\n")
  
  collapse_likelihoods <- colSums(likelihood_matrix)/k
  best_index <- which.min(collapse_likelihoods)
  best_lambda <- penaltymatrixsetup(lambdalist[[best_index]],zero.diagonal.penalty,l,basisdistancematrix)
  
  mylist <- list(nugget_variance=tau_sq,Phi=basis.setup$Phi,best_lambda=best_lambda,best_lambda_index=best_index,likelihood_matrix=likelihood_matrix)
  
  if(final.guess)
  {
    finalguess <- BGL_DC(best_lambda,Phi_Dinv_Phi=Phi_Phi/tau_sq,Phi_Dinv_S_Dinv_Phi=Phi_S_Phi/(tau_sq^2),guess=guess,outer_tol,MAX_ITER,MAX_RUNTIME_SECONDS)
    mylist$Q=finalguess
  }
  
  return(mylist)
  
}

###############################################################################

#the full search space considered, commented for runtime, will take very long
#CVguess <- BGL_CV(kfolds=2, y=tmin$data, locs=tmin$lon.lat.proj, lambdalist=1:30, basis="LatticeKrig",
#outer_tol=0.05, MAX_ITER=50, MAX_RUNTIME_SECONDS=86400, NC=30, nlevel=1,distance.penalty=TRUE)

############################################################################

BGL <- function(y,locs,lambda,zero.diagonal.penalty=TRUE,basis="LatticeKrig",Phi=NULL,guess=NULL,outer_tol=NULL,MAX_ITER=NULL,MAX_RUNTIME_SECONDS=NULL,tau_sq=NULL,verbose=TRUE,distance.penalty=FALSE,...)
{
  # Set up basis
  basis.setup <- BGLBasisSetup(y=dat,locs=locs,basis=basis,Phi=Phi,crossvalidation=FALSE,distance.penalty=distance.penalty,...)
  Phi_Phi <- basis.setup$Phi_Phi
  Phi_S_Phi <- basis.setup$Phi_S_Phi
  trS <- basis.setup$trS
  basiscenters <- basis.setup$basiscenters
  basisdistancematrix <- basis.setup$basisdistancematrix
  
  # Estimate nugget effect variance
  if(is.null(tau_sq))
  {
    cat('Estimating nugget effect variance: ')
    tau_sq <- nugget_estimate(Phi_Phi,Phi_S_Phi,trS,n=dim(locs)[1])
    cat(tau_sq, '\n', sep="")
  }
  
  # Set up identity starting matrix for QUIC algorithm if not prespecified
  l <- dim(Phi_Phi)[1]
  if(is.null(guess))
  {
    guess <- diag(l)
  }
  
  penaltymatrix <- penaltymatrixsetup(lambda=lambda,zero.diagonal.penalty=TRUE,l=l,basisdistancematrix=basisdistancematrix)
  
  # Execute DC algorithm
  cat('Running basis graphical lasso algorithm:\n')
  finalguess <- BGL_DC(penaltymatrix,Phi_Dinv_Phi=Phi_Phi/tau_sq,Phi_Dinv_S_Dinv_Phi=Phi_S_Phi/(tau_sq^2),guess=guess,outer_tol,MAX_ITER,MAX_RUNTIME_SECONDS)
  
  # Return estimated precision matrix and nugget effect variance
  mylist <- list(Q=finalguess,nugget_variance=tau_sq,Phi=basis.setup$Phi)
  
  return(mylist)
}

precision.fit <- BGL(y=tmin$data, locs=tmin$lon.lat.proj, lambda=7, basis="LatticeKrig",
                         distance.penalty=TRUE,outer_tol=5e-2, MAX_ITER=50,
                         MAX_RUNTIME_SECONDS=86400, NC=20, nlevel=1)
                         
# The estimate in paper, with the tolerance set to 1e-2 not 5e-2 and NC=30 not NC=20.
# Takes a few minutes longer than the version above
# precision.fit <- BGL(y=tmin$data, locs=tmin$lon.lat.proj, lambda=7, basis="LatticeKrig",
# distance.penalty=TRUE,outer_tol=1e-2, MAX_ITER=50,
# MAX_RUNTIME_SECONDS=86400, NC=30, nlevel=1)
# Plot standard errors
cholQ <- chol(precision.fit$Q)
dim(cholQ) 
Phi <- precision.fit$Phi
dim(Phi)
marginal_variances <- rep(NA,dim(Phi)[1])
for(i in 1:dim(Phi)[1])
{
     marginal_variances[i] <- norm(backsolve(cholQ,Phi[i,],transpose=TRUE) ,type="F")^2
}
quilt.plot(tmin$lon.lat.proj,sqrt(marginal_variances), 
            main="Estimated process standard deviation")    
# A (noisy) simulation
c_coef <- backsolve(cholQ,rnorm(dim(Phi)[2])) 
sim <- Phi %*% c_coef + sqrt(precision.fit$nugget_variance)*rnorm(dim(Phi)[1])

grid_list <- list(
  x = seq(0, 1, length.out = 90),   # 4 x-axis breaks
  y = seq(0, 1, length.out = 90)    # 3 y-axis breaks
)

quilt.plot(tmin$lon.lat.proj,sim,main="Simulation")

install.packages('huge')
library(huge)
L = huge.generator(n=900, d=90, graph = "hub")
huge.plot(L$theta)

################################################################################
#Take range parameter = 1
library(LatticeKrig)
library(QUIC)
library(FRK)
library(psych) 
library(sp)
library(INLA)
library(glasso)

rho = 0.01*exp(-1)
l = 90
Sigma = matrix(0, nrow = l, ncol = l)
for(i in 1:l){
  for(j in 1:l){
    Sigma[i,j] = rho^abs(i-j)
  }
}
Q <- chol2inv(chol(Sigma))
cholQ <- chol(Q) 

x = seq(1,30,1)
locs = expand.grid(x,x)
d <- data.frame(locs)
G <- auto_basis(nres = 2, data = d,
                  type = "bisquare")
Phi <- as.matrix(eval_basis(G,as.matrix(locs)))

#the noise to signal ratio τ^2/(tr(τQ^−1T)/n) is fixed at 0.1
tao_sq <-  0.1*(tr(Phi%*%Sigma%*%t(Phi))/900) 
# nugget variance is fixed at .4225

# Generating m number of realizations on n=900 locations
sim_data <- function(n, m, cholQ, Phi, tao_sq){
  data <- matrix(0, ncol = m, nrow = n) 
  for(i in 1:m){
    c_coef <- backsolve(cholQ,rnorm(dim(Phi)[2]))
    sim <- Phi %*% c_coef + sqrt(tao_sq)*rnorm(dim(Phi)[1])
    data[,i] <- sim
  }
  return(data)
}
dat <- sim_data(dim(locs)[1], 1, cholQ, Phi, tao_sq)
quilt.plot(locs,dat,main="Simulation")


# Regression method
trial = 100
m = 5
rg_norm = numeric(trial)
for(i in 1:trial){
  dat <- sim_data(dim(locs)[1], m, cholQ, Phi, tao_sq)
  c_cap <- matrix(0, ncol = m, nrow = dim(Phi)[2]) 
  s <- matrix(0, ncol = dim(Phi)[2], nrow = dim(Phi)[2])
  for(j in 1:m){
    c_cap[,j] <- solve(t(Phi)%*%Phi)%*%t(Phi)%*%dat[,j]
    s <- s + c_cap[,j]%*%t(c_cap[,j]) 
  }
  gl <- glasso(s/m, rho = 0.25)
  rg_norm[i] <- norm(gl$wi - Q,type="F")/norm(Q,type="F")
}
avg_rgnorm <- mean(rg_norm)

# BGL Method
trial = 10
m = 5
fb_norm = numeric(trial)
for(i in 1:trial){
  dat <- sim_data(dim(locs)[1], m, cholQ, Phi, tao_sq)
  precision.fit <- BGL(y=dat, locs=locs, lambda=0.25, basis="LatticeKrig",
                       distance.penalty=TRUE,outer_tol=5e-2, MAX_ITER=50,
                       MAX_RUNTIME_SECONDS=86400, NC=4, NC.buffer = 3, nlevel=1)
  fb_norm[i] <- norm(precision.fit$Q - Q,type="F")/norm(Q,type="F")
}
avg_fbnorm <- mean(fb_norm)

# FRK method
trial = 10
m = 5
frk_norm = numeric(trial)
for(i in 1:trial){
  Phi_Q <- qr.Q(qr(Phi))
  Phi_R <- qr.R(qr(Phi))
  dat <- sim_data(dim(locs)[1], m, cholQ, Phi, tao_sq)
  s <- matrix(0, ncol = dim(Phi)[1], nrow = dim(Phi)[1])
  for(j in 1:m){
    s <- s + dat[,j]%*%t(dat[,j]) 
  }
  K <- solve(Phi_R)%*%t(Phi_Q)%*%((s/m) - tao_sq*diag(dim(Phi)[1]))%*%Phi_Q%*%t(solve(Phi_R))   
  frk_Q <- solve(K)
  frk_norm[i] <- norm(frk_Q - Q,type="F")/norm(Q,type="F")
}
avg_frknorm = mean(frk_norm)



##################### Timing study ######################
locs = data.frame(x = runif(2500), y = runif(2500))

sDomain <- apply(locs, 2, "range")

sim_data <- function(n, m, cholQ, Phi, tao_sq){
  data <- matrix(0, ncol = m, nrow = n) 
  for(i in 1:m){
    c_coef <- backsolve(cholQ,rnorm(dim(Phi)[2]))
    sim <- Phi %*% c_coef + sqrt(tao_sq)*rnorm(dim(Phi)[1])
    data[,i] <- sim
  }
  return(data)
}

set.seed(100)
NC<- seq(2,28,2)
obs <- c(50,100,250,500)
tm <- matrix(0, nrow = length(obs), ncol = length(NC))
for(i in 1:length(obs)){
  for(j in 1:length(NC)){
    LKinfo <- NULL
    LKinfo <- LKrigSetup(sDomain, a.wght=4.05, nu=0.5, NC = NC[j], NC.buffer = 3, nlevel = 1)
    Phi <- as.matrix(LKrig.basis(locs, LKinfo))
    l <- dim(Phi)[2]
    Q <- diag(l)
    cholQ <- chol(Q) 
    Sigma <- solve(Q)
    #the noise to signal ratio τ^2/(tr(τQ^−1T)/n) is fixed at 0.1
    tao_sq <-  0.1*(tr(Phi%*%Sigma%*%t(Phi))/2500)
    dat <- sim_data(dim(locs)[1], obs[i], cholQ, Phi, tao_sq)
    precision.fit <- NULL
    start <- Sys.time()
    precision.fit <- BGL(y=dat, locs=locs, lambda=0.25, basis="LatticeKrig",
        distance.penalty=TRUE,outer_tol=5e-2, MAX_ITER=50,
        MAX_RUNTIME_SECONDS=86400, NC= NC[j], NC.buffer = 3, nlevel=1)
    tm[i,j] <- Sys.time() - start
    print(i+j) 
  }
}


plot(NC, tm[1,], type = 'l', lwd = 2, col = 'red', xlab = 'Number of Basis Functions', ylab = 'Elapsed Time')
lines(NC, tm[2,], lwd = 2, col = 'blue')
lines(NC, tm[3,], lwd = 2, col = 'black')
lines(NC, tm[4,], lwd = 2, col = 'green')
legend("topleft", 
       legend = c("50", "100", "250", "500"),
       col = c("red", 'blue', 'black', 'green'),
       lwd = 2,
       title = "Number of Realizations")

#library(ggplot2)
#library(tidyr)
#m = matrix(seq(1,56,1) , ncol = 14)
#data = data.frame(NC = NC, obs1 = m[1,], obs2 = m[2,], obs3 = m[3,], obs4 = m[4,])

#df_long<- tidyr::gather(data, key = "Variable", value ="Value", -NC)

#ggplot(df_long, aes(x = NC, y = Value, color = Variable)) + 
#  geom_line() +
#  labs(x = "Response", y = "Value", title = "ABCD") +
#  scale_color_manual(values = c("obs1" = 'red','obs2' = 'blue','obs3' = 'green', 'obs4' = 'black'))+
#  theme_minimal()


# LKinfo = A list with components that give the information describing a multi-resolution basis with a Markov random field used for the covariance of the basis coefficients.


###############################################################################3
sim_data <- function(n, m, cholQ, Phi, tao_sq){
  data <- matrix(0, ncol = m, nrow = n) 
  for(i in 1:m){
    c_coef <- backsolve(cholQ,rnorm(dim(Phi)[2]))
    sim <- Phi %*% c_coef + sqrt(tao_sq)*rnorm(dim(Phi)[1])
    data[,i] <- sim
  }
  return(data)
}

set.seed(100)
n = c(1e4, 22500, 40000)
l = c(10.5, 15.5, 20.5) 
avg_fbnorm <- matrix(0, nrow = length(n), ncol = length(l)) 
avg_KLnorm <- matrix(0, nrow = length(n), ncol = length(l)) 
avg_percentage_MZ <- matrix(0, nrow = length(n), ncol = length(l))  
avg_percentage_NMZ <- matrix(0, nrow = length(n), ncol = length(l)) 
avg_diff_nugget_effect <- matrix(0, nrow = length(n), ncol = length(l)) 
avg_loglike_ratio <- matrix(0, nrow = length(n), ncol = length(l)) 
for(i in 1:length(n)){
  locs = as.matrix(data.frame(x = runif(n[i]), y = runif(n[i])))
  for(j in 1:length(l)){      
    LKinfo <- LKrigSetup(locs, a.wght=4.05, nu = 0.5, NC = l[j], NC.buffer = 0, nlevel = 1) #l = 100
    Q <- as.matrix(LKrig.precision(LKinfo))
    cholQ <- chol(Q)
    Sigma <- chol2inv(chol(Q))
    Phi <- as.matrix(LKrig.basis(locs, LKinfo))
    tao_sq <-  0.1*(tr(Phi%*%Sigma%*%t(Phi))/n[i]) 
    S_omit_neg.loglike <- log(det(Q + ((t(Phi)%*%Phi)/tao_sq))) - log(det(Q)) + n[i]*log(tao_sq)

    trial = 2
    m = 500
    fb_norm = numeric(trial)
    KL_norm = numeric(trial)
    percentage_MZ = numeric(trial)
    percentage_NMZ = numeric(trial)
    diff_nugget_effect = numeric(trial)
    loglike_ratio = numeric(trial) 
    true_neg_loglike = numeric(trial)
    est_neg_loglike = numeric(trial)
    dat = NULL
    precision.fit = NULL
    for(k in 1:trial){
      dat <- sim_data(dim(locs)[1], m, cholQ, Phi, tao_sq)
      CVguess <- BGL_CV(kfolds=5, y=dat, locs=locs, lambdalist=seq(0.005, 0.1, length = 8), basis="LatticeKrig",
                        outer_tol=0.05, MAX_ITER=50, MAX_RUNTIME_SECONDS=86400, NC=l[j], NC.buffer=0, nlevel=1, distance.penalty=TRUE)
      precision.fit <- BGL(y=dat, locs=locs, lambda=CVguess$best_lambda, basis="LatticeKrig",
                           distance.penalty=TRUE,outer_tol=5e-2, MAX_ITER=50,
                           MAX_RUNTIME_SECONDS=86400, NC=l[j], NC.buffer = 0, nlevel=1)
      fb_norm[k] <- norm(precision.fit$Q - Q,type="F")/norm(Q,type="F")
      KL_norm[k] <- tr(precision.fit$Q%*%Q) - log(det(precision.fit$Q%*%Q)) - dim(precision.fit$Q)[1]
      zeros_Qest <- precision.fit$Q == 0
      zeros_Q <- Q == 0
      missed_zeros <- zeros_Q & !zeros_Qest
      percentage_MZ[k] <- sum(missed_zeros) / sum(zeros_Q) * 100
      missed_nonzeros <- !zeros_Q & zeros_Qest
      percentage_NMZ[k] <- sum(missed_nonzeros) / sum(!zeros_Q) * 100
      diff_nugget_effect[k] <-  precision.fit$nugget_variance - tao_sq
      S <- matrix(0, ncol = dim(Phi)[1], nrow = dim(Phi)[1])
      for(a in 1:m){
        S <- S + dat[,a]%*%t(dat[,a]) 
      }
      S <- S/m
      true_neg_loglike[k] <- S_omit_neg.loglike - tr((t(Phi)%*%S%*%Phi%*%solve(Q + ((t(Phi)%*%Phi)/tao_sq)))/(tao_sq^2)) + tr(S)/tao_sq
      est_neg_loglike[k] <-  log(det(precision.fit$Q + ((t(Phi)%*%Phi)/precision.fit$nugget_variance))) - log(det(precision.fit$Q)) - tr((t(Phi)%*%S%*%Phi%*%solve(precision.fit$Q + ((t(Phi)%*%Phi)/precision.fit$nugget_variance)))/(precision.fit$nugget_variance^2)) + n[i]*log(precision.fit$nugget_variance) + tr(S)/precision.fit$nugget_variance
      loglike_ratio[k] <- est_neg_loglike[k]/true_neg_loglike[k]
    }
    avg_fbnorm[i,j] <- mean(fb_norm)
    avg_KLnorm[i,j] <- mean(KL_norm)
    avg_percentage_MZ[i,j] <- mean(percentage_MZ) 
    avg_percentage_NMZ[i,j] <- mean(percentage_NMZ)
    avg_diff_nugget_effect[i,j] <- mean(diff_nugget_effect)
    avg_loglike_ratio[i,j] <- mean(loglike_ratio)
    print(n[i]+l[j])
  }
}

my_list = list(avg_fbnorm, avg_KLnorm, avg_percentage_MZ, avg_percentage_NMZ, avg_diff_nugget_effect, avg_loglike_ratio) 
save(my_list, file = "single_level_local_basis.Rdata")





