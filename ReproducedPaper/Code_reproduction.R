library(BasisGraphicalLasso)  
library(QUIC)
#####Despcription#####
#This package contains code to perform the basis graphical lasso analysis of:
#Krock, M., Kleiber, W., and Becker, S. (2021), “Nonstationary modeling with sparsity 
#for spatial data via the basis graphical lasso,”
library(sp)
library(LatticeKrig)
#####Despcription#####
#Methods for the interpolation of large spatial datasets. This package follows a ``fixed rank Kriging'' approach but
#provides a surface fitting method that can approximate standard spatial data models.
library(FRK)
library(INLA)
#####Despcription#####
#A tool for spatial/spatio-temporal modelling and prediction with large datasets. The approach
#models the field, and hence the covariance function, using a set of basis functions.
library(psych) 
#####Despcription#####
#used to determine trace of a matrix
library(glasso)
# Estimation of a sparse inverse covariance matrix using a lasso (L1) penalty.
library(huge)
library(sf)
####################### Code to generate a realization from the basis model.#########################
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
# nugget variance is tao_sq

# Generating m number of realizations on n locations
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
dat <- sim_data(dim(locs)[1], 1, cholQ, Phi, tao_sq)
quilt.plot(locs,dat,main="Simulation")

# Graph structure
L = huge.generator(n=900, d=90, graph = "hub")
huge.plot(L$theta)

#################### Comparison with other methods #############################
set.seed(221282)
# Regression method
trial = 5
m = 1
rg_norm = numeric(trial)
KL_norm = numeric(trial) 
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
  KL_norm[i] <- tr(gl$wi%*%Q) - log(det(gl$wi%*%Q)) - dim(gl$wi)[1]
}
avg_rgnorm <- mean(rg_norm)
avg_klnorm <- mean(KL_norm)

# BGL Method
trial = 10
m = 1
fb_norm = numeric(trial)
KL_norm = numeric(trial) 
for(i in 1:trial){
  dat <- sim_data(dim(locs)[1], m, cholQ, Phi, tao_sq)
  precision.fit <- BGL(y=dat, locs=locs, lambda=0.25, basis="LatticeKrig",
                       distance.penalty=TRUE,outer_tol=5e-2, MAX_ITER=50,
                       MAX_RUNTIME_SECONDS=86400, NC=2, NC.buffer = 1, nlevel=3)
  fb_norm[i] <- norm(precision.fit$Q - Q,type="F")/norm(Q,type="F")
  KL_norm[i] <- tr(precision.fit$Q%*%Q) - log(det(precision.fit$Q%*%Q)) - dim(precision.fit$Q)[1]
}
avg_fbnorm <- mean(fb_norm)
avg_klnorm <- mean(KL_norm)

# FRK method
# For m = 1 use the following code
trial = 10
m = 1
frk_norm = numeric(trial)
KL_norm = numeric(trial) 
if(m == 1){ 
  for(i in 1:trial){
      dat <- sim_data(dim(locs)[1], m, cholQ, Phi, tao_sq)
      df <- data.frame(locs, dat)  
      coordinates(df) = ~Var1+Var2
      basis <- auto_basis(manifold = plane(), data = df, nres = 2)
      frk.fit <- FRK(f = as.vector(dat)~1, data = df, basis = basis, K_type = 'unstructured')
      frk_norm[i] <- norm(frk.fit@Q_eta - Q,type="F")/norm(Q,type="F")
      KL_norm[i] <- tr(precision.fit$Q%*%Q) - log(det(precision.fit$Q%*%Q)) - dim(precision.fit$Q)[1]
  }
  avg_frknorm = mean(frk_norm)
  avg_klnorm <- mean(KL_norm)
}else if(m > 1){
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
    KL_norm[i] <- tr(precision.fit$Q%*%Q) - log(det(precision.fit$Q%*%Q)) - dim(precision.fit$Q)[1]
  }
  avg_frknorm = mean(frk_norm)
  avg_klnorm <- mean(KL_norm)
}

###################### Timing Study for simulation #############################
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
NC <- seq(2,28,2)
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
save(tm, file = "Timing_study.Rdata")

plot(NC, tm[1,], type = 'l', lwd = 2, col = 'red', xlab = 'Number of Basis Functions', ylab = 'Elapsed Time')
lines(NC, tm[2,], lwd = 2, col = 'blue')
lines(NC, tm[3,], lwd = 2, col = 'black')
lines(NC, tm[4,], lwd = 2, col = 'green')
legend("topleft", 
       legend = c("50", "100", "250", "500"),
       col = c("red", 'blue', 'black', 'green'),
       lwd = 2,
       title = "Number of Realizations")

##################### Local Basis: Single level case ###########################
# Please run the following code carefully. It may take 2 days to complete the
# execution. So, it is preferable to run the part inside the for loop individually 
# for a particular value of n and l, to match the given output.
set.seed(100)
n = c(1e4, 22500, 40000)
l = c(10, 15, 20) 
avg_fbnorm <- matrix(0, nrow = length(n), ncol = length(l)) 
avg_KLnorm <- matrix(0, nrow = length(n), ncol = length(l)) 
avg_percentage_MZ <- matrix(0, nrow = length(n), ncol = length(l))  
avg_percentage_NMZ <- matrix(0, nrow = length(n), ncol = length(l)) 
avg_diff_nugget_effect <- matrix(0, nrow = length(n), ncol = length(l)) 
avg_loglike_ratio <- matrix(0, nrow = length(n), ncol = length(l)) 
for(i in 1:length(n)){
  locs = as.matrix(data.frame(x = runif(n[i]), y = runif(n[i])))
  for(j in 1:length(l)){      
    LKinfo <- LKrigSetup(locs, a.wght=4.05, nu = 0.5, NC = l[j], NC.buffer = 0.25, nlevel = 1) #l = 100
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
                        outer_tol=0.05, MAX_ITER=50, MAX_RUNTIME_SECONDS=86400, NC=l[j], NC.buffer=0.25, nlevel=1, distance.penalty=TRUE)
      precision.fit <- BGL(y=dat, locs=locs, lambda=CVguess$best_lambda, basis="LatticeKrig",
                           distance.penalty=TRUE,outer_tol=5e-2, MAX_ITER=50,
                           MAX_RUNTIME_SECONDS=86400, NC=l[j], NC.buffer = 0.25, nlevel=1)
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

my_list1 = list(avg_fbnorm, avg_KLnorm, avg_percentage_MZ, avg_percentage_NMZ, avg_diff_nugget_effect, avg_loglike_ratio) 
save(my_list1, file = "single_level_local_basis.Rdata")

##################### Local Basis: Multiple level case ###########################
# Please run the following code carefully. It may take 2 days to complete the
# execution. So, it is preferable to run the part inside the for loop individually 
# for a particular value of n and NC & nlevel, to match the given output.
# NC = For regular grids of lattice points the maximum number of lattice grid points 
#     for a spatial coordinate and at the coarsest level of resolution. 
#     For a example, for a square region, (and LKGeometry = "LKRectangle") NC=5 results 
#      in a 5X5 = 25 lattice points at the first level. Note that the default is that 
#      lattice points are also taken to have the same spacing in every dimension.
# nlevel = Number of levels in multi-resolution. Note that each subsequent level increases 
#          number of basis functions within the spatial domain size by a factor of roughly 4.
set.seed(100)
n = c(1e4, 22500, 40000)
NC = c(2,4,3)
nres = c(4,3,4) 
avg_fbnorm <- matrix(0, nrow = length(n), ncol = length(l)) 
avg_KLnorm <- matrix(0, nrow = length(n), ncol = length(l)) 
avg_percentage_MZ <- matrix(0, nrow = length(n), ncol = length(l))  
avg_percentage_NMZ <- matrix(0, nrow = length(n), ncol = length(l)) 
avg_diff_nugget_effect <- matrix(0, nrow = length(n), ncol = length(l)) 
avg_loglike_ratio <- matrix(0, nrow = length(n), ncol = length(l)) 
for(i in 1:length(n)){
  locs = as.matrix(data.frame(x = runif(n[i]), y = runif(n[i])))
  for(j in 1:length(l)){      
    LKinfo <- LKrigSetup(locs, a.wght=4.05, nu = 0.5, NC = NC[j], NC.buffer = 0.25, nlevel = nres[j]) #l = 100
    Q <- as.matrix(LKrig.precision(LKinfo))
    dim(Q)
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
                        outer_tol=0.05, MAX_ITER=50, MAX_RUNTIME_SECONDS=86400, NC=NC[j], NC.buffer=0.25, nlevel=nres[j], distance.penalty=TRUE)
      precision.fit <- BGL(y=dat, locs=locs, lambda=CVguess$best_lambda, basis="LatticeKrig",
                           distance.penalty=TRUE,outer_tol=5e-2, MAX_ITER=50,
                           MAX_RUNTIME_SECONDS=86400, NC=NC[j], NC.buffer = 0.25, nlevel=nres[j])
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
    print(n[i]+NC[j])
  }
}

my_list2 = list(avg_fbnorm, avg_KLnorm, avg_percentage_MZ, avg_percentage_NMZ, avg_diff_nugget_effect, avg_loglike_ratio) 
save(my_list2, file = "multi_level_local_basis.Rdata")

#################################  Application  ########################################

load('tmin.Rdata') 
locs <- tmin$lon.lat 
LKinfo <- LKrigSetup(locs, a.wght=4.05, nu = 0.5, NC = 30, nlevel = 1)
look<- LKrigLatticeCenters(LKinfo, Level = 1)
look<- make.surface.grid(look)
look<- as.data.frame(look)
plot(look, cex = .5)
par(new = T)
points(tmin$lon.lat[, 1],tmin$lon.lat[, 2], col = tmin$data[, 1]+50, pch = 16, xaxt = 'n', yaxt = 'n')
############ Couldn't generate the plot in the paper

########### Performed BGL
precision.fit <- BGL(y=tmin$data, locs=tmin$lon.lat.proj, lambda=7, basis="LatticeKrig",
                     distance.penalty=TRUE,outer_tol=1e-2, MAX_ITER=50,
                     MAX_RUNTIME_SECONDS=86400, NC=30, nlevel=1)

# Plot standard errors
cholQ <- chol(precision.fit$Q)
Phi <- precision.fit$Phi
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
quilt.plot(tmin$lon.lat.proj,sim,main="Simulation")





