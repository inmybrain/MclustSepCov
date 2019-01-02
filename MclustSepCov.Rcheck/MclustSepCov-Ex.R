pkgname <- "MclustSepCov"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('MclustSepCov')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Mclust_SEP_cpp")
### * Mclust_SEP_cpp

flush(stderr()); flush(stdout())

### Name: Mclust_SEP_cpp
### Title: The model-based clustering for longitudinal data
### Aliases: Mclust_SEP_cpp

### ** Examples

# Gaussian mixture model with two components
K <- 2
p <- 2
q <- 3
U <- lapply(1:K, function(noarg) getCovariance(p, 0.3, "AR"))
V <- lapply(1:K, function(noarg) getCovariance(q, 0.2, "CS"))
Sigma <- Map(kronecker, U, V) # separable covariance matrix
mu <- list(rep(0, p * q), 5 / sqrt(p*q) * rep(1, p * q)) # distinct mean vectors
Y <- vector(mode = "list", length = K)
for(i in 1:K){
  Y[[i]] <- mvtnorm::rmvnorm(n = 20, mean = mu[[i]], sigma = Sigma[[i]])
}
fit <- Mclust_SEP_cpp(Y = Reduce(rbind, Y), p = p, q = q, Ks = 2, type_cov = "EEE-ECS")



cleanEx()
nameEx("Mclust_SEP_each_cpp")
### * Mclust_SEP_each_cpp

flush(stderr()); flush(stdout())

### Name: Mclust_SEP_each_cpp
### Title: The maximum likelihood estimation of the mixture distribution
### Aliases: Mclust_SEP_each_cpp

### ** Examples

# Gaussian mixture model with two components
K <- 2
p <- 2
q <- 3
U <- lapply(1:K, function(noarg) getCovariance(p, 0.3, "AR"))
V <- lapply(1:K, function(noarg) getCovariance(q, 0.2, "CS"))
Sigma <- Map(kronecker, U, V) # separable covariance matrix
mu <- list(rep(0, p * q), 5 / sqrt(p*q) * rep(1, p * q)) # distinct mean vectors
Y <- vector(mode = "list", length = K)
for(i in 1:K){
  Y[[i]] <- mvtnorm::rmvnorm(n = 20, mean = mu[[i]], sigma = Sigma[[i]])
}
fit <- Mclust_SEP_each_cpp(Y = Reduce(rbind, Y), p = p, q = q, K = 2, type_vari = "EEE", type_time = "ECS")



cleanEx()
nameEx("Optimization")
### * Optimization

flush(stderr()); flush(stdout())

### Name: Optimization
### Title: Newton-Raphson's algorithm to find the optimal temporal
###   correlation
### Aliases: Optimization LB_algorithm_cpp

### ** Examples

q <- 10
# AR model
set.seed(6)
Y <- mvtnorm::rmvnorm(100, rep(0, q), getCovariance(q, 0.3, "AR"))
LB_algorithm_cpp(a = nrow(Y), Z = Y, rho0 = 1e-3, type = "AR")

# CS model
set.seed(6)
Y <- mvtnorm::rmvnorm(100, rep(0, q), getCovariance(q, 0.3, "CS"))
LB_algorithm_cpp(a = nrow(Y), Z = Y, rho0 = 1e-3, type = "CS")



cleanEx()
nameEx("getCovariance")
### * getCovariance

flush(stderr()); flush(stdout())

### Name: getCovariance
### Title: Generate temporal covariance matrices
### Aliases: getCovariance

### ** Examples

getCovariance(3, 0.3, "AR") # AR
getCovariance(3, 0.3, "CS") # CS

# AR structure with heterogeneous variances
hvar <- c(1, 2, 3) # variances
diag(sqrt(hvar)) %*% getCovariance(3, 0.3, "AR") %*% diag(sqrt(hvar))



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
