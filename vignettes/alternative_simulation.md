alternative_simulation
================

``` r
library(irregMFPCA)
library(tidyverse)
library(fda)
library(fdapace)
library(MFPCA)
set.seed(17)
```

``` r
rm(list = ls())
```

## Function to Normalize Eigenscores

``` r
normalizer = function(mat, normal = NULL){
  if(is.null(normal)){normal = nrow(mat)}
  
  lmat = as.list(data.frame(mat))
  constants = sqrt(normal/sapply(lmat, function(x) sum(x^2)))
  
  return(matrix(c(constants[1]*mat[, 1],
                  constants[2]*mat[, 2],
                  constants[3]*mat[, 3]), 
                nrow = nrow(mat),
                ncol = ncol(mat)))
}
```

## Generate Mean Functions

$$
\mu_{X_1}(t) = -2(t-0.5)^2 + 5 \text{ for }t \in [0,1]
$$

$$
\mu_{X_2}(t) = 3(t-0.75)^3 -5 \text{ for }t \in [0,1]
$$

``` r
TT = 100
t = seq(0, 1, length.out = TT)
mu_X1 = -2*(t-0.5)^2 + 5
mu_X2 = 3*(t-0.75)^3 - 5
```

## Generate Eigenfunctions and Eigenscores

$$
\phi^{X_1}_1(t) = \sqrt{2}\sin(\pi t), \ \phi^{X_1}_2(t) = \sqrt{2}\sin(7\pi t), \ \phi^{X_1}_3(t) = \sqrt{2}\cos(7\pi t) 
$$

$$
\phi^{X_2}_1(t) = \sqrt{2}\cos(5\pi t), \ \phi^{X_2}_2(t) = \sqrt{2}\cos(\pi t), \ \phi^{X_2}_3(t) = \sqrt{2}\sin(7\pi t) 
$$

``` r
n = 199;
lambda11 = 5; lambda21 = 3; lambda31 = 1
D = diag(c(lambda11, lambda21, lambda31))
xi_X = matrix(0, nrow = n, ncol = 3)
xi_X[,1] = rnorm(n, 0, 1)
xi_X[,2] = rnorm(n, 0, 1)
xi_X[,3] = rnorm(n, 0, 1)

xi_X = apply(xi_X, 2, scale)
xi_X = normalizer(xi_X)

phi_X1 = matrix(0, nrow = TT, ncol = 3)
phi_X1[,1] = sin(2*pi*t)
phi_X1[,2] = sin(4*pi*t)
phi_X1[,3] = sin(6*pi*t)

phi_X1 = normalizer(phi_X1, normal = TT)

phi_X2 = matrix(0, nrow = TT, ncol = 3)
phi_X2[,1] = cos(3*pi*t)
phi_X2[,2] = cos(pi*t)
phi_X2[,3] = cos(5*pi*t)

phi_X2 = normalizer(phi_X2, normal = TT)

phi_X = rbind(phi_X1, phi_X2)

## Check for orthogonality

tol = 1e-1
check = 1/TT*t(phi_X) %*% phi_X
if(check[1, 2] > tol | check[1, 2] < -tol){
  cat("The pair of functions were not orthogonal. \nPlease choose a different pair.\n")
}

## Resplit into phi_X1 and phi_X2
phi_X1 = phi_X[1:TT,]
phi_X2 = phi_X[(TT+1):(2*TT),]

X = xi_X %*% D %*% t(phi_X) + rep(1, n) %*% t(c(mu_X1, mu_X2))

# X1 = xi_X %*% D %*% t(phi_X1) + rep(1, n) %*% t(mu_X1)
# X2 = xi_X %*% D %*% t(phi_X2) + rep(1, n) %*% t(mu_X2)

phi_X %>%
  as.data.frame() %>%
  ggplot(aes(x = 1:200)) +
  geom_line(aes(y = V1), col = "red") +
  geom_line(aes(y = V2), col = "blue") +
  geom_line(aes(y = V3), col = "green") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-335-1.png)<!-- -->

``` r

plot(mu_X1)
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-335-2.png)<!-- -->

``` r


plot(mu_X2)
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-335-3.png)<!-- -->

``` r
# V = normalizer(phi_X1)
# 
# t(V) %*% V / 100
```

## Generate Response

$$
\mu_{Y_1}(t) = 6\exp(-(t-1)^2)
$$

$$
\mu_{Y_2}(t) = -2\cdot 14^{(t-0.5)}
$$

$$
\phi^{Y_1}_1(t) = \sqrt{2}\cos(9\pi t), \ \phi^{Y_1}_2(t) = \sqrt{2}\sin(5\pi t), \ \phi^{Y_1}_3(t) = \sqrt{2}\cos(2\pi t)
$$

$$
\phi^{Y_2}_1(t) = \sqrt{2}\sin(3\pi t), \ \phi^{Y_2}_2(t) = \sqrt{2}\cos(\pi t), \ \phi^{Y_2}_3(t) = \sqrt{2}\sin(7\pi t)
$$

``` r
mu_Y1 = 6*exp(-(t-1)^2)
mu_Y2 = -2*14^(t-0.5)

phi_Y1 = matrix(0, nrow = TT, ncol = 3)
phi_Y1[,1] = cos(9*pi*t)
phi_Y1[,2] = cos(5*pi*t)
phi_Y1[,3] = cos(2*pi*t)
phi_Y1 = normalizer(phi_Y1, normal = TT/2)

phi_Y2 = matrix(0, nrow = TT, ncol = 3)
phi_Y2[,1] = sin(3*pi*t)
phi_Y2[,2] = sin(5*pi*t)
phi_Y2[,3] = sin(7*pi*t)
phi_Y2 = normalizer(phi_Y2, normal = TT/2)

phi_Y = rbind(phi_Y1, phi_Y2)

check = 1/TT*t(phi_Y) %*% phi_Y
if(check[1, 2] > tol | check[1, 2] < -tol){
  cat("The pair of functions were not orthogonal. \nPlease choose a different pair.\n")
}

phi_Y1 = phi_Y[1:TT,]
phi_Y2 = phi_Y[(TT+1):(2*TT),]

phi_Y %>%
  as.data.frame() %>%
  ggplot(aes(x = 1:200)) +
  geom_line(aes(y = V1), col = "red") +
  geom_line(aes(y = V2), col = "blue") +
  geom_line(aes(y = V3), col = "green") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-337-1.png)<!-- -->

``` r

plot(mu_Y1)
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-337-2.png)<!-- -->

``` r

# phi_Y2 %>%
#   as.data.frame() %>%
#   ggplot(aes(x = t)) +
#   geom_line(aes(y = V1), col = "red") +
#   geom_line(aes(y = V2), col = "blue") +
#   geom_line(aes(y = V3), col = "green") +
#   theme_bw()

plot(mu_Y2)
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-337-3.png)<!-- -->

## Combine Components

``` r
# X = cbind(X1, X2)
# phi_X = cbind(phi_X1, phi_X2)
# xi_X = cbind(xi_X1, xi_X2)

# phi_Y = cbind(phi_Y1, phi_Y2)
# xi_Y = cbind(xi_Y1, xi_Y2)
mu_Y = c(mu_Y1, mu_Y2)
```

## Noise

``` r
sigma = 0.001
n_eig_X = 3
n_eig_Y = 3
```

## Regression

``` r
# Generate B matrix from Unif(-U, U)
# U = 3
# B = matrix(sample(seq(-U, U, 0.01),n_eig_X*n_eig_Y),
#            nrow = n_eig_Y,
#            ncol = n_eig_X)

# sample(-9:9, 9, replace = T)
B = matrix(c(-1, 2, 3, -1, -9, 3, 5, 5, -3),
           nrow = n_eig_Y,
           ncol = n_eig_X)

# epsilon = matrix(rnorm(n_eig_Y*n, mean = 0, sd = 0.01), 
#                  n_eig_Y, 
#                  n)

Y1 = matrix(nrow = n, 
            ncol = TT)
Y2 = matrix(nrow = n, 
            ncol = TT)
```

### Introduce Outliers (DISCARDED FOR NOW)

``` r
# outliers = 20
# non_outliers = n - outliers
# 

# This is equal to xi_X1 %*% t(B) + mu_Y1

# for (i in 1:n){
#   Y1[i,] = t(B%*%xi_X[i,])%*%t(phi_Y1) + mu_Y1 # + t(phi_Y1%*%epsilon[,i])
#   Y2[i,] = t(B%*%xi_X[i,])%*%t(phi_Y2) + mu_Y2 # + t(phi_Y2%*%epsilon[,i])
# }

# Y1 = xi_X %*% B %*% D %*% t(phi_Y1) + rep(1, n) %*% t(mu_Y1)
# Y2 = xi_X %*% B %*% D %*% t(phi_Y2) + rep(1, n) %*% t(mu_Y2)

Y = xi_X %*% B %*% D %*% t(phi_Y) + rep(1, n) %*% t(c(mu_Y1, mu_Y2))
```

## Scenarios

``` r
# scenario = 1
# 
# if (scenario==1){
#   outlier_var = 0.5
#   B_out = B + matrix(rnorm(n_eig_X*n_eig_Y, mean = 0, sd = outlier_var),
#                      nrow = n_eig_Y,
#                      ncol = n_eig_X)
#   
#   for (i in (non_outliers+1):n){
#     Y1[i,] = t(B_out%*%xi_X1[i,])%*%t(phi_Y1) + mu_Y1 + t(phi_Y1%*%epsilon[,i])
#     Y2[i,] = t(B_out%*%xi_X1[i,])%*%t(phi_Y2) + mu_Y2 + t(phi_Y2%*%epsilon[,i])
#   }
# }

# if (scenario==2){
#   gaittime = seq(0, 1, len=TT)
#   gaitrange = c(0,1)
#   nord = 10
#   gaitbasis = create.bspline.basis(c(0,1), nbasis = nord, norder = 3)
#   evaluated_basis = eval.basis(gaitbasis, gaittime)
#   normalised_basis = eval.basis(gaitbasis, gaittime)
#   
#   values = diag(inprod(gaitbasis,gaitbasis)*TT)
#   intu = sample(nord, 1)
#   
#   out_basis = normalised_basis[, intu]
#   
#   B_out = rbind(B, rnorm(3,3,1))
#   
#   for (i in (non_outliers+1):n){
#     Y[i,] = t(cbind(phi_Y, out_basis/values[intu]) %*% (B_out %*% xi_X[i,]) + mu_Y) + t(phi_Y%*%epsilon[,i])
#   }
# }
```

``` r
# Y = cbind(Y1, Y2)

E = matrix(rnorm(2*TT*n, mean = 0, sd = sigma), n, 2*TT)

Y =  Y + E
```

``` r
#### plot Predictior and Response curves ####

matplot(t(X), 
        type='l', 
        ylab='X(t)', 
        xlab='time', 
        main='Plot of predictor curves', 
        col=rgb(0,0,0,alpha=0.4))
matlines(apply(t(X), 1, mean),
         type='l',
         lwd=3,
         lty=1,
         col="red")
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-344-1.png)<!-- -->

``` r

matplot(t(Y), 
        type='l', 
        ylab='Y(t)', 
        xlab='time', 
        main='Plot of response curves', 
        col=rgb(0,0,0,alpha=0.6))
matlines(apply(t(Y), 1, mean),
         type='l',
         lwd=3,
         lty=1)
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-344-2.png)<!-- -->

### MFPCA

``` r
X_miss = irregMFPCA::missing_data_simulator(n_obs = 99,
                                            dat = X,
                                            t = seq(0, 1, length.out = TT),
                                            seed = 16)


df = X_miss %>% tibble_format(tindex = t) %>% fpca_format()
#> New names:
#> • `v` -> `v...3`
#> • `v` -> `v...4`
```

``` r
res1 = FPCA(df$Component1,
            df$Time,
            list(dataType='Sparse',
                 error=FALSE,
                 kernel='epan',
                 verbose=TRUE,
                 nRegGrid = 100))
#> No binning is needed!
#> At most 48 number of PC can be selected, thresholded by `maxK` = 20.

res2 = FPCA(df$Component2,
            df$Time,
            list(dataType='Sparse',
                 error=FALSE,
                 kernel='epan',
                 verbose=TRUE,
                 nRegGrid = 100))
#> No binning is needed!
#> At most 59 number of PC can be selected, thresholded by `maxK` = 20.
```

``` r
# plot(res)
```

### Compare to Actual: First Component

``` r
act = data.frame(act1 = mu_X1,
                 act2 = mu_X2)
hat = data.frame(hat1 = res1$mu,
                 hat2 = res2$mu)

if(t(res1$phi[, 1]) %*% phi_X1[, 1] < 0){
  res1$phi[, 1] = -res1$phi[, 1]
}

if(t(res1$phi[, 2]) %*% phi_X1[, 2] < 0){
  res1$phi[, 2] = -res1$phi[, 2]
}

if(t(res1$phi[, 3]) %*% phi_X1[, 3] < 0){
  res1$phi[, 3] = -res1$phi[, 3]
}

hat %>%
  ggplot() +
  geom_line(aes(x = 1:100, y = hat1)) +
  geom_line(data = act, aes(x = seq(1, 100, length.out = TT), y = act1), linetype = "dashed") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-348-1.png)<!-- -->

``` r

phi_X1_df = phi_X1 %>% as.data.frame()

res1$phi[, 1:3] %>%
  as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = 1:100, y = V1), col = "red") +
  geom_line(data = phi_X1_df, aes(x = 1:100, y = V1), col = "red", linetype = "dashed") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-348-2.png)<!-- -->

``` r

res1$phi[, 1:3] %>%
  as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = 1:100, y = V2), col = "blue") +
  geom_line(data = phi_X1_df, aes(x = seq(1, 100, length.out = TT), y = V2), col = "blue", linetype = "dashed") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-348-3.png)<!-- -->

``` r

res1$phi[, 1:3] %>%
  as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = 1:100, y = V3), col = "green") +
  geom_line(data = phi_X1_df, aes(x = seq(1, 100, length.out = TT), y = V3), col = "green", linetype = "dashed") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-348-4.png)<!-- -->

### Compare to Actual: Second Component

``` r
hat %>%
  ggplot() +
  geom_line(aes(x = 1:100, y = hat2)) +
  geom_line(data = act, aes(x = seq(1, 100, length.out = TT), y = act2), linetype = "dashed") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-349-1.png)<!-- -->

``` r

phi_X2_df = phi_X2 %>% as.data.frame()

if(t(res2$phi[, 1]) %*% phi_X2[, 1] < 0){
  res2$phi[, 1] = -res2$phi[, 1]
}

if(t(res2$phi[, 2]) %*% phi_X2[, 2] < 0){
  res2$phi[, 2] = -res2$phi[, 2]
}

if(t(res2$phi[, 3]) %*% phi_X2[, 3] < 0){
  res2$phi[, 3] = -res2$phi[, 3]
}

res2$phi[, 1:3] %>%
  as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = 1:100, y = V1), col = "red") +
  geom_line(data = phi_X2_df, aes(x = seq(1, 100, length.out = TT), y = V1), col = "red", linetype = "dashed") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-349-2.png)<!-- -->

``` r

res2$phi[, 1:3] %>%
  as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = 1:100, y = V2), col = "blue") +
  geom_line(data = phi_X2_df, aes(x = seq(1, 100, length.out = TT), y = V2), col = "blue", linetype = "dashed") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-349-3.png)<!-- -->

``` r

res2$phi[, 1:3] %>%
  as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = 1:100, y = V3), col = "green") +
  geom_line(data = phi_X2_df, aes(x = seq(1, 100, length.out = TT), y = V3), col = "green", linetype = "dashed") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-349-4.png)<!-- -->

## Irregular MFPCA

``` r
mfpca_res = irregMFPCA::irregMFPCA(components = 3, #TODO: Ambiguous parameter name. Number of principal components to consider
                                   split = T,
                                   res1,
                                   res2)
# TODO: Check normalizing function and nRegGrid
```

### Check Results

``` r
mfpca_eigenf = mfpca_res$unstackpsi
colnames(mfpca_eigenf) = c("var_1_1", "var_1_2", "var_1_3", "var_2_1", "var_2_2", "var_2_3")
mfpca_eigens = mfpca_res$rho
colnames(mfpca_eigens) = c("function_1", "function_2", "function_3")

normalized_eigens = mfpca_eigens %*% solve(mfpca_res$Dhat)
```

``` r
tplot = seq(0, 1, length.out = 100)
mfpca_eigenf %>%
  as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = tplot, y = var_1_1), col = "red") +
  geom_line(data = phi_X1_df, aes(x = seq(0, 1, length.out = 100), y =V1), col = "black", linetype = "dotted") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-352-1.png)<!-- -->

``` r

mfpca_eigenf %>%
  as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = tplot, y = -var_1_2), col = "red") +
  geom_line(data = phi_X1_df, aes(x = seq(0, 1, length.out = 100), y =V2), col = "black", linetype = "dotted") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-352-2.png)<!-- -->

``` r

mfpca_eigenf %>%
  as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = tplot, y = -var_1_3), col = "red") +
  geom_line(data = phi_X1_df, aes(x = seq(0, 1, length.out = 100), y =V3), col = "black", linetype = "dotted") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-352-3.png)<!-- -->

``` r


mfpca_eigenf %>%
  as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = tplot, y = var_2_1), col = "blue") +
  geom_line(data = phi_X2_df, aes(x = seq(0, 1, length.out = 100), y =V1), col = "black", linetype = "dotted") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-352-4.png)<!-- -->

``` r

mfpca_eigenf %>%
  as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = tplot, y = -var_2_2), col = "blue") +
  geom_line(data = phi_X2_df, aes(x = seq(0, 1, length.out = 100), y =V2), col = "black", linetype = "dotted") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-352-5.png)<!-- -->

``` r

mfpca_eigenf %>%
  as.data.frame() %>%
  ggplot() +
  geom_line(aes(x = tplot, y = -var_2_3), col = "blue") +
  geom_line(data = phi_X2_df, aes(x = seq(0, 1, length.out = 100), y =V3), col = "black", linetype = "dotted") +
  theme_bw()
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-352-6.png)<!-- -->

``` r
eigens = data.frame(est1 = -mfpca_eigens[,1]/5,
                    est2 = mfpca_eigens[,2]/3,
                    est3 = mfpca_eigens[,3],
                    act1 = xi_X[, 1],
                    act2 = xi_X[, 2],
                    act3 = xi_X[, 3])

eigens %>%
  ggplot(aes(x = est1, y = act1)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  geom_smooth(method = "lm", se = F) +
  theme_bw()
#> `geom_smooth()` using formula = 'y ~ x'
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-353-1.png)<!-- -->

``` r

eigens %>%
  ggplot(aes(x = est2, y = act2)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  geom_smooth(method = "lm", se = F) +
  theme_bw()
#> `geom_smooth()` using formula = 'y ~ x'
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-353-2.png)<!-- -->

``` r

eigens %>%
  ggplot(aes(x = est3, y = act3)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  geom_smooth(method = "lm", se = F) +
  theme_bw()
#> `geom_smooth()` using formula = 'y ~ x'
```

![](alternative_simulation_files/figure-gfm/unnamed-chunk-353-3.png)<!-- -->

## Prepare for Regression

### MFPCA for Y

``` r
Y_miss = irregMFPCA::missing_data_simulator(n_obs = 99,
                                            dat = Y,
                                            t = seq(0, 1, length.out = TT),
                                            seed = 16)
df_y = Y_miss %>% tibble_format(tindex = seq(0, 1, length.out = TT)) %>% fpca_format()
#> New names:
#> • `v` -> `v...3`
#> • `v` -> `v...4`

res1_y = FPCA(df_y$Component1,
            df_y$Time,
            list(dataType='Sparse',
                 error=FALSE,
                 kernel='epan',
                 verbose=TRUE,
                 nRegGrid = 100))
#> No binning is needed!
#> 
#> At most 58 number of PC can be selected, thresholded by `maxK` = 20.

res2_y = FPCA(df_y$Component2,
            df_y$Time,
            list(dataType='Sparse',
                 error=FALSE,
                 kernel='epan',
                 verbose=TRUE,
                 nRegGrid = 100))
#> No binning is needed!
#> 
#> At most 53 number of PC can be selected, thresholded by `maxK` = 20.

mfpca_resy = irregMFPCA::irregMFPCA(components = 3,
                                    split = F,
                                    res1_y,
                                    res2_y)
```

### Multivariate Multiple Regression

``` r
response = mfpca_resy$rho
predictor = mfpca_res$rho

mod = lm(response ~ -1 + predictor)

B; mod$coefficients
#>      [,1] [,2] [,3]
#> [1,]   -1   -1    5
#> [2,]    2   -9    5
#> [3,]    3    3   -3
#>                  [,1]       [,2]        [,3]
#> predictor1 -0.7817527  0.7300717  0.34046584
#> predictor2  6.8925923  1.0337879  0.06637909
#> predictor3 -4.3358331 11.4329957 -0.60267123
```
