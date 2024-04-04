#' Internal function to generate orthogonal eigenfunctions
#'
#' @param orthfuns choice of orthogonal functions
#' @param periods choice of periods to adjust orthogonal functions by
#' @param t time vector
#'
#' @return A pair of orthogonal vectors
#' @noRd
#' @importFrom dplyr %>%
#'
genOrthPair = function(orthfuns = c(sin, cos),
                       periods = c(pi, pi),
                       t = seq(0, 1, length.out = 51)){
  pt = sapply(periods, function(p, x) p * x, x = t)
  mat = matrix(NA, nrow = length(t), ncol = length(orthfuns))

  for(i in 1:length(orthfuns)){
    mat[,i] = sapply(pt[,i], orthfuns[[i]])
  }

  return(mat)
}


# s.t. t(M) %*% M / normal = Identity
#' Internal function to normalize randomly generated eigenscores and functions
#'
#' @param mat matrix to normalize
#' @param normal normalizing constant
#'
#' @return a normalized matrix
#' @noRd
#'
normalizer = function(mat, normal = NULL){
  if(is.null(normal)){normal = nrow(mat)}

  lmat = as.list(data.frame(mat))
  constants = sqrt(normal/sapply(lmat, function(x) sum(x^2)))

  return(matrix(c(constants[1]*mat[, 1],
                  constants[2]*mat[, 2]),
                nrow = nrow(mat),
                ncol = ncol(mat)))
}

#' Function to simulate missingness in simulated data
#'
#' @param n_obs number of observations to keep
#' @param dat original data
#' @param t time vector
#' @param seed random seed
#'
#' @return a data with artificial missing value
#' @export
#'
missing_data_simulator = function(n_obs,
                                  dat,
                                  t = seq(0, 1, length.out = 51),
                                  seed = 51){
  set.seed(seed)

  N = nrow(dat)
  l = length(t)
  c = ncol(dat)/l
  for(i in 1:N){
    s = sample(1:l, l - n_obs)
    S = s
    if(c > 1){
      for(j in 1:(c-1)){
        S = c(S, s + j*l)
      }
    }
    dat[i, ][S] = NA
  }

  # Add time row
  dat = rbind(dat, rep(t, c))

  # Transpose and change to data frame
  # dat = dat %>% t() %>% as.data.frame()

  return(dat)
}

#' Function to simulate irregular functional data
#'
#' @param seed random seed for replication
#' @param D variance matrix
#' @param N number of subjects
#' @param c number of components
#' @param t time vector; by default, a sequence from 0 to 1 with length 51
#' @param n_obs number of non-missing observations; must be shorter than t
#'
#'
#' @importFrom stats rnorm
#' @return entirely simulated data
#' @export
#'
simMFPCA = function(seed = 16,
                    D = matrix(c(3, 0, 0, 1), nrow = 2, ncol = 2, byrow = T),
                    N = 73,
                    c = 2,
                    t = seq(0, 1, length.out = 51),
                    n_obs = 10){
  # Simulate XiEst (i.e. eigen scores for each 73 individuals)
  # Ensure orthogonality
  set.seed(seed)
  U = matrix(rnorm(N*c), nrow = N, ncol = c) %>% normalizer()

  # Simulate Phi
  V1 = genOrthPair(orthfuns = c(sin, cos),
                   periods = c(pi, pi),
                   t = t) %>%
    normalizer(normal = 51/2)

  V2 = genOrthPair(orthfuns = c(sin, sin),
                   periods = c(pi, 2*pi),
                   t = t) %>%
    normalizer(normal = 51/2)
  V = rbind(V1, V2)

  # Check to ensure orthogonality
  tol = 1e-1
  check = 1/length(t)*t(V) %*% V
  if(check[1, 2] > tol | check[1, 2] < -tol){
    cat("The pair of functions were not orthogonal. \nPlease choose a different pair.\n")
  }

  X = U %*% D %*% t(V)
  full_data = X + matrix(rnorm(N*c*length(t), mean = 0, sd = 0.1),
                         nrow = N, ncol = c*length(t))

  miss_data = missing_data_simulator(n_obs, full_data)

  return(list("exact" = X,
              "obs.error" = full_data,
              "missing" = miss_data))
}
