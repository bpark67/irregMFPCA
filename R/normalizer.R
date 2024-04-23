#' Function to normalize eigenscores and eigenfunctions
#'
#' @param mat matrix of eigenscores or eigenfunctions
#' @param normal normalizing constant; by default set to the number of rows
#'
#' @return
#' @export
#'
#' @examples
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
