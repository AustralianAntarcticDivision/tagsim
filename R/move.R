#' Create movement matrix
#'
#' Create a transition matrix of movement probabilities
#' @param method movement method
#' @param n_rows number of rows in population matrix
#' @param n_cols number of columns in population matris
#' @param prob movement probability
#' @param matrix user supplied movement matrix
#' @export
move_matrix <- function(method, n_rows, n_cols, prob=NULL, matrix=NULL){
  ## add some checks
  ## for no movement
  if(method=="none"){
    obj <- diag(1, nrow=n_rows*n_cols, ncol=n_rows*n_cols)
  }else if(method=="complete"){
    prob = 1/(n_rows*n_cols)
    obj <- matrix(prob, nrow=n_rows*n_cols, ncol=n_rows*n_cols)
  }else if (method=="adjacent"){
    obj <- adjacted_move(n_rows=n_rows, n_cols=n_cols, prob=prob)
  }else if (method=="user"){
    obj <- matrix
  }
  ## add a class
  class(obj) <- "move"
  ## return the movement matrix
  obj
}

#' Create adjacent movement matrix
#'
#' Create a transition matrix of movement to adjacent cells
#'
#' @param n_rows number of rows in population matrix
#' @param n_cols number of columns in population matrix
#' @param prob probability of movement
#' @export
adjacted_move <- function(n_rows, n_cols, prob){
  ## check the probailiy of movement is <0.125
  if(prob>0.125) stop("Movement probability must be between 0 and 0.125")
  ## create an object to store the results
  obj <- matrix(0, nrow = n_rows*n_cols, ncol = n_rows*n_cols)
  ## need to loop over the rows and columns
  for(i in 1:n_rows){
    for(j in 1:n_cols){
      ## define array of ij's
      temp_ij <- data.frame(matrix(NA, nrow=8, ncol=2))
      names(temp_ij) <- c("rows", "cols")
      ## defined all possible positions
      temp_ij[1,] <- c(i-1, j-1)
      temp_ij[2,] <- c(i-1, j)
      temp_ij[3,] <- c(i-1, j+1)
      temp_ij[4,] <- c(i, j-1)
      temp_ij[5,] <- c(i, j+1)
      temp_ij[6,] <- c(i+1, j-1)
      temp_ij[7,] <- c(i+1, j)
      temp_ij[8,] <- c(i+1, j+1)
      ## remove those indices less than 1 or above n_rows
      temp_ij <- temp_ij[!(apply(temp_ij, 1, function(y) any(y == 0 | y > n_rows))),]
      ## loop over temp_ij
      for(k in 1:nrow(temp_ij)){
        ## assign to appropriate cells
        obj[i + (j-1)*n_rows, temp_ij[k,1] + ((temp_ij[k,2]-1)*n_rows)] <- prob
      }
      ## assign the rate for cell i,j
      obj[i + (j-1)*n_rows, i + (j-1)*n_rows] <- 1 - prob*nrow(temp_ij)
    }
  }
  ## return the matrix
  return(obj)
}

#' Move a fish population
#'
#' Implement movement of fish in a population array
#' @param d vector of starting population by region
#' @param move_mat movement transition matrix of
#' @export
move_N <- function(d, move_mat){
  ## checks
  if(!class(move_mat) %in% c("move")) stop("The movement matrix is not of class 'move'")
  if(any(d<0)) stop("vector of population d contains negative numbers")
  ## vector to store the final population
  obj <- rep(0, length(d))
  ## loop over the regions
  for(i in 1:length(d)){
    obj <- obj + d[i]*move_mat[i,]
  }
  return(obj)
}
