# Mallows kernel function ----

countTies <- function(x)
{
  # @param x Must be a matrix!
  # @return A vector representing number of tied pairs in each rows of "x".
  
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  } else if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  
  n1 <- apply(x, 1, function(u){
    tab <- table(u)
    sum(tab * (tab-1) / 2)
  })
  
  return(n1)
}

M_AllKendall = function (bandwidth, r, seqs = NULL) 
{
  if (is.vector(r)) {
    r <- matrix(r, nrow = 1)
  } else if (!is.matrix(r)) {
    r <- as.matrix(r)
  }
  kmat <- kendall_total(x = r, y = seqs)
  if (is.null(seqs)) {seqs <- r}
  stopifnot(ncol(r) == ncol(seqs))
  n0 <- choose(ncol(r), 2)
  v1 <- sqrt(n0 - countTies(r))
  v2 <- sqrt(n0 - countTies(seqs))
  
  stopifnot(nrow(kmat) == length(v1))
  stopifnot(ncol(kmat) == length(v2))
  kmat <- sweep(sweep(kmat, 1, v1, "*"), 2, v2, "*")
  dmat <- sweep(sweep(-2 * kmat, 1, v1 * v1, "+"), 2, v2 * 
                  v2, "+")
  return(exp(-bandwidth * (0.25 * dmat)))
}


