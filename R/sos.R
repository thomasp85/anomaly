#' Calculates the perplexity and affinity for a single point at a given precision
#' 
#' This function is equivalent to the Python \code{get_perplexity()}. It takes
#' the distances from one point to all others in the dataset and calculates the
#' affinites for that point at a specific precision.
#' 
#' @param d A vector of distances from one point to all the rest
#' 
#' @param beta The precision of the gaussian distribution used
#' 
#' @return A list with the following entries:
#' \describe{
#'  \item{A}{A vector of affinities to the points given in d}
#'  \item{H}{The perplexity of the point}
#' }
#' 
#' @noRd
#' 
getPerplexity <- function(d, beta) {
    A <- exp(-d*beta)
    H <- log(sum(A)) + beta*sum(d*A)/sum(A)
    list(A=A, H=H)
}
#' Calculates the affinites for a given distance matrix
#' 
#' This function is equivalent to the Python \code{d2a()}. It takes a distance
#' matrix and a desired perplexity (optionally a tolerance) and calculates the
#' affinities for all points by individually adjusting the gaussian kernel to 
#' satisfy the desired perplexity for each point.
#' 
#' @param dist A distance matrix; usually as returned by \code{\link[stats]{dist}}
#' 
#' @param perplexity The desired perplexity
#' 
#' @param tolerance The acceptable difference between the obtained and desired
#' perplexity
#' 
#' @return A matrix matching the dimensions of dist with the affinities for each
#' point. Contrary to dist is it not symmetrical. It should be read row-wise in
#' the sense that at row i you will have sample i's affinity to all other 
#' points. Conversily at column i you will have al other points affinity toward
#' sample i.
#' 
#' @noRd
#' 
affinity <- function(dist, perplexity, tolerance=1e-5) {
    maxIter <- 50
    dist <- as.matrix(dist)
    n <- nrow(dist)
    ans <- matrix(0, ncol=n, nrow=n)
    beta <- rep(1, n)
    logU <- log(perplexity)
    
    for(i in 1:n) {
        Di <- dist[i, -i]
        betamin = -Inf
        betamax = Inf
        perp <- getPerplexity(Di, beta[i])
        Hdiff <- perp$H - logU
        tries <- 0
        while(is.nan(Hdiff) || (abs(Hdiff) > tolerance && tries < maxIter)) {
            if(is.nan(Hdiff)) {
                beta[i] <- beta[i]/10
            } else if(Hdiff > 0) {
                betamin <- beta[i]
                if(betamax == Inf || betamax == -Inf) {
                    beta[i] <- beta[i] * 2
                } else {
                    beta[i] <- (beta[i]+betamax) / 2
                }
            } else {
                betamax <- beta[i]
                if(betamin == Inf || betamin == -Inf) {
                    beta[i] <- beta[i] / 2
                } else {
                    beta[i] <- (beta[i]+betamin) / 2
                }
            }
            perp <- getPerplexity(Di, beta[i])
            Hdiff <- perp$H - logU
            tries <- tries + 1
        }
        ans[i, -i] <- perp$A
    }
    ans
}
#' Calculate the binding probabilities given an affinity matrix
#' 
#' This function is equivalent to the Python \code{a2b()}. It calculates the
#' probability that a point binds to another point in a stochastic neighbor
#' graph. The binding probability is proportional to the affinity and is in fact
#' just the affinity normalised to the sum of all affinities for a point.
#' 
#' @param affinity An affinity matrix as returned by \code{\link{affinity}}
#' 
#' @return A matrix matching the dimensions of affinity, with the binding 
#' probabilities for each directed edge
#' 
#' @noRd
#' 
bindingProb <- function(affinity) {
    affinity/apply(affinity, 1, sum)
}
#' Calculate the outlier probabilities given a binding probabilty matrix
#' 
#' This function is equivalent to the Python \code{b2o()}. It calculates for 
#' each point the probabilty that its in-degree in a stochastic neighbor graph 
#' is 0.
#' 
#' @param bindingP A matrix with binding probabilities as returned by 
#' \code{\link{bindingProb}}
#' 
#' @return A vector of length equal to the number of rows in bindingP, with the
#' outlier probability for each point.
#' 
#' @noRd
#' 
outlierProb <- function(bindingP) {
    apply(1-bindingP, 2, prod)
}
#' Perform stochastic outlier selection on a distance matrix
#' 
#' This function calculates for each point in a distance matrix its probabilty
#' of being an outlier, given a certain perplexity. The perplexity parameter can
#' be seen as the size of the neighborhood taken into account when assessing the
#' outlierness of a given point, though contrary to k in KNNDD and LOF it does 
#' not have to be an integer but can take any positive value (should be lower 
#' than the number of points though). By setting a probability threshold the 
#' output can be converted to a hard outlier classification.
#' 
#' @param dist A distance matrix or a dist object as returned by \code{\link[stats]{dist}}
#' 
#' @param perplexity The size of the neighborhood to be considered
#' 
#' @return A list with the following entries:
#' \describe{
#'  \item{outlierProbabilty}{A vector of outlier probabilities}
#'  \item{bindingProbability}{A matrix with the binding probabilities between all point}
#' }
#' 
#' @examples
#' irisDist <- dist(iris[, 1:4])
#' irisSos <- sos(irisDist, 65)
#' irisMDS <- cmdscale(irisDist)
#' plot(irisMDS, pch=19, col=ifelse(irisSos$outlierProbabilty > 0.55, 'red', 'black'))
#' 
#' @export
#' 
sos <- function(dist, perplexity) {
    A <- affinity(dist, perplexity)
    B <- bindingProb(A)
    O <- outlierProb(B)
    list(outlierProbabilty=O, bindingProbability=B)
}