#' Internal function: function for community dynamics
#' @param n Population size
#' @param r Population growth rate
#' @param r0 Maximum population growth rate at the optimal temperature
#' @param interaction Interaction matrix
#' @param k Carrying capacity
#' @param z Power exponent in the Hassel model. Beverton-Holt (z = 1), under- (z < 1), or over-compensation (z > 1)
#' @export

hassel <- function(n, r, r0, interaction, k, z) {

      (n * r) / (1 + (r0 - 1) * ((interaction %*% n) / k)^z)

    }
