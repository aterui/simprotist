#' Internal function: dispersal
#'
#' @param x Matrix of population sizes
#' @param m_p_dispersal Matrix of dispersal probability
#' @param m_dispersal Dispersal matrix
#' @param boundary_condition Define boundary condition for dispersal. \code{retain} has not loss, \code{loss} induces net loss out of the network.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

fun_dispersal <- function(x,
                          m_p_dispersal,
                          m_dispersal,
                          boundary_condition) {

  if (any(diag(m_dispersal) != 0)) stop("error in m_dispersal")

  # x is n_species (row) x n_patch (column) matrix
  if (boundary_condition == "retain") {
    # m_e_hat: expected number of emigrants from each habitat patch
    m_e_hat <- x * m_p_dispersal

    # m_e_sum: summed across patches
    v_e_sum <- rowSums(m_e_hat)

    # immigration potential for each patch = m_e_hat x m_dispersal (unit: individuals
    m_i_raw <- m_e_hat %*% m_dispersal

    # v_i_sum: summed across patches
    v_i_sum <- rowSums(m_i_raw)

    # insert 1 if v_i_sum = 0 (to avoid NaN for division)
    v_i_sum[v_i_sum == 0] <- 1

    # immigration probability = patch-specific potential / summed potential across patches
    m_i_prob <- m_i_raw / v_i_sum

    # expected immigrants: immigration prob. x total emigrants across patches
    m_i_hat <- m_i_prob * v_e_sum
  }

  if (boundary_condition == "loss") {

    m_e_hat <- x * m_p_dispersal
    m_i_hat <- x %*% m_dispersal

  }

  m_n_prime <- x + m_i_hat - m_e_hat

  return(m_n_prime)
}
