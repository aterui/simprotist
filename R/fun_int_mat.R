#' Internal function: to vector
#'
#' @param n_species number of species
#' @param alpha interaction strength
#' @param min_alpha minimum interaction strength
#' @param max_alpha maximum interaction strength
#' @param interaction_type type of interaction
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export
#'

fun_int_mat <- function(n_species,
                        alpha,
                        min_alpha,
                        max_alpha,
                        interaction_type) {

  if (interaction_type == "constant") {

    if (is.list(alpha)) {

      v_dim <- unlist(lapply(alpha, function(x) unique(dim(x))))
      if (any(v_dim != n_species)) stop("invalid dimension of alpha")
      if (any(unlist(alpha) < 0)) stop("invalid value of alpha")
      m_interaction <- alpha

    } else {

      if (any(alpha < 0)) stop("invalid value of alpha")

      if (length(alpha) == 1) {

        m_interaction <- matrix(alpha,
                                nrow = n_species,
                                ncol = n_species)

      }

      if (is.matrix(alpha)) {

        if (any(dim(alpha) != n_species)) stop("invalid dimension of alpha")
        m_interaction <- alpha

      }

    }

  } else {

    if (interaction_type != "random") stop("invalid interaction_type")
    if (is.null(min_alpha) | is.null(max_alpha)) stop("provide min_alpha and max_alpha")
    if (min_alpha < 0 | max_alpha < 0) stop("invalid values of min_alpha and/or max_alpha - values must be positive values")
    if (min_alpha > max_alpha) stop("max_alpha must exceed min_alpha")

    alpha <- runif(n = n_species * n_species,
                   min = min_alpha,
                   max = max_alpha)

    m_interaction <- matrix(alpha,
                            nrow = n_species,
                            ncol = n_species)

    diag(m_interaction) <- 1

  }

  return(m_interaction)

}
