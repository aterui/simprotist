#' Simulate metacommunity dynamics
#'
#' @param n_species Number of species in a metacommunity.
#' @param n_patch Number of patches in a metacommunity.
#' @param n_warmup Number of time-steps for warm-up. Default \code{200}.
#' @param n_burnin Number of time-steps for burn-in. Default \code{200}.
#' @param n_timestep Number of time-steps to be saved. Default \code{1000}.
#' @param type Model type. Either "bh" or "ricker".
#' @param propagule_interval Time interval for propagule introduction during warm-up. If \code{NULL}, a value of \code{ceiling(n_warmup / 10)} will be used.
#' @param propagule_seed Mean number of propagules
#' @param carrying_capacity Carrying capacity at each patch. Length should be one or equal to \code{n_patch}. Default \code{100}.
#' @param interaction_type Species interaction type. \code{"constant"} or \code{"random"}. \code{"constant"} assumes the unique interaction strength of \code{alpha} for all pairs of species. \code{"random"} draws random numbers from a uniform distribution with \code{min_alpha} and \code{max_alpha}.
#' @param alpha Species interaction strength. Enabled if \code{interaction_type = "constant"}. Default \code{0}.
#' @param min_alpha Minimum value of a uniform distribution that generates species interaction strength. Enabled if \code{interaction_type = "random"}. Default \code{NULL}.
#' @param max_alpha Maximum value of a uniform distribution that generates species interaction strength. Enabled if \code{interaction_type = "random"}. Default \code{NULL}.
#' @param r0 Maximum population growth rate at the optimal temperature
#' @param a Normalization constant
#' @param e_a Activation energy
#' @param e_d Deactivation energy
#' @param tmp_h Temperature at which r is half of r_peak.
#' @param scale_factor Scale factor.
#' @param xy_coord Each row should correspond to an individual patch, with x and y coordinates (columns). Must be give as a data frame. Defualt \code{NULL}.
#' @param distance_matrix Distance matrix indicating distance between habitat patches. If provided, the distance matrix will be used to generate dispersal matrix and to calculate distance decay of environmental correlations. Default \code{NULL}.
#' @param dispersal_matrix Dispersal matrix to be used to simulate dispersal process. Override distance_matrix. Default \code{NULL}.
#' @param boundary_condition Define boundary condition for dispersal. \code{retain} has not loss, \code{loss} induces net loss out of the network.
#' @param landscape_size Length of a landscape on a side. Enabled if \code{dispersal_matrix = NULL}.
#' @param mean_env Temperature at each patch. Length should be one or equal to \code{n_patch}.
#' @param sd_env Standard deviation of temporal temperature variation at each patch.
#' @param spatial_env_cor If \code{TRUE}, spatial autocorrelation in temporal environmental fluctuation is considered. Default \code{FALSE}.
#' @param phi Parameter describing the distance decay of spatial autocorrelation in temporal environmental fluctuation. Enabled if \code{spatial_env_cor = TRUE}.
#' @param p_dispersal Dispersal probability. Length should be one or equal to \code{n_species}.
#' @param theta Dispersal parameter describing dispersal capability of species.
#' @param dispersal_interval Interval for dispersal event.
#' @param plot If \code{TRUE}, five sample patches and species of \code{df_dynamics} are plotted.
#'
#' @return \code{df_dynamics} data frame containing simulated metacommunity dynamics.
#' @return \code{df_species} data frame containing species attributes.
#' @return \code{df_patch} data frame containing patch attributes.
#' @return \code{df_diversity} data frame containing diversity metrics. alpha and gamma diversities were averaged over the simulation time steps (\code{n_timestep}). beta diversity is calculated as gamma / alpha.
#'
#' @importFrom dplyr %>% filter
#' @importFrom ggplot2 ggplot vars labeller geom_line aes scale_color_viridis_c labs facet_grid label_both
#' @importFrom stats dist rbinom rgeom rnorm rpois runif
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom rlang .data
#'
#' @section Reference: see \href{https://github.com/aterui/mcbrnet}{github page} for instruction
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @examples
#' \dontrun{
#' # not run
#' mcsimp(n_patch = 5, n_species = 5)
#' }
#'
#' @export
#'
mcsimp <- function(n_species = 5,
                   n_patch = 5,
                   n_warmup = 200,
                   n_burnin = 200,
                   n_timestep = 1000,
                   type = "bh",
                   propagule_interval = NULL,
                   propagule_seed = 0.5,
                   carrying_capacity = 100,
                   interaction_type = "constant",
                   alpha = 0,
                   min_alpha = NULL,
                   max_alpha = NULL,
                   r0,
                   a,
                   e_a,
                   e_d,
                   tmp_h,
                   scale_factor,
                   xy_coord = NULL,
                   distance_matrix = NULL,
                   dispersal_matrix = NULL,
                   boundary_condition = "retain",
                   landscape_size = 10,
                   mean_env = 0,
                   sd_env = 0,
                   spatial_env_cor = FALSE,
                   phi = 1,
                   p_dispersal = 0.1,
                   theta = 1,
                   dispersal_interval = 1,
                   plot = FALSE
) {

  # parameter setup ---------------------------------------------------------

  # patch attributes
  ## internal function; see "fun_to_m.R"

  ## carrying capacity
  list_k <- fun_to_m(x = carrying_capacity,
                     n_species = n_species,
                     n_patch = n_patch,
                     param_attr = "patch")

  m_k <- list_k$m_x

  ## mean patch environment
  list_mu_z <- fun_to_m(x = mean_env,
                        n_species = n_species,
                        n_patch = n_patch,
                        param_attr = "patch")

  v_mu_z <- list_mu_z$v_x

  # species attributes
  ## internal function; see "fun_to_m.R"

  ## maximum growth at optimum
  if (any(r0 < 1)) stop("r0 must be >= one")

  list_r0 <- fun_to_m(x = r0,
                      n_species = n_species,
                      n_patch = n_patch,
                      param_attr = "species")

  v_r0 <- list_r0$v_x
  m_r0 <- list_r0$m_x

  ## temperature specific growth
  list_param <- list(a = a,
                     e_a = e_a,
                     e_d = e_d,
                     tmp_h = tmp_h)

  list_param_v <- lapply(list_param, function(x) {
    y <- fun_to_m(x = x,
                  n_species = n_species,
                  n_patch = n_patch,
                  param_attr = "species")
    return(y$v_x)
  })

  list_param_m <- lapply(list_param, function(x) {
    y <- fun_to_m(x = x,
                  n_species = n_species,
                  n_patch = n_patch,
                  param_attr = "species")
    return(y$m_x)
  })

  ## optimal temperature
  niche_optim <- (e_d * tmp_h) / (e_d + 8.6*1E-5 * tmp_h * log((e_d / e_a) - 1))

  # interaction matrix
  ## internal function; see "fun_int_mat.R"
  m_interaction <- fun_int_mat(n_species = n_species,
                               alpha = alpha,
                               min_alpha = min_alpha,
                               max_alpha = max_alpha,
                               interaction_type = interaction_type)

  # dispersal matrix
  ## internal function; see "fun_disp_mat.R"
  list_dispersal <- fun_disp_mat(n_patch = n_patch,
                                 landscape_size = landscape_size,
                                 theta = theta,
                                 xy_coord = xy_coord,
                                 distance_matrix = distance_matrix,
                                 dispersal_matrix = dispersal_matrix)

  m_distance <- list_dispersal$m_distance
  m_dispersal <- list_dispersal$m_dispersal
  df_xy_coord <- list_dispersal$df_xy_coord

  ## internal function; see "fun_to_v.R"
  list_p_dispersal <- fun_to_m(x = p_dispersal,
                               n_species = n_species,
                               n_patch = n_patch,
                               param_attr = "patch")

  m_p_dispersal <- list_p_dispersal$m_x

  # model
  fun_dyn <- fun_dyn_set(type)

  # dynamics ----------------------------------------------------------------

  ## object setup ####

  ## number of replicates
  n_sim <- n_warmup + n_burnin + n_timestep
  n_discard <- n_warmup + n_burnin

  ## environment
  var_env <- sd_env^2

  if (spatial_env_cor == TRUE) {

    if (is.null(m_distance)) stop("Provide distance matrix to model spatial environmental autocorrelation")
    m_sigma <- var_env * exp(-phi * m_distance)

  }else{

    m_sigma <- var_env * diag(x = 1, nrow = n_patch, ncol = n_patch)

  }

  m_z <- mvtnorm::rmvnorm(n = n_sim,
                          mean = v_mu_z,
                          sigma = m_sigma)

  ## community object
  colname <- c("timestep",
               "patch_id",
               "mean_env",
               "env",
               "carrying_capacity",
               "species",
               "niche_optim",
               "r_xt",
               "abundance")

  m_dynamics <- matrix(NA,
                       nrow = n_species * n_patch * n_timestep,
                       ncol = length(colname))

  colnames(m_dynamics) <- colname

  st_row <- seq(from = 1,
                to = nrow(m_dynamics),
                by = n_species * n_patch)

  ## propagule
  if (is.null(propagule_interval)) {

    propagule_interval <- ceiling(n_warmup / 10)

  }

  propagule <- seq(from = max(c(1, propagule_interval)),
                   to = max(c(1, n_warmup)),
                   by = max(c(1, propagule_interval)))

  ## dispersal
  if (dispersal_interval > n_sim) {

    psi <- rep(0, n_sim)

  } else {

    dispersal <- seq(from = dispersal_interval,
                     to = n_sim,
                     by = dispersal_interval)

    psi <- rep(0, n_sim)
    psi[dispersal] <- 1

  }

  ## initial values
  m_n <- matrix(rpois(n = n_species * n_patch,
                      lambda = propagule_seed),
                nrow = n_species,
                ncol = n_patch)

  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)

  ## temporal process ####
  for (n in seq_len(n_sim)) {

    if (n_warmup > 0) {

      if (n %in% propagule) {
        m_n <- m_n + matrix(rpois(n = n_species * n_patch,
                                  lambda = propagule_seed),
                            nrow = n_species,
                            ncol = n_patch)
      }

    }

    ## time-specific local environment
    m_z_xt <- matrix(rep(x = m_z[n, ],
                         each = n_species),
                     nrow = n_species,
                     ncol = n_patch)

    ## time-specific intrinsic growth rate
    m_r_xt <- fun_r_tmp(a = list_param_m$a,
                        e_a = list_param_m$e_a,
                        e_d = list_param_m$e_d,
                        k_b = 8.6*1E-5,
                        tmp = m_z_xt,
                        tmp_r = 298.15,
                        tmp_h = list_param_m$tmp_h,
                        scale_factor = scale_factor)

    ## internal community dynamics with competition
    m_n_hat <- fun_dyn(n = m_n,
                       r = m_r_xt,
                       r0 = m_r0,
                       k = m_k,
                       interaction = m_interaction)

    ## dispersal, internal function: see "fun_dispersal.R"
    m_n_bar <- fun_dispersal(x = m_n_hat,
                             m_p_dispersal = m_p_dispersal,
                             m_dispersal = m_dispersal,
                             boundary_condition = boundary_condition)

    m_n_prime <- psi[n] * m_n_bar + (1 - psi[n]) * m_n_hat

    ## demographic stochasticity
    m_n <- matrix(rpois(n = n_species * n_patch,
                        lambda = m_n_prime),
                  nrow = n_species,
                  ncol = n_patch)

    if (n > n_discard) {

      row_id <- seq(from = st_row[n - n_discard],
                    to = st_row[n - n_discard] + n_species * n_patch - 1,
                    by = 1)

      m_dynamics[row_id, ] <- cbind(# timestep
                                    I(n - n_discard),
                                    # patch_id
                                    rep(x = seq_len(n_patch), each = n_species),
                                    # mean_env
                                    rep(x = v_mu_z, each = n_species),
                                    # env
                                    c(m_z_xt),
                                    # carrying_capacity
                                    c(m_k),
                                    # species
                                    rep(x = seq_len(n_species), times = n_patch),
                                    # optimal niche
                                    rep(x = niche_optim, times = n_patch),
                                    # r_xt
                                    c(m_r_xt),
                                    # abundance
                                    c(m_n))

    }

    setTxtProgressBar(pb, n)

  }

  close(pb)


  # visualization -----------------------------------------------------------

  if (plot == TRUE) {

    sample_patch <- sample(seq_len(n_patch),
                           size = min(c(n_patch, 5)),
                           replace = FALSE)

    sample_species <- sample(seq_len(n_species),
                             size = min(c(n_species, 5)),
                             replace = FALSE)

    g <- dplyr::as_tibble(m_dynamics) %>%
      dplyr::filter(.data$patch_id %in% sample_patch,
                    .data$species %in% sample_species) %>%
      ggplot() +
      facet_grid(rows = vars(.data$species),
                 cols = vars(.data$patch_id),
                 labeller = labeller(.rows = label_both,
                                     .cols = label_both)) +
      geom_line(mapping = aes(x = .data$timestep,
                              y = .data$abundance,
                              color = abs(.data$niche_optim - .data$env))) +
      scale_color_viridis_c(alpha = 0.8) +
      labs(color = "Deviation \nfrom niche optimum")

    print(g)

  }


  # return ------------------------------------------------------------------

  # dynamics
  df_dyn <- dplyr::as_tibble(m_dynamics)

  # species attributes
  df_species <- df_dyn %>%
    dplyr::group_by(.data$species) %>%
    dplyr::summarise(mean_abundance = mean(.data$abundance)) %>%
    dplyr::right_join(dplyr::tibble(species = seq_len(n_species),
                                    r0 = v_r0,
                                    a = list_param_v$a,
                                    e_a = list_param_v$e_a,
                                    e_d = list_param_v$e_d,
                                    tmp_h = list_param_v$tmp_h,
                                    p_dispersal = rowMeans(m_p_dispersal)),
                      by = "species")

  # patch attributes
  df_patch <- df_dyn %>%
    dplyr::group_by(.data$patch_id) %>%
    dplyr::summarise(alpha_div = sum(.data$abundance > 0) / n_timestep) %>%
    dplyr::left_join(dplyr::tibble(patch_id = seq_len(n_patch),
                                   mean_env = mean_env,
                                   carrying_capacity = colMeans(m_k),
                                   connectivity = colSums(m_dispersal)),
                     by = "patch_id")

  # diversity metrics
  alpha_div <- sum(df_dyn$abundance > 0) / (n_timestep * n_patch)

  gamma_div <- df_dyn %>%
    dplyr::group_by(.data$timestep) %>%
    dplyr::summarise(gamma_t = dplyr::n_distinct(.data$species[.data$abundance > 0])) %>%
    dplyr::pull(.data$gamma_t) %>%
    mean()

  if (gamma_div == 0 | alpha_div == 0) {
    beta_div <- NA
  } else {
    beta_div <- gamma_div / alpha_div
  }

  # return
  return(list(df_dynamics = df_dyn,
              df_species = df_species,
              df_patch = df_patch,
              df_diversity = dplyr::tibble(alpha_div, beta_div, gamma_div),
              df_xy_coord = df_xy_coord,
              distance_matrix = m_distance,
              interaction_matrix = m_interaction))
}
