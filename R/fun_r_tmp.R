#' Internal function: function for temperature-dependent population growth
#' @param a normalization constant.
#' @param e_a activation energy
#' @param e_d deactivation energy
#' @param k_b Boltzmann's constant
#' @param tmp temperature
#' @param tmp_r reference temperature
#' @param tmp_h half-r temperature
#' @export

fun_r_tmp <- function(a,
                      e_a,
                      e_d,
                      k_b = 8.6*1E-5,
                      tmp,
                      tmp_r = 298.15,
                      tmp_h,
                      scale_factor = 0) {

  log_r <- a +
    (e_a / k_b) * ((1 / tmp_r) - (1 / tmp)) -
    log(1 + exp((e_d / k_b) * ((1 / tmp_h) - (1 / tmp))))

  r <- exp(log_r) - scale_factor
  return(r)
}
