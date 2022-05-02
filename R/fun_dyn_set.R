#' Internal function: function for community dynamics
#' @param type Function type.
#' @export

fun_dyn_set <- function(type) {

  if (type == "bh") {

    fun_dyn <- function(n,
                        r,
                        r0 = NULL,
                        interaction,
                        k,
                        fix_k = TRUE) {

      if (fix_k == TRUE) {

        r_max <- r

      } else {

        r_max <- r0

      }

      (n * r) / (1 + ((r_max - 1) / k) * (interaction %*% n))

    }

  }

  if (type == "ricker") {

    fun_dyn <- function(n,
                        r,
                        r0 = NULL,
                        interaction,
                        k,
                        fix_k = TRUE) {

      if (fix_k == TRUE) {

        r_max <- r

      } else {

        r_max <- r0

      }

      n * r_max * exp(1 - (interaction %*% n) / k)

    }

  }

return(fun_dyn)

}
