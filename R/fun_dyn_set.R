#' Internal function: function for community dynamics
#' @param type Function type.
#' @export

fun_dyn_set <- function(type) {

  if (type == "bh") {

    fun_dyn <- function(n,
                        r,
                        r0 = NULL,
                        interaction,
                        k) {

      (n * r) / (1 + ((r0 - 1) / k) * (interaction %*% n))

    }

  }

  if (type == "ricker") {

    fun_dyn <- function(n,
                        r,
                        r0 = NULL,
                        interaction,
                        k) {

      n * r * exp(1 - (interaction %*% n) / k)

    }

  }

return(fun_dyn)

}
