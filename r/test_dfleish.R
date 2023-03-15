dfleish_test <- function(x = 0, mean = 0, sd = 1, skew = 2) {
  stopifnot(requireNamespace("fGarch"))
  dsnorm_out <- dsnorm(x = x, mean = mean, sd = sd, xi = skew)
  dfleish_out <- dfleish(x = x, mean = mean, sd = sd, skew = skew)
  print(paste0("The results of dsnorm: ", dsnorm_out))
  print(paste0("The results of dfleish: ", dfleish_out))
}
