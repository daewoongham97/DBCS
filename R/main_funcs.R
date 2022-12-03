# helper functions
cum_avg = function(x) {
  cumsum(x)/seq_along(x)
}

get_eta = function(t_star, alpha) {
  eta = sqrt((-lamW::lambertWm1(-alpha^2*exp(1)) - 1)/t_star)
  return(eta)
}

# main function to be exported
#' Constructs design based confidence sequence
#'
#' This function takes a dataset (input) with treatment and response and returns the
#' sequential confidence interval for the running mean of the average treatment effect.
#'
#' @param data A dataframe containing treatment and response (W, Y)
#' @param treatment Character string specifying the column name of data that contains treatment
#' Assume treatment to be numeric. E.g. (0,1) for binary.
#' @param response Character string specifying the column name of data that contains response
#' @param PS Optional vector of length nrow(data) that specifies the propensity scores.
#' If left blank, it is assumed to 0.5 for all units in the data.
#' @param alpha Numeric scalar for desired type 1 error control. Default of 0.05.
#' @param t_opt Integer scalar for choosing eta such that the width is smallest at t_opt.
#' Default at 10.
#' @param proxy_outcomes Optional vector of length nrow(data) where users can supply
#' own proxy outcomes to reduce variance. Default is NULL (no proxy outcomes).
#' @param nonasymp Boolean that specifies to return an asymptotic confidence sequence (FALSE)
#' or a non-asymptotic confidence sequence (TRUE). Default is FALSE.
#' @param M only supplied if nonasymp is TRUE. Hyperparameter corresponding to maximum response.
#' @param comparison Optional vector of length 2 used for treatments with more than 2 groups.
#' Example: comparison = c(0,1) obtains ATE of Y(1) - Y(0) and comparison = c(0, 2)
#' obtains ATE of Y(2) - Y(0). Comparison must match the same value of supplied treatment
#'
#' @return A list containing: \item{lower}{A vector containing the lower confidence sequence.}
#' \item{upper}{A vector containing the upper confidence sequence.}
#' \item{center}{A vector containing the center of the confidence sequence.}
#' \item{param}{A list specifying parameters related to the confidence sequence,
#' e.g., alpha, eta, etc.}
#'
#' @export
#' @references Ham, D., Bojinov, I., Lindon, M., Tingley, M. (2022)
#' Design-Based Confidence Sequences for Anytime-valid Causal Inference
DB_CS = function(data, treatment, response, PS = rep(0.5, nrow(data)), alpha = 0.05,
                 t_opt = 10, proxy_outcomes = NULL, nonasymp = FALSE, M = NULL, comparison = c(0, 1)) {
  # checks
  if (!(treatment %in% colnames(data))) stop("treatment should be in data")
  if (!(class(data[, treatment]) == "numeric")) stop("treatment should be numeric")
  if (!(response %in% colnames(data))) stop("response should be in data")
  if (nonasymp & is.null(M)) stop("M is required if nonasymp is TRUE")
  if (!(all(comparison %in%  unique(data[, treatment])))) stop("treatment does not match comparison")

  # setup
  y = data[, response]
  W = data[, treatment]
  t = nrow(data)
  in_idx = as.numeric(W %in% comparison)
  signW = W
  signW[signW == comparison[1]] = 0
  signW[signW == comparison[2]] = 1
  # main proposed asymptotic CS
  if (nonasymp == FALSE) {
    eta = get_eta(t_star = t_opt, alpha = alpha)
    if (is.null(proxy_outcomes)) {
      # without proxy outcomes
      tau_hat = y/PS*(-1)^(1-signW)*in_idx
      sig2_hat = y^2/PS^2*in_idx

    } else {
      # with proxy outcomes
      res_y = y - proxy_outcomes

      tau_hat = res_y/PS*(-1)^(1-signW)*in_idx
      sig2_hat = res_y^2/PS^2*in_idx

    }

    Sn = cumsum(sig2_hat)
    width = 1/(1:t)*sqrt(((Sn*eta^2 + 1)/eta^2)*log( (Sn*eta^2 + 1)/alpha^2))

    center = cum_avg(tau_hat)


  } else{
    p_min = min(PS)
    tau_hat = y/PS*(-1)^(1-signW)*in_idx
    sig2_hat = y^2/PS^2*in_idx
    Sn = cumsum(sig2_hat)
    m = M/p_min
    width = 1/(1:t)*(m*(m+1)*log(2/alpha) + Sn*((m+1)/m *log(1 + 1/m) - 1/m))
    center = cum_avg(tau_hat)

  }

  # Returning l
  lower = center - width
  upper = center + width
  param = list(nonasymp = nonasymp, alpha = alpha, PS = PS)
  l = list(lower = lower, upper = upper, center = center, param = param)

  return(l)

}
