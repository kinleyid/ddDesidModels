% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonlinear-model-fitting.R
\name{dd_prob_model}
\alias{dd_prob_model}
\title{Probabilistic delay discounting model}
\usage{
dd_prob_model(
  data,
  discount_function = "all",
  absval = "none",
  dplus = T,
  param_ranges = NULL,
  silent = T
)
}
\arguments{
\item{data}{A data frame with columns `val_imm` and `val_del` for the values of the immediate and delayed rewards, `del` for the delay, and `imm_chosen` (Boolean) for whether the immediate reward was chosen}

\item{discount_function}{A vector of strings specifying the name of the discount functions to use. Options are `'hyperbolic'`, `'exponential'`, `'inverse-q-exponential'`, `'nonlinear-time-hyperbolic'`, `'scaled-exponential'`, `'dual-systems-exponential'`, and `'nonlinear-time-exponential'`. Default is `'all'`, meaning every discount function is tested and the one with the best AIC is selected.}

\item{absval}{A string specifying how the absolute value of the delayed reward should be accounted for. Defaults to `'none'`. Other options are `'identity'` (linear scaling) and `'varsigma'` (flexible nonlinear scaling)}

\item{dplus}{A Boolean specifying whether the model should satisfy the desiderata that subjects should always prefer something over nothing (i.e., nonzero delayed reward over nothing) and the same reward sooner rather than later}

\item{param_ranges}{A list containing the starting values to try for each parameter. Defaults to `c(-5, 0, 5)` for most parameters}

\item{silent}{A Boolean specifying whether the call to `optim` (which occurs in a `try` block) should be silent on error}
}
\value{
A list from `optim` with additional components specifying the AIC, the discount function, and the probabilistic model
}
\description{
Compute a probabilistic model for a single subject's delay discounting
}
\note{
The `par` component of the output list is for internal use. For statistical analyses, use the `untransformed_parameters`
}
\examples{
df <- data.frame(val_imm = seq(1, 99, length.out = 10), val_del = 100, del = rep(exp(1:10), each=10))
logistic <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x / (1 - x))
prob <- logistic(logit(df$val_imm / df$val_del) - logit(1 / (1 + 0.001*df$del)))
df$imm_chosen <- runif(nrow(df)) < prob
mod <- dd_prob_model(df)
print(mod$discount_function_name)
}
