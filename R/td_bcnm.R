
#' Temporal discounting binary choice model
#'
#' Compute a binary choice model for a single subject
#' @param data A data frame with columns \code{val_imm} and \code{val_del} for the values of the immediate and delayed rewards, \code{del} for the delay, and \code{imm_chosen} (Boolean) for whether the immediate reward was chosen. Other columns can also be present but will be ignored.
#' @param discount_function A string specifying the name of the discount functions to use, or an object of class \code{td_fn} (used for creating custom discount functions), or a list of objects of class \code{td_fn}.
#' @param choice_rule A string specifying whether the \code{'logistic'} (default), \code{'probit'}, or \code{'power'} choice rule should be used.
#' @param fixed_ends A Boolean (false by default) specifying whether the model should satisfy the desiderata that subjects should always prefer something over nothing (i.e., nonzero delayed reward over nothing) and the same reward sooner rather than later.
#' @param fit_err_rate A Boolean (false by default) specifying whether the model should include an error rate (parameterized by "eps"). See Eq. 5 here: https://doi.org/10.3758/s13428-015-0672-2
#' @param gamma_par_starts A vector of starting values to try for the "gamma" parameter (which controls the steepness of the choice rule) during optimization.
#' @param eps_par_starts A vector of starting values to try for the "eps" parameter (which controls the error rate) during optimization. Ignored if `fit_err_rate = FALSE`.
#' @param optim_args Additional arguments to pass to \code{optim()}. Default is \code{list(silent = T)}.
#' @param silent Boolean (true by default). The call to \code{optim()} occurs within a \code{try()} wrapper. The value of \code{silent} is passed along to \code{try()}.
#' @param na.action Action to take when data contains \code{NA} values. Default is \code{na.omit}.
#' @param ... Additional arguments to provide finer-grained control over the model configuration.
#' @return An object of class \code{td_bcnm} with components \code{data} (containing the data used for fitting), \code{config} (containing the internal configuration of the model, including the \code{discount_function}), and \code{optim} (the output of \code{optim()}).
#' @examples
#' \dontrun{
#' data("td_bc_single_ptpt")
#' mod <- td_bcnm(td_bc_single_ptpt, discount_function = "hyperbolic", fixed_ends = T)
#' # Custom discount function
#' custom_discount_function <- td_fn(
#'   name = 'custom',
#'   fn = function(data, p) (1 - p['b'])*exp(-p['k']*data$del) + p['b'],
#'   par_starts = list(k = c(0.001, 0.1), b = c(0.001, 0.1)),
#'   par_lims = list(k = c(0, Inf), b = c(0, 1)),
#'   ED50 = function(p) 'non-analytic'
#' )
#' mod <- td_bcnm(td_bc_single_ptpt, discount_function = custom_discount_function, fit_err_rate = T)
#' }
#' @export
td_bcnm <- function(
    data,
    discount_function = c('all',
                          'hyperbolic',
                          'exponential',
                          'inverse-q-exponential',
                          'nonlinear-time-hyperbolic',
                          'scaled-exponential',
                          'dual-systems-exponential',
                          'nonlinear-time-exponential',
                          'model-free',
                          'noise'),
    choice_rule = c('logistic', 'probit', 'power'),
    fixed_ends = F,
    fit_err_rate = F,
    fit_bias = F,
    # robust = F,
    # param_ranges = NULL,
    gamma_par_starts = c(0.01, 1, 100),
    eps_par_starts = c(0.01, 0.25),
    bias_par_starts = c(-0.1, 0, 0.1),
    silent = T,
    optim_args = list(),
    na.action = na.omit,
    ...) {
  
  # From a user's POV, it's easier to specify `choice_rule` and `fixed_ends`
  # Internally, it makes more sense to use `noise_dist`, `gamma_scale`, and `transform`
  # The user can control the latter through the `...` argument
  config <- list(...)
  if (length(config) == 0) {
    choice_rule <- match.arg(choice_rule)
    config <- list()
    # From `choice_rule` and `fixed_ends`, get `noise_dist`, `gamma_scale`, and `transform`
    if (choice_rule == 'logistic') {
      config$noise_dist <- 'logis'
      config$gamma_scale <- 'linear'
      if (fixed_ends) {
        config$transform <- 'noise_dist_quantile'
      } else {
        config$transform <- 'identity'
      }
    } else if (choice_rule == 'probit') {
      config$noise_dist <- 'norm'
      config$gamma_scale <- 'linear'
      if (fixed_ends) {
        config$transform <- 'noise_dist_quantile'
      } else {
        config$transform <- 'identity'
      }
    } else if (choice_rule == 'power') {
      config$noise_dist = 'logis'
      config$gamma_scale = 'none'
      if (fixed_ends) {
        config$transform = 'noise_dist_quantile'
      } else {
        config$transform = 'log'
      }
    }
    
  } else {
    
    req_args <- c('noise_dist', 'gamma_scale', 'transform')
    if (!setequal(req_args, names(config))) {
      stop(sprintf('If choice_rule is not specified, then the following must all be specified:\n%s', paste('- ', req_args, collapse = '\n')))
    
    }
  }
  
  # Set discount function(s)
  # This can't be as simple as match.arg because the user must also be able to provide a custom discount function
  if ('all' %in% discount_function) {
    # If "all" is used, replace discount_function with a vector of all the options
    discount_function <- eval(formals(td_fn)$predefined)
  }
  
  # Required data columns
  req_cols <- c('val_imm', 'val_del', 'del', 'imm_chosen')
  require_columns(data, req_cols)
  data <- data[req_cols]
  data <- na.action(data)
  
  # Ensure imm_chosen is logical
  data$imm_chosen <- as.logical(data$imm_chosen)
  
  # Attention checks
  attention_checks(data, warn = T)
  
  # All immediate chosen or all delayed chosen?
  invariance_checks(data, warn = T)
  
  # Valid discount function name
  validate_discount_function(discount_function)
  
  # Get arguments as a list
  args <- c(
    config,
    list(
      fit_err_rate = fit_err_rate,
      fit_bias = fit_bias
    )
  )
  
  # Get a list of candidate td_fn objects to test
  if (is.character(discount_function)) {
    # Vector of characters
    cand_fns <- list()
    for (fn_name in discount_function) {
      cand_fns <- c(cand_fns, list(td_fn(fn_name)))
    }
  } else {
    if (is(discount_function, 'td_fn')) {
      cand_fns <- list(discount_function) # List containing only one td_fn object
    } else {
      cand_fns <- discount_function # assumed to already be a list of td_fn objects, ensured by validation above
    }
  }
  
  # Run optimization for each discount function
  best_crit <- Inf
  best_mod <- list()
  cand_mod <- list(data = data)
  class(cand_mod) <- c('td_bcnm', 'td_um')
  for (cand_fn in cand_fns) {
    config <- args
    config$discount_function <- cand_fn
    
    cand_mod$config <- config
    
    # Get prob. model with the given settings but parameter values unspecified
    prob_mod_frame <- do.call(get_prob_mod_frame, config)
    # Get function to compute negative log likelihood
    robust <- F # Maybe in the future
    if (robust) {
      # nll_fn <- get_rob_fn(data, prob_mod_frame)
    } else {
      nll_fn <- get_nll_fn(data, prob_mod_frame)
    }
    
    # Get parameter starting values
    if (is.function(config$discount_function$par_starts)) {
      par_starts <- config$discount_function$par_starts(data)
    } else {
      par_starts <- config$discount_function$par_starts
    }
    # Add gamma start values
    par_starts <- c(
      par_starts,
      list(
        gamma = gamma_par_starts
      )
    )
    # Add epsilon start values, if fitting error rate
    if (config$fit_err_rate) {
      par_starts <- c(
        par_starts,
        list(
          eps = eps_par_starts
        )
      )
    }
    # Same for bias
    if (config$fit_bias) {
      par_starts <- c(
        par_starts,
        list(
          bias = bias_par_starts
        )
      )
    }
    
    # Get parameter bounds
    if (is.function(config$discount_function$par_lims)) {
      par_lims <- config$discount_function$par_lims(data)
    } else {
      par_lims <- config$discount_function$par_lims
    }
    # Add gamma limits
    par_lims <- c(
      par_lims,
      list(
        gamma = c(0, Inf)
      )
    )
    # Add epsilon limits, if fitting error rate
    if (config$fit_err_rate) {
      par_lims <- c(
        par_lims,
        list(
          eps = c(0, 0.5)
        )
      )
    }
    # Same for bias
    if (config$fit_bias) {
      par_lims <- c(
        par_lims,
        list(
          bias = c(-Inf, Inf)
        )
      )
    }
    # Run optimization
    optimized <- run_optimization(nll_fn,
                                  par_starts,
                                  par_lims,
                                  optim_args,
                                  silent)
    if (length(optimized) == 0) {
      stop('Optimization failed; optim() returned an error for every choice of initial parameter values')
    }
    
    cand_mod$optim <- optimized
    
    curr_crit <- BIC(cand_mod)
    
    if (curr_crit < best_crit) {
      best_crit <- curr_crit
      if ('par_chk' %in% names(config$discount_function)) {
        optimized$par <- config$discount_function$par_chk(optimized$par)
      }
      best_mod <- cand_mod
    }
  }
  best_mod$data <- data
  
  return(best_mod)
}


get_score_func_frame <- function(...) {
  # Get the function for computing scores without parameter values
  # specified, but with "structural" aspects (such as the discount
  # function) specified
  
  args <- list(...)
  
  # What is the discount function?
  discount_function <- args$discount_function$fn
  # Get the transform applied to v_i/v_d and f(t)
  transform <- get_transform(args)
  delta <- function(data, par) {
    transform(data$val_imm/data$val_del) - transform(laplace_smooth(discount_function(data, par)))
  }

  # Get the factor by which gamma is scaled
  gamma_scale <- switch(args$gamma_scale,
                        'none' = function(data, par) 1,
                        'linear' = function(data, par) data$val_del)
  
  if (args$fit_bias) {
    bias <- function(par) {par['bias']}
  } else {
    bias <- function(par) {0}
  }
  
  # Get the final frame
  frame <- function(data, par) {
    par['gamma'] * gamma_scale(data, par) * delta(data, par)} + bias(par)
  
  return(frame)
  
}

get_prob_func_frame <- function(...) {
  # Get probabilistic model without parameter values specified, but with
  # "structural" aspects (such as the CDF) specified
  
  args <- list(...)
  
  # What is the cdf?
  cdf <- get(sprintf('p%s', args$noise_dist))

  # Are we fitting error rate?
  if (args$fit_err_rate) {
    err_rate <- function(e, p) e + (1 - 2*e)*p
  } else {
    err_rate <- function(e, p) p
  }
  
  frame <- function(scores, par) err_rate(par['eps'], cdf( scores ))
  
  return(frame)
}

get_prob_mod_frame <- function(...) {
  # Get decision model without parameter values specified, but with
  # "structural" aspects (such as the discount function) specified
  score_func <- get_score_func_frame(...)
  prob_func <- get_prob_func_frame(...)
  frame <- function(data, par) {
    scores <- score_func(data, par)
    probs <- prob_func(scores, par)
    return(probs)
  }
}

get_nll_fn <- function(data, prob_mod_frame) {
  # Get negative log-likelihood function, given a set of data and a model
  # "frame" with structural aspects specified but parameters unspecified
  
  nll_fn <- function(par) {
    p <- laplace_smooth(prob_mod_frame(data, par))
    return(sum(-ll(p, data$imm_chosen)))
  }
  return(nll_fn)
}

# Robust stuff---may incorporate later
# huber <- function(t, c) {
#   idx <- t > c
#   t[idx] <- 2*sqrt(t[idx]*c) - c
#   # t[idx] <- c*log(c*t[idx]) - c*(2*log(c) - 1) # Same idea but with logarithm
#   return(t)
# }
# 
# get_rob_fn <- function(data, prob_mod_frame) {
#   # Get robust loss function, given a set of data and a model
#   # "frame" with structural aspects specified but parameters unspecified
#   
#   rob_fn <- function(par) {
#     p <- laplace_smooth(prob_mod_frame(data, par))
#     nll <- -ll(p, data$imm_chosen)
#     return(sum(huber(nll, 1)))
#   }
#   return(rob_fn)
# }
