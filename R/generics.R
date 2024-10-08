
#' @export
print.td_bcnm <- function(x, ...) {
  cat(sprintf('\nTemporal discounting binary choice model\n\n'))
  cat(sprintf('Discount function: %s, with coefficients:\n\n', x$config$discount_function$name))
  print(coef(x))
  cat(sprintf('\nConfig:\n'))
  for (comp_name in c('noise_dist', 'gamma_scale', 'transform')) {
    cat(sprintf(' %s: %s\n', comp_name, x$config[[comp_name]]))
  }
  cat(sprintf('\nED50: %s\n', ED50(x)))
  cat(sprintf('AUC: %s\n', AUC(x, verbose = F)))
  cat(sprintf('BIC: %s\n', BIC(x)))
}

#' @export
print.td_bclm <- function(x, ...) {
  cat(sprintf('\nTemporal discounting binary choice linear model\n\n'))
  cat(sprintf('Discount function: %s from model %s, with coefficients:\n\n',
              x$config$discount_function$name,
              x$config$model))
  print(coef(x))
  NextMethod()
}

#' @export
print.td_ipm <- function(x, ...) {
  cat(sprintf('\nTemporal discounting indifference point model\n\n'))
  cat(sprintf('Discount function: %s, with coefficients:\n\n', x$config$discount_function$name))
  print(coef(x))
  cat(sprintf('\nED50: %s\n', ED50(x)))
  cat(sprintf('AUC: %s\n', AUC(x, verbose = F)))
}

#' @export
print.td_fn <- function(x, ...) {
  obj <- x
  cat(sprintf('\n"%s" temporal discounting function\n\n', obj$name))
  
  code <- deparse(body(obj$fn), width.cutoff = 500)
  code <- gsub('p\\["([^"]+)"\\]', '\\1', code)
  code <- gsub('data\\$', '', code)
  cat(sprintf('Indifference points:\n'))
  cat(paste(code, collapse = "\n"))
  
  cat(sprintf('\n\nParameter limits:\n'))
  for (par in names(obj$par_lims)) {
    cat(sprintf('%s < %s < %s\n', obj$par_lims[[par]][1], par, obj$par_lims[[par]][2]))
  }
}

#' Model Predictions
#'
#' Generate predictions from a temporal discounting binary choice model
#' @param object A temporal discounting binary choice model. See \code{td_bcnm}.
#' @param newdata Optionally, a data frame to use for prediction. If omitted, the data used to fit the model will be used for prediction.
#' @param type The type of prediction required. As in predict.glm, \code{"link"} (default) and \code{"response"} give predictions on the scales of the linear predictors and response variable, respectively. \code{"indiff"} gives predicted indifference points. In this case, \code{newdata} needs only a \code{del} column.
#' @param ... Additional arguments currently not used.
#' @return A vector of predictions
#' @examples
#' \dontrun{
#' data("td_bc_single_ptpt")
#' mod <- td_bcnm(td_bc_single_ptpt, discount_function = 'hyperbolic')
#' indiffs <- predict(mod, newdata = data.frame(del = 1:100), type = 'indiff')
#' }
#' @export
predict.td_bcnm <- function(object, newdata = NULL, type = c('link', 'response', 'indiff'), ...) {
  
  if (is.null(newdata)) {
    newdata <- object$data
  }
  
  type <- match.arg(type)
  
  if (type == 'link') {
    
    score_func <- do.call(get_score_func_frame, object$config)
    scores <- score_func(newdata, coef(object))
    return(scores)
    
  } else if (type == 'response') {
    
    prob_mod <- do.call(get_prob_mod_frame, object$config)
    probs <- prob_mod(newdata, coef(object))
    return(probs)
    
  } else if (type == 'indiff') {
    
    indiff_func <- object$config$discount_function$fn
    indiffs <- indiff_func(newdata, coef(object))
    names(indiffs) <- NULL
    return(indiffs)
    
  }
}

#' Model Predictions
#'
#' Generate predictions from a temporal discounting binary choice linear model
#' @param object A temporal discounting binary choice linear model. See \code{td_bclm}.
#' @param newdata Optionally, a data frame to use for prediction. If omitted, the data used to fit the model will be used for prediction.
#' @param type The type of prediction required. For \code{'indiff'} (default) gives predicted indifference points. In this case, \code{newdata} needs only a \code{del} column. For all other values (e.g. \code{"link"}, \code{"response"}), this function is just a wrapper to \code{predict.glm()}
#' @param ... Additional arguments passed to predict.glm if type != \code{'indiff'}
#' @return A vector of predictions
#' @examples
#' \dontrun{
#' data("td_bc_single_ptpt")
#' mod <- td_bclm(td_bc_single_ptpt, model = 'hyperbolic.1')
#' indiffs <- predict(mod, newdata = data.frame(del = 1:100), type = 'indiff')
#' }
#' @export
predict.td_bclm <- function(object, newdata = NULL, type = c('indiff', 'link', 'response', 'terms'), ...) {

  if (is.null(newdata)) {
    newdata <- object$data
  }
  
  type <- match.arg(type)
  if (type == 'indiff') {
    
    return(predict.td_bcnm(object, newdata = newdata, type = type))
    
  } else {
    
    newdata <- add_beta_terms(newdata, model = object$config$model)
    preds <- predict.glm(object, newdata = newdata, type = type, ...)
    return(preds)
    
  }
}

#' Model Predictions
#'
#' Generate predictions from a temporal discounting indifference point model
#' @param object A temporal discounting indifference point model. See \code{td_ipm}.
#' @param newdata A data frame to use for prediction. If omitted, the data used to fit the model will be used for prediction.
#' @param type Type of prediction, either \code{'indiff'} (indifference points) or \code{'response'} (whether the participants would is predicted to choose the immediate (1) or delayed reward (0))
#' @param ... Additional arguments currently not used.
#' @return A vector of predictions
#' @examples
#' \dontrun{
#' data("td_ip_simulated_ptpt")
#' mod <- td_ipm(td_ip_simulated_ptpt, discount_function = 'hyperbolic')
#' indiffs <- predict(mod, del = 1:100)
#' indiffs <- predict(mod, newdata = data.frame(del = 1:100))
#' }
#' @export
predict.td_ipm <- function(object, newdata = NULL, type = c('indiff', 'response'), ...) {
  
  if (is.null(newdata)) {
    if (length(list(...)) > 0) {
      newdata <- data.frame(...) # to enable predict(mod, del = 1:100) type syntax
    } else {
      newdata <- object$data
    }
  }
  
  indiff_func <- object$config$discount_function$fn
  indiffs <- indiff_func(newdata, coef(object))
  names(indiffs) <- NULL
  
  type <- match.arg(type)
  if (type == 'indiff') {
    out <- indiffs
  } else if (type == 'response') {
    require_columns(newdata, c('val_imm', 'val_del'))
    out <- as.numeric((newdata$val_imm / newdata$val_del) > indiffs)
  }
  
  return(out)
  
}

#' Get fitted values
#' 
#' Get fitted values of a temporal discounting binary choice model
#' @param object An object of class \code{td_bcnm}
#' @param ... Additional arguments currently not used.
#' @return A named vector of fitted values
#' @export
fitted.td_bcnm <- function(object, ...) {predict(object, type = 'response')}

#' Get fitted values
#' 
#' Get fitted values of a temporal discounting indifference point model
#' @param object An object of class \code{td_ipm}
#' @param ... Additional arguments currently not used.
#' @return A named vector of fitted values
#' @export
fitted.td_ipm <- function(object, ...) {predict(object)}

#' Extract model coefficients
#' 
#' Get coefficients of a temporal discounting binary choice model
#' @param object An object of class \code{td_bcnm}
#' @param ... Additional arguments currently not used.
#' @return A named vector of coefficients
#' @export
coef.td_bcnm <- function(object, ...) {object$optim$par}

#' Extract model coefficients
#' 
#' Get coefficients of a temporal discounting binary choice model
#' @param object An object of class \code{td_bcnm}
#' @param df_par Boolean specifying whether the coefficients returned should be the parameters of a discount function (versus the beta parameters from the regression)
#' @param ... Additional arguments currently not used.
#' @return A named vector of coefficients
#' @export
coef.td_bclm <- function(object, df_par = T, ...) {
  if (df_par) {
    # In terms of discount function parameters
    p <- object$coefficients
    d <- object$config$model
    if (d == 'hyperbolic.1') {
      cf <- c('k' = unname(p['.B2']/p['.B1']))
    } else if (d == 'hyperbolic.2') {
      cf <- c('k' = unname(exp(p['.B2']/p['.B1'])))
    } else if (d == 'exponential.1') {
      cf <- c('k' = unname(p['.B2']/p['.B1']))
    } else if (d == 'exponential.2') {
      cf <- c('k' = unname(exp(p['.B2']/p['.B1'])))
    } else if (d == 'scaled-exponential') {
      cf <- c('k' = unname(p['.B2']/p['.B1']),
              'w' = unname(exp(-p['.B3']/p['.B1'])))
    } else if (d == 'nonlinear-time-hyperbolic') {
      cf <- c('k' = unname(exp(p['.B3']/p['.B1'])),
              's' = unname(p['.B2']/p['.B1']))
    } else if (d == 'nonlinear-time-exponential') {
      cf <- c('k' = unname(exp(p['.B3']/p['.B1'])),
              's' = unname(p['.B2']/p['.B1']))
    } else if (d == 'itch') {
      cf <- object$coefficients
    } else if (d == 'naive') {
      cf <- object$coefficients
    }
  } else {
    cf <- object$coefficients
  }
  return(cf)
}

#' Extract model coefficients
#' 
#' Get coefficients of a temporal discounting indifference point model
#' @param object An object of class \code{td_ipm}
#' @param ... Additional arguments currently not used.
#' @return A named vector of coefficients
#' @export
coef.td_ipm <- function(object, ...) {object$optim$par}

#' Residuals from temporal discounting model
#'
#' Get residuals from a temporal discounting binary choice model
#' @param object A temporal discounting binary choice model. See \code{td_bcnm}.
#' @param type The type of residuals to be returned. See \code{residuals.glm}.
#' @param ... Additional arguments currently not used.
#' @return A vector of residuals
#' @export
residuals.td_bcnm <- function(object, type = c('deviance', 'pearson', 'response'), ...) {
  
  # args <- list(...)
  # type <- args$type
  # type <- match.arg(type, choices = c('deviance', 'pearson', 'response'))
  type <- match.arg(type)
  
  y <- object$data$imm_chosen
  yhat <- fitted(object)
  e <- y - yhat
  r <- switch (type,
               'deviance' = sign(e)*sqrt(-2*(y*log(yhat) + (1 - y)*log(1 - yhat))),
               'pearson' = e / sqrt(yhat[1]*(1 - yhat[1])),
               'response' = e
  )
  
  return(r)
}

#' Residuals from temporal discounting model
#'
#' Get residuals from a temporal discounting indifference point model
#' @param object A temporal discounting model. See \code{td_bcnm}.
#' @param type The type of residuals to be returned. See \code{residuals.nls}.
#' @param ... Additional arguments currently not used.
#' @return A vector of residuals
#' @export
residuals.td_ipm <- function(object, type = c('response', 'pearson'), ...) {
  
  # args <- list(...)
  # type <- match.arg(args$type, choices = c('response', 'pearson'))
  type <- match.arg(type)
  
  data <- object$data
  if ('indiff' %in% names(data)) {
    if (type == 'response') {
      y <- data$indiff
      yhat <- fitted(object)
      val <- y - yhat
    } else if (type == 'pearson') {
      # From residuals.nls
      val <- residuals(object, type = 'response')
      std <- sqrt(sum(val^2)/(length(val) - length(coef(object))))
      val <- val/std
    }
  } else {
    stop('Data was not fit directly on indifference points, so residuals cannot be computed.')
  }
  
  return(val)
}

#' Extract log-likelihood
#' 
#' Compute log-likelihood for a temporal discounting binary choice model.
#' @param mod An object of class \code{td_bcnm}
#' @export
logLik.td_bcnm <- function(mod) {
  p <- laplace_smooth(predict(mod, type = 'response'))
  x <- mod$data$imm_chosen
  val <- sum(ll(p, x))
  attr(val, "nobs") <- nrow(mod$data)
  attr(val, "df") <- length(coef(mod))
  class(val) <- "logLik"
  return(val)
}

#' Extract log-likelihood
#' 
#' Compute log-likelihood for a temporal discounting indifference point model.
#' @param mod An object of class \code{td_ipm}
#' @export
logLik.td_ipm <- function(mod) {
  
  # From logLik.nls
  res <- residuals(mod)
  N <- length(res)
  w <- rep_len(1, N) # Always unweighted
  ## Note the trick for zero weights
  zw <- w == 0
  val <-  -N * (log(2 * pi) + 1 - log(N) - sum(log(w + zw)) + log(sum(w*res^2)))/2
  ## the formula here corresponds to estimating sigma^2.
  attr(val, "df") <- 1L + length(coef(mod))
  attr(val, "nobs") <- attr(val, "nall") <- sum(!zw)
  class(val) <- "logLik"

  return(val)
}

#' Model deviance
#' 
#' Compute deviance for a temporal discounting binary choice model.
#' @param mod An object of class \code{td_bcnm}
#' @export
deviance.td_bcnm <- function(mod) return(-2*logLik.td_bcnm(mod))

#' Plot models
#'
#' Plot delay discounting models
#' @param x A delay discounting model. See \code{dd_prob_model} and \code{dd_det_model}
#' @param type Type of plot to generate
#' @param del Plots data for a particular delay
#' @param val_del Plots data for a particular delayed value
#' @param legend Logical: display a legend? Ignored if \code{type != 'summary'}
#' @param verbose Whether to print info about, e.g., setting del = ED50 when type == 'endpoints'
#' @param ... Additional arguments to \code{plot()}
#' @examples
#' \dontrun{
#' data("td_bc_single_ptpt")
#' mod <- td_bclm(td_bc_single_ptpt, model = 'hyperbolic.1')
#' plot(mod, type = 'summary')
#' plot(mod, type = 'endpoints')
#' }
#' @export
plot.td_um <- function(x, type = c('summary', 'endpoints', 'link'), legend = T, verbose = T, del = NULL, val_del = NULL, ...) {
  
  type <- match.arg(type)

  if (type == 'summary') {
    
    # Plot of binary choices or indifference points overlaid with discount function
    
    data <- x$data
    max_del <- max(data$del)
    min_del <- min(data$del)
    plotting_delays <- seq(min_del, max_del, length.out = 1000)
    if (is.null(val_del) & ('val_del' %in% names(x$data))) {
      val_del = mean(x$data$val_del)
      if (verbose) {
        cat(sprintf('Plotting indifference curve for val_del = %s (mean of val_del from data used to fit model). Override this behaviour by setting the `val_del` argument to plot() or set verbose = F to suppress this message.\n', val_del))
      }
    }
    pred_indiffs <- predict(x,
                            newdata = data.frame(del = plotting_delays,
                                                 val_del = val_del %def% NA),
                            type = 'indiff')
    
    # Set up axes
    plot(NA, NA,
         xlim = c(min_del, max_del), ylim = c(0, 1),
         xlab = 'Delay',
         ylab = 'val_imm / val_del',
         ...)
    
    # Plot indifference curve
    lines(pred_indiffs ~ plotting_delays)
    
    # Visualize stochasticity---goal for later. For now, don't know how to do this for td_bclm
    # if (x$config$gamma_scale != 'none') {
    #   if (verbose) {
    #     cat(sprintf('gamma parameter (steepness of curve) is scaled by val_del.\nThus, the curve will have different steepness for a different value of val_del.\nDefaulting to val_del = %s (mean of val_del from data used to fit model).\nUse the `val_del` argument to specify a custom value.\n\n', val_del))
    #   }
    # }
    # p_range <- args$p_range %def% c(0.4, 0.6)
    # lower <- invert_decision_function(x, prob = p_range[1], del = plotting_delays)
    # upper <- invert_decision_function(x, prob = p_range[2], del = plotting_delays)
    # lines(lower ~ plotting_delays, lty = 'dashed')
    # lines(upper ~ plotting_delays, lty = 'dashed')
    
    if ('indiff' %in% colnames(data)) {
      # Plot empirical indifference points
      points(indiff ~ del, data = data)
    } else {
      # Plot binary choices, immediate in red, delayed in blue
      data$rel_val <- data$val_imm / data$val_del
      points(rel_val ~ del, col = 'red',
             data = data[data$imm_chosen, ])
      points(rel_val ~ del, col = 'blue',
             data = data[!data$imm_chosen, ])
      if (legend) {
        legend("topright",
               inset = 0.02,
               title = 'Choices',
               legend = c("Imm.", "Del."),
               col = c("red", "blue"),
               pch = 1,
               box.lty = 0, # No border
               bg = rgb(1, 1, 1, 0.5)) # Background color with transparency
      }
    }
    
    title(x$config$discount_function$name)
    
  } else {
    if (is(x, 'td_ipm')) {
      
      stop('Only the "summary" plot type is applicable for models of type td_ipm.')
      
    } else {
      if (type == 'endpoints') {
        
        # Plot of psychometric curve
        
        if (is.null(val_del)) {
          val_del = mean(x$data$val_del)
          if (x$config$gamma_scale %def% 'none' != 'none') {
            if (verbose) {
              cat(sprintf('gamma parameter (steepness of psychometric curve curve) is scaled by val_del.\nThus, the curve will have different steepness for a different value of val_del.\nDefaulting to val_del = %s (mean of val_del from data used to fit model).\nUse the `val_del` argument to specify a custom value or use verbose = F to suppress this message.\n', val_del))
            }
          }
        }
        
        if (is.null(del)) {
          del <- ED50(x, val_del = val_del)
          if (del == 'none') {
            del <- mean(c(min(x$data$del), max(c(x$data$del))))
            if (verbose) {
              cat(sprintf('ED50 is undefined. Therefore setting del=%s (halfway between min. delay and max. delay from fitting data).\nThis can be specified manually with the `del` argument.\n\n', del))
            }
          } else {
            if (verbose) {
              cat(sprintf('Setting del = %s (ED50) to center the curve.\nThis can be changed using the `del` argument.\n\n', del))
            }
          }
        }
        
        plotting_rs <- seq(0, 1, length.out = 1000)
        newdata <- data.frame(
          del = del,
          val_del = val_del,
          val_imm = plotting_rs*val_del
        )
        p <- predict(x, newdata = newdata, type = 'response')
        plot(p ~ plotting_rs, type = 'l',
             ylim = c(0, 1),
             xlab = 'val_imm/val_del',
             ylab = 'Prob. Imm',
             ...)
        
        # If applicable, plot the choices at the given delay
        if (del %in% x$data$del) {
          sdf <- x$data[x$data$del == del, ]
          sdf$R <-sdf$val_imm / sdf$val_del
          points(imm_chosen ~ R, data = sdf)
        }
        
        title(sprintf('del = %s, val_del = %s', del, val_del))
        
      } else if (type == 'link') {
        
        # Plot of probabilities and data against linear predictors
        
        # Get score range
        if (is(x, 'td_bcnm')) {
          score_func <- do.call(get_score_func_frame, x$config)
          scores <- score_func(x$data, coef(x))
        } else if (is(x, 'td_bclm')) {
          scores <- x$linear.predictors
        }
        lim <- max(abs(min(scores)), abs(max(scores)))
        # Plot choices
        plot(x$data$imm_chosen ~ scores,
             ylim = c(0, 1),
             xlim = c(-lim, lim),
             ylab = 'imm_chosen',
             xlab = 'Linear predictor',
             ...)
        # Plot probabilities
        plotting_scores <- seq(-lim, lim, length.out = 1000)
        if (is(x, 'td_bcnm')) {
          prob_func <- do.call(get_prob_func_frame, x$config)
          p <- prob_func(plotting_scores, coef(x))
        } else if (is(x, 'td_bclm')) {
          p <- x$family$linkinv(plotting_scores)
        }
        lines(p ~ plotting_scores)
        
      }
    }
  }
}

