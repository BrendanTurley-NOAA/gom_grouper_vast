### https://github.com/James-Thorson-NOAA/VAST/wiki/Visualize-covariate-response

load("~/Desktop/professional/projects/Postdoc_FL/data/grouper/2022-05-17_vmod14_results.RData")

#####################
# Effects package
#####################

library(effects)  # Used to visualize covariate effects

# Must add data-frames to global environment (hope to fix in future)
covariate_data_full = fit$effects$covariate_data_full
catchability_data_full = fit$effects$catchability_data_full

# Plot 1st linear predictor, but could use `transformation` to apply link function
pred = Effect.fit_model( fit,
                         focal.predictors = c("bot_do"),
                         which_formula = "X1",
                         xlevels = 100)
                         # transformation = list(link=identity, inverse=identity) )
plot(pred)

#####################
# pdp package
#####################

library(pdp)

units_options(allow_mixed = TRUE)

# Make function to interface with pdp
pred.fun = function( object, newdata ){
  predict( x=object,
           Lat_i = object$data_frame$Lat_i,
           Lon_i = object$data_frame$Lon_i,
           t_i = object$data_frame$t_i,
           a_i = object$data_frame$a_i,
           what = "P1_iz",
           new_covariate_data = newdata,
           do_checks = FALSE )
}

# Run partial
Partial = partial( object = fit,
                   pred.var = "bot_do",
                   pred.fun = pred.fun,
                   train = fit$covariate_data )

# Make plot using ggplot2
library(ggplot2)
autoplot(Partial)
