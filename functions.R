# Yichao Wang (s2056237)

# Place your function definitions that may be needed in the report.Rmd, including function documentation.
# You can also include any needed library() calls here
library(StatCompLab)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(knitr)
library(dplyr)

###################

#' P-value and confidence interval
#' Construct the 95% confidence interval for each station
#'
#' The result will return a table includes the Name, p-value, 
#' Monte Carlo standard deviation, and the CI for each station.
#'
#' @param data the desired data to be put into model
#' @param alpha the value of alpha, here we take 0.05 to construct a 95% CI
#' @param N the number of permutations
#' 

p_value_CI <- function(data, alpha=0.05, N=10000) {
  # Initialize the count to be zero for each test statistic
  count <- rep(0, length = length(T_stats))
  for (n in 1:N) {
    # Generate a random permutation of the indices
    permutation <- sample(data$Summer, length(data$Summer), replace = FALSE)
    # Apply the permutation to the data and group by weather station
    ghcnd_permutation <- data %>%
      mutate("Summer" = permutation) %>%
      group_by(Name)
    # Compute the test statistic for the permuted data
    ghcnd_permutated <- ghcnd_permutation %>%
      summarise(Value = abs(mean(Value[Summer]) - mean(Value[!Summer])), .groups = "drop")
    # add one to the count if the test statistic is less than or equal to the observed value
    count <- count + (ghcnd_permutated$Value >= T_stats)
  }
  
  # Compute the p-value and Monte Carlo standard deviation
  p_value <- count / N
  monte_carlo_sd <- sqrt(p_value * (1 - p_value) / N)
  
  # Compute confidence interval for p-value
  lower_CI <- numeric(length(p_value))
  upper_CI <- numeric(length(p_value))
  for (i in 1:length(p_value)) {
    if (p_value[i] > 0) { # When x > 0 (p_value>0)
      lower_CI[i] <- p_value[i] - qnorm(1 - alpha/2) * monte_carlo_sd[i]
      upper_CI[i] <- p_value[i] + qnorm(1 - alpha/2) * monte_carlo_sd[i]
    } else { # When x = 0, CI=(0, 1-0.0025^(1/N))
      lower_CI[i] <- 0
      upper_CI[i] <- 1 - 0.025^(1/N)
    }
  }
  
  # Store the result in a dataframe and show them as a table
  result_table <- data.frame(Name = unique(winter_avg$Name), 
                             p_value = p_value,
                             monte_carlo_sd = monte_carlo_sd,
                             lower_CI = lower_CI,
                             upper_CI = upper_CI)
  return(kable(result_table, align = "lcccc"))
}




#' Estimate models for the square root of the 
#' monthly averaged precipitation values in Scotland
#'
#' The result will return the fitted model for five different 
#' model M0, M1, M2, M3, M4 as different value of k entered
#'
#' @param data the desired data to be put into model
#' @param k the model M_k(k from 0 to 4)
#' 

fit_model <- function(data, k) {
  # Define the formula for the model
  if(k == 0) { # When k = 0 
    model <- Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear
  } else{ # When k = 1,2,3,4
    model <- formula(paste0("Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear + ",
                            paste0("cos(2*pi*", seq(1, k), "*DecYear) + sin(2*pi*", seq(1, k), "*DecYear)", 
                                   collapse = " + ")))
  }
  # Fit the linear model
  fit <- lm(model, data = data)
  
  # Return the model object
  return(fit)
}




#' Stratified cross validation function
#' 
#' The result will return a 96*5 dataframe shows the average SE and DS scores for each 
#' month by station after the corss validation process, including the Name, ID, Month,
#' and its corresponding standard error scores(SE) and Dawid-Sebastiani scores(DS). 
#'
#' since there are 8 stations and 12 months, we got a total of 12*8=96 rows.
#' 
#' @param data the desired data to do the corss validation
#' @param k the model M_k(k from 0 to 4)

stratified_cv <- function(data, k) {
  pred2 <- data %>% mutate(mean = NA_real_, sd = NA_real_)
  for (id in unique(data$ID)) {
    if(k == 0) { # When k = 0 
      model <- Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear
    } else{ # When k = 1,2,3,4
      model <- formula(paste0("Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear + ",
                              paste0("cos(2*pi*", seq(1, k), "*DecYear) + sin(2*pi*", seq(1, k), "*DecYear)", 
                                     collapse = " + ")))
    }
    # fit and predict the model for the data
    fit <- lm(model, data %>% filter(ID != id))
    pred <- predict(fit, newdata = data %>% filter(ID == id), se.fit = TRUE)
    # evaluate the mean and sd
    pred2$mean[(pred2$ID == id)] <- pred$fit
    pred2$sd[(pred2$ID == id)] <- sqrt((pred$se.fit)^2 + (pred$residual.scale)^2)
    # mutate the table
    pred2 <- pred2 %>%
      group_by(ID, Name, Year, Month) %>%
      mutate(se = proper_score("se", Value_sqrt_avg, mean = mean),
             ds = proper_score("ds", Value_sqrt_avg, mean = mean, sd = sd))
  }
  # Calculate the average score for each month
  Average_score <- pred2 %>%
    group_by(ID, Name, Month) %>%
    summarise(se = mean(se), ds = mean(ds), .groups = "drop")
  return(Average_score)
}

