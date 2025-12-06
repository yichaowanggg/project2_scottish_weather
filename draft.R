library(ggplot2)

data(ghcnd_stations, package = "StatCompLab")
data(ghcnd_values, package = "StatCompLab")

ghcnd <- left_join(ghcnd_values, ghcnd_stations, by = "ID")

# Filter data for selected one year(2000)
ghcnd_year <- filter(ghcnd, Year == 2015)

# Create plot
ggplot(ghcnd_year, aes(x = DecYear, y = Value, color = Element)) +
  geom_point() +
  facet_wrap(~ Name) +
  labs(x = "Month",
       y = "Temperature (°C) / Precipitation (mm)",
       title = paste0("Temperature and Precipitation for ", 2000)) +
  scale_x_continuous(breaks = 1:12,
                     labels = c(1:12))




# Compute separate p-values for each weather station 
# and their respective Monte Carlo standard deviations. 
# Construct a function p_value_CI to construct the interval.
library(dplyr)

# Define the function to compute the p-value and confidence interval
p_value_CI <- function(T, alpha = 0.05, n = 1000) {
  # Compute the observed test statistic
  obs_T <- abs(mean(T[T > 0]) - mean(T[T <= 0]))
  
  # Generate null distribution using Monte Carlo simulation
  T_sim <- replicate(n, {
    # Shuffle the data randomly
    T_perm <- sample(T, replace = TRUE)
    
    # Compute the test statistic for the permuted data
    abs(mean(T_perm[T_perm > 0]) - mean(T_perm[T_perm <= 0]))
  })
  
  # Compute the p-value
  p_value <- mean(T_sim >= obs_T) + mean(T_sim <= -obs_T)
  
  # Compute the Monte Carlo standard deviation
  mc_sd <- sd(T_sim)
  
  # Compute the confidence interval
  ci <- c(obs_T - qnorm(1 - alpha/2) * mc_sd, obs_T + qnorm(1 - alpha/2) * mc_sd)
  
  # Return the results
  list(p_value = p_value, mc_sd = mc_sd, ci = ci)
}

# Compute winter and summer averages for each station
# Winter
ghcnd_winter <- ghcnd_selected %>%
  filter(Summer == FALSE)

winter_avg <- ghcnd_winter %>%
  group_by(Name, ID) %>%
  summarise(winter_average = mean(Value), .groups = "drop")

# Summer
ghcnd_summer <- ghcnd_selected %>%
  filter(Summer == TRUE)

summer_avg <- ghcnd_summer %>%
  group_by(Name, ID) %>%
  summarise(summer_average = mean(Value), .groups = "drop")

# Calculate the test statistic
T_stats <- winter_avg$winter_average - summer_avg$summer_average

# Compute p-values and Monte Carlo standard deviations for each station
p_values <- sapply(T_stats, function(T) {
  p_value_CI(T)$p_value
})
mc_sds <- sapply(T_stats, function(T) {
  p_value_CI(T)$mc_sd
})

# Combine the results into a table
result_table <- cbind(winter_avg, summer_avg, T_stats, p_values, mc_sds)

# Display the table
knitr::kable(result_table, align = "lcccccc")




##########################################
##########################################

monte_carlo_permutation_test <- function(x, y, n_permutations = 1000) {
  # Compute the observed test statistic
  obs_stat <- abs(mean(x) - mean(y))
  # Combine the data into a single vector
  xy <- c(x, y)
  # Initialize the vector to store the permutation test statistics
  perm_stats <- numeric(n_permutations)
  # Permute the data and compute the test statistic for each permutation
  for (i in 1:n_permutations) {
    # Generate a random permutation of the indices
    perm <- sample(length(xy))
    # Compute the test statistic for the permuted data
    perm_x <- xy[perm[1:length(x)]]
    perm_y <- xy[perm[(length(x) + 1):length(xy)]]
    perm_stats[i] <- abs(mean(perm_x) - mean(perm_y))
  }
  # Compute the p-value
  p_value <- sum(perm_stats >= obs_stat) / n_permutations
  # Compute the standard deviation of the test statistic
  sd_stat <- sd(perm_stats)
  # Return the results as a named list
  return(list(p_value = p_value, sd = sd_stat))
}



p_value_CI <- function(data, summer_months, alpha=0.05, num_perm=1000) {
  # Compute winter and summer averages for each station
  winter_data <- data %>%
    filter(Month %in% setdiff(1:12, summer_months))
  summer_data <- data %>%
    filter(Month %in% summer_months)
  
  winter_avg <- winter_data %>%
    group_by(ID, Name) %>%
    summarise(winter_average = mean(Value), .groups = "drop")
  
  summer_avg <- summer_data %>%
    group_by(ID, Name) %>%
    summarise(summer_average = mean(Value), .groups = "drop")
  
  # Calculate the test statistic for each station
  T_stats <- abs(winter_avg$winter_average - summer_avg$summer_average)
  
  # Generate permutations and compute p-values
  p_values <- numeric(length(T_stats))
  for (i in seq_along(T_stats)) {
    permuted_data <- data %>%
      filter(ID == winter_avg$ID[i]) %>%
      mutate(Summer = sample(Summer, replace=FALSE)) 
    permuted_winter <- permuted_data %>%
      filter(Month %in% setdiff(1:12, summer_months)) %>%
      summarise(winter_average = mean(Value))
    permuted_summer <- permuted_data %>%
      filter(Month %in% summer_months) %>%
      summarise(summer_average = mean(Value))
    permuted_T <- abs(permuted_winter$winter_average - permuted_summer$summer_average)
    p_values[i] <- mean(permuted_T >= T_stats[i])
  }
  
  # Compute confidence intervals
  z_score <- qnorm(1 - alpha/2)
  montecarlo_sd <- sqrt(p_values * (1 - p_values) / num_perm) 
  ci_upper <- T_stats + z_score * montecarlo_sd
  ci_lower <- T_stats - z_score * montecarlo_sd
  
  # Combine results into a table
  result <- winter_avg %>%
    select(ID, Name, winter_average) %>%
    left_join(summer_avg, by=c("ID", "Name")) %>%
    mutate(p_value=p_values, ci_lower=ci_lower, ci_upper=ci_upper)
  
  return(result)
}



#######################################################       
####################################################### 
####################################################### 
p_value_CI <- function(data, alpha=0.05, N=1000) {
  # Compute the p-value and Monte Carlo sd
  count <- rep(0, length = length(T_stats))
  for (n in 1:N) {
    # Generate a random permutation of the indices
    permutation <- sample(data$Summer, length(data$Summer), replace = FALSE)
    # Apply the permutation to the data and group by weather station
    ghcnd_permutation <- data %>%
      mutate("Summer" = permutation) %>%
      group_by(Name)
    # Compute the test statistic for the permuted data
    ghcnd_permutation_1 <- ghcnd_permutation %>%
      summarise(Value = abs(mean(Value[Summer]) - mean(Value[!Summer])), .groups = "drop")
    # Increment the count if the test statistic is greater than or equal to the observed value
    count <- count + (ghcnd_permutation_1$Value >= T_stats)
  }
  
  p_value <- count / N
  monte_carlo_sd <- sqrt(p_value * (1 - p_value) / N)
  
  # Compute confidence interval
  lower_CI <- T_stats - qnorm(1 - alpha / 2) * monte_carlo_sd
  upper_CI <- T_stats + qnorm(1 - alpha / 2) * monte_carlo_sd
  
  # Create result table
  result_table <- data.frame(Name = unique(winter_avg$Name),
                             T_stats = T_stats,
                             p_value = p_value,
                             monte_carlo_sd = monte_carlo_sd,
                             lower_CI = lower_CI,
                             upper_CI = upper_CI)
  return(result_table)
}




p_value_CI <- function(data, alpha=0.05, N=1000) {
  # Compute the p-value and Monte Carlo sd
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
    # Increment the count if the test statistic is greater than or equal to the observed value
    count <- count + (ghcnd_permutated$Value >= T_stats)
  }
  
  p_value <- count / N
  monte_carlo_sd <- sqrt(p_value * (1 - p_value) / N)
  
  # Compute confidence interval for p-value
  if (p_value == 0) {
    lower_CI <- 0
    upper_CI <- 1 - 0.025^(1/N)
  } else {
    lower_CI <- p_value - qnorm(1 - alpha/2) * monte_carlo_sd
    upper_CI <- p_value + qnorm(1 - alpha/2) * monte_carlo_sd
  }
  
  
  # Compute confidence interval for p-value
  CIp <- p_value + qnorm(1 - alpha/2) * monte_carlo_sd
  CIp[p_value == 0] <- c(0, 1 - 0.025^(1/N))
  
  
  # Create result table
  result_table <- data.frame(Name = unique(winter_avg$Name),
                             T_stats = T_stats,
                             p_value = p_value,
                             monte_carlo_sd = monte_carlo_sd,
                             lower_CI = lower_CI,
                             upper_CI = upper_CI)
  return(kable(result_table, align = "lccccc"))
}

p_value_CI(ghcnd_selected)



#######################################
#######################################
fit_models <- function(data, k) {
  models <- list()
  for (i in 0:k) {
    if(k == 0) {
      # Model M0
      formula <- Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear
      fit <- lm(formula, data = data)
      model <- summary(fit)
      return(model)
    } else{
      # Construct the formula for the model
      formula <- as.formula(paste("Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear", 
                                  paste0("+ cos(2*pi*", 1:i, "*DecYear) + sin(2*pi*", 1:i, "*DecYear)", collapse = ""), 
                                  sep = " "))
      # Fit the linear model
      model <- lm(formula, data = data)
      # Add the model to the list
      models[[i+1]] <- model
    }
  }
  # Return the summary of the models
  lapply(models, summary)
}

fit_models(ghcnd_spatial, 0)
fit_models(ghcnd_spatial, 1)
fit_models(ghcnd_spatial, 2)
fit_models(ghcnd_spatial, 3)
fit_models(ghcnd_spatial, 4)




fit_models <- function(data, k) {
  models <- list()
  for (i in 0:k) {
    if(k == 0) {
      # Model M0
      formula <- Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear
      fit <- lm(formula, data = data)
      model <- summary(fit)
      return(model)
    } else{
      # Construct the formula for the model
      formula <- as.formula(paste("Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear", 
                                  paste0("+ cos(2*pi*", i, "*DecYear) + sin(2*pi*", i, "*DecYear)", collapse = ""), 
                                  sep = " "))
      # Fit the linear model
      model <- lm(formula, data = data)
      # Add the model to the list
      models[[i+1]] <- model
    }
  }
  # Return the summary of the model with k = i
  for (i in 1:k) {
    print(summary(models[[i+1]]))
  }

}



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
  return(summary(fit))
}

# Fit the model for M0, M1, M2, M3, and M4
fit_model(data = ghcnd_spatial, k = 0)
fit_model(data = ghcnd_spatial, k = 1)
fit_model(data = ghcnd_spatial, k = 2)
fit_model(data = ghcnd_spatial, k = 3)
fit_model(data = ghcnd_spatial, k = 4)

##################################
#corss validation
##################################
# Stratified cross-validation function
stratified_cv <- function(data, model, n_folds = 5) {
  # Create a vector of unique station IDs
  stations <- unique(data$station_id)
  
  # Initialize a vector to store the cross-validation scores for each station and month
  cv_scores <- vector("list", length = length(stations))
  
  # Loop over each station
  for (i in seq_along(stations)) {
    # Subset the data for the current station
    station_data <- data[data$station_id == stations[i],]
    
    # Calculate the number of folds based on the number of observations for the current station
    n_obs <- nrow(station_data)
    fold_size <- floor(n_obs/n_folds)
    folds <- rep(1:n_folds, each = fold_size)
    if (n_obs %% n_folds != 0) {
      # Add remaining observations to last fold
      folds <- c(folds, rep(n_folds, n_obs - length(folds)))
    }
    
    # Initialize a vector to store the cross-validation scores for each fold
    fold_scores <- vector("list", length = n_folds)
    
    # Loop over each fold
    for (j in seq_along(fold_scores)) {
      # Subset the data for the current fold
      fold_data <- station_data[folds == j,]
      
      # Fit the model to the training data
      fit <- lm(model, data = station_data[!folds %in% j,])
      
      # Make predictions on the test data
      pred <- predict(fit, newdata = fold_data, se.fit = TRUE)
      
      # Calculate the cross-validation score for the current fold
      fold_score <- proper_score("ds", fold_data$y, pred$fit, sqrt((pred$se.fit)^2 + (pred$residual.scale)^2))
      
      # Store the fold score
      fold_scores[[j]] <- fold_score
    }
    
    # Average the fold scores and store the result for the current station and month
    cv_scores[[i]] <- tapply(unlist(fold_scores), station_data$month, mean)
  }
  
  # Compute the overall cross-validated average scores aggregated to the 12 months of the year
  overall_cv_scores <- tapply(unlist(cv_scores), factor(rep(1:12, length(cv_scores))), mean)
  
  # Return a list containing the cross-validation scores for each station and month, as well as the overall cross-validated average scores
  list(cv_scores = cv_scores, overall_cv_scores = overall_cv_scores)
}


###############################
##############################
library(tidyverse)
library(splitstackshape)
library(caret)
library(randomForest)

data=iris

## split data into train and test using stratified sampling
d <- rownames_to_column(data, var = "id") %>% mutate_at(vars(id), as.integer)
training <- d %>% stratified(., group = "Species", size = 0.90)
dim(training)

## proportion check
prop.table(table(training$Species)) 

testing <- d[-training$id, ]
dim(testing)
prop.table(table(testing$Species)) 


## Modelling

set.seed(123)

tControl <- trainControl(
  method = "cv", #cross validation
  number = 10, #10 folds
  search = "random" #auto hyperparameter selection
)




cross_validated_predictions <- function(data, form) { result <- data.frame(
  mean = numeric(nrow(data)),
  sd = numeric(nrow(data))
)
for (id in unique(data$ID)) {
  fit <- lm(form, data %>% filter(ID != id))
  pred <- predict(fit, newdata = data %>% filter(ID == id), se.fit = TRUE) result$mean[data$ID == id] <- pred$fit
  result$sd[data$ID == id] <- sqrt(pred$se.fitˆ2 + pred$residual.scaleˆ2)
}
result }


##################################
##################################
# fit model
fit_models <- function(k, data) {
  # Create a formula string for the initial model with no Fourier terms
  formula <- "Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear"
  
  # Loop through n and fit models with increasing numbers of Fourier terms
  for (i in 0:k) {
    if (i > 0) {
      # Add Fourier terms to the formula for models with more than 0 terms
      formula <- paste0(formula, "+ cos(", 2*i, "*pi*DecYear) + sin(", 2*i, "*pi*DecYear)")
    }
  }
  model <- lm(formula, data = data)
  return(model)
}


# formula define 
formula_model <- function(order = 1, extra = NULL) {
  form <- "Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear"
  if (order > 0) {
    form <- paste(c(form, paste0("cos", seq_len(order), "+sin", seq_len(order))), collapse = "+")
  }
  as.formula(paste0(c(form, extra), collapse = "+"))
}

add_harmonic <- function(dtad, order = 1) {
  for (k in seq_len(order)) {
    data[[paste0("cos", k)]] <- cos(2*pi*k*data$DecYear)
    data[[paste0("sin", k)]] <- sin(2*pi*k*data$DecYear)
  }
  data
}
# cross validation
stratified_cross_validation <- function(data, form) {
  result <- data.frame(mean = numeric(nrow(data)),
                       sd = numeric(nrow(data))
  )
  for (id in unique(data$ID)) {
    fit <- lm(form, data %>% filter(ID != id))
    pred <- predict(fit, newdata = data %>% filter(ID == id), se.fit = TRUE) 
    result$mean[data$ID == id] <- pred$fit
    result$sd[data$ID == id] <- sqrt(pred$se.fitˆ2 + pred$residual.scaleˆ2)
  }
  result
}

# calculate the score
all_score <- function(data, pred) {
  se <- proper_score("se", data$Value, mean = pred$mean)
  ds <- proper_score("ds", data$Value, mean = pred$mean, sd = pred$sd)
  return(data.frame(se, ds))
}

pred1 <- stratified_cross_validation(ghcnd_data, form)
all_score(ghcnd_data, pred1)



#####################################################
########## STRATIFIED CROSS VALIDATION ##############
#####################################################
stratified_cv <- function(data, model) {
  pred2 <- data %>% mutate(mean = NA_real_, sd = NA_real_)
  for (id in unique(data$ID)) {
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

# try some models to check whether the function defined for models below are correct
model0 <- Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear
score_result0 <- stratified_cv(ghcnd_data, model0)

model1 <- Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear + 
  sin(2*pi*DecYear) + cos(2*pi*DecYear)
score_result1 <- stratified_cv(ghcnd_data, model1)

model2 <- Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear + 
  sin(2*pi*DecYear) + cos(2*pi*DecYear) + sin(4*pi*DecYear) + cos(4*pi*DecYear)
score_result2 <- stratified_cv(ghcnd_data, model2)

# define the function for 5 models
modelform <- function(k){
  if(k == 0) { # When k = 0 
    model <- Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear
  } else{ # When k = 1,2,3,4
    model <- formula(paste0("Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear + ",
                            paste0("cos(2*pi*", seq(1, k), "*DecYear) + sin(2*pi*", seq(1, k), "*DecYear)", 
                                   collapse = " + ")))
  }
  return(model)
}

score_result222 <- stratified_cv(ghcnd_data, modelform(2))

# check the average of all ith month score
here222 <- score_result222 %>%
  group_by(Name, ID) %>%
  summarise(average_se = mean(se), average_ds = mean(ds), .groups = "drop")

#######################################################
#######################################################
# combine two functions together
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
# check correct BINGO
score_together <- stratified_cv(ghcnd_data, 2)

##############################################
# Joint the data together

### Prediction scores for EACH station
# Create an empty list to store the data frames

cv_list <- list()
# Loop through the values of k and run the stratified_cv function
for (k in 0:4) {
  cv_data <- stratified_cv(ghcnd_data, k)
  cv_data$k <- k
  cv_list[[k+1]] <- cv_data
}
# Combine the data frames
cv_all <- bind_rows(cv_list)

### OVERALL cross-validated average scores
overall_average <- cv_all %>%
  group_by(k, Name, ID) %>%
  summarise(average_se = mean(se), average_ds = mean(ds), .groups = "drop")

#############################################
################## PLOT #####################
#############################################
# Plot the average corss validation scores for each station.
# Is the prediction accuracy the same across the whole year?

# Define a vector of colors
colors <- c("darkblue", "darkturquoise", "deeppink", "#F0E442", "springgreen3")
# DS docre
ggplot(data = cv_all) +
  geom_point(aes(x = Month, y = ds, color = factor(k))) +
  facet_wrap(~Name, ncol = 4) +
  labs(title = "Dawid-Sebastiani scores", x = "Month", y = "DS", color = "Model") +
  scale_color_manual(labels = c("M0", "M1", "M2", "M3", "M4"), values = colors) +
  scale_x_continuous(breaks = 1:12) +
  theme(axis.text.x = element_text(size = 6))
# SE socre
ggplot(data = cv_all) +
  geom_point(aes(x = Month, y = se, color = factor(k))) +
  facet_wrap(~Name, ncol = 4) +
  labs(title = "Standard Errors scores", x = "Month", y = "SE", color = "Model") +
  scale_color_manual(labels = c("M0", "M1", "M2", "M3", "M4"), values = colors) +
  scale_x_continuous(breaks = 1:12) +
  theme(axis.text.x = element_text(size = 6))

### Is the model equally good at predicting the different stations?
# Create a plot to compare the average DS scores across stations and models
DS_plot <- ggplot(data = overall_average, aes(x = k, y = average_ds, color = Name)) +
  geom_point() +
  geom_line() +
  labs(title = "Overall Average DS", x = "Model", y = "DS score", color = "Station") +
  scale_x_continuous(breaks = 0:4, labels = paste0("M", 0:4)) +
  theme(plot.title = element_text(size = 12),
        legend.text=element_text(size=6),
        axis.text.x=element_text(size=7))
# Create a plot to compare the average SE scores across stations and models
SE_plot <- ggplot(data = overall_average, aes(x = k, y = average_se, color = Name)) +
  geom_point() +
  geom_line() +
  labs(title = "Overall Average SE", x = "Model", y = "SE score", color = "Station") +
  scale_x_continuous(breaks = 0:4, labels = paste0("M", 0:4)) +
  theme(plot.title = element_text(size = 12),
        legend.text=element_text(size=6),
        axis.text.x=element_text(size=7))
DS_plot + SE_plot + plot_layout(ncol = 2, guides = "collect")



#######################################################
#Techer Ding's format
#######################################################
PRCP_month <- rainfall_data %>%
  group_by(ID, Year, Month) %>%
  summarise(ID=ID, Year=Year, Month=Month, DecYear_avg=mean(DecYear), Value_sqrt_avg=sqrt(mean(Value)),
            Longitude=Longitude, Latitude=Latitude, Elevation=Elevation, .groups = "drop") %>% 
  distinct()

fit_models <- function(n, df) {
  # Create a formula string for the initial model with no Fourier terms
  formula <- "Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear_avg"
  
  # Loop through n and fit models with increasing numbers of Fourier terms
  for (i in 0:n) {
    if (i > 0) {
      # Add Fourier terms to the formula for models with more than 0 terms
      formula <- paste0(formula, "+ cos(", 2*i, "*pi*DecYear_avg) + sin(", 2*i, "*pi*DecYear_avg)")
    }
  }
  model <- lm(formula, data = df)
  return(model)
}

summary(fit_models(0, PRCP_month))
summary(fit_models(1, PRCP_month))
summary(fit_models(2, PRCP_month))
summary(fit_models(3, PRCP_month))
summary(fit_models(4, PRCP_month))


single_score <- function(fit, pred, v_data, n, score_type){
  residual_var <- sum(fit$residuals^2)/fit$df.residual
  sd_pred <- sqrt(pred$se.fit^2 + residual_var)
  if (score_type=="se"){
    score <- proper_score(type = "se", obs = v_data$Value_sqrt_avg,
                          mean = pred$fit)
  }
  else if (score_type=="ds"){
    score <- proper_score(type = "ds", obs = v_data$Value_sqrt_avg,
                          mean = pred$fit, sd = sd_pred)
  }
  return(score)
}

score_result <- function(K, score_type){
  # create an empty data frame to store the results
  score_df <- data.frame(ID = c(),
                         Year = integer(),
                         Month = integer(),
                         score_se = numeric())
  
  # loop through each unique combination of ID, Year, and Month
  for (i in unique(PRCP_month$ID)){
    PRCP_station <- PRCP_month %>% filter(ID == i)
    t_data <- PRCP_month %>% filter(!ID == i)
    fit <- fit_models(K, df = t_data)
    #monthly_score_se <- c()
    for (j in unique(PRCP_station$Year)){
      PRCP_station_year <- PRCP_station %>% filter(Year == j)
      for (k in unique(PRCP_station_year$Month)){
        v_data <- PRCP_station_year %>% filter(Month == k)
        pred <- predict(fit, newdata = v_data, se.fit = TRUE)
        score_se <- single_score(fit, pred, v_data, K, score_type)
        score_df <- rbind(score_df, data.frame(ID = i,
                                               Year = j,
                                               Month = k,
                                               score_se = score_se))
      }
    }
  }
  avg_score <- score_df %>%
    group_by(ID, Month) %>%
    summarise(ID=ID, Month=Month, score_se = mean(score_se), .groups = "drop") %>% 
    distinct()
  return(avg_score)
}

ds1 <- score_result(1,"ds")






pred2 <- ghcnd_data %>% mutate(mean = NA_real_, sd = NA_real_)
for (id in unique(ghcnd_data$ID)) {
  if(k == 0) { # When k = 0 
    model <- Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear
  } else{ # When k = 1,2,3,4
    model <- formula(paste0("Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear + ",
                            paste0("cos(2*pi*", seq(1, k), "*DecYear) + sin(2*pi*", seq(1, k), "*DecYear)", 
                                   collapse = " + ")))
  }
  # fit and predict the model for the data
  fit <- lm(model, ghcnd_data %>% filter(ID != id))
  pred <- predict(fit, newdata = ghcnd_data %>% filter(ID == id), se.fit = TRUE)
  # evaluate the mean and sd
  pred2$mean[(pred2$ID == id)] <- pred$fit
  pred2$sd[(pred2$ID == id)] <- sqrt((pred$se.fit)^2 + (pred$residual.scale)^2)
  # mutate the table
  pred2 <- pred2 %>%
    group_by(ID, Name, Year, Month) %>%
    mutate(se = proper_score("se", Value_sqrt_avg, mean = mean),
           ds = proper_score("ds", Value_sqrt_avg, mean = mean, sd = sd)) }
Asss<- pred2 %>%
  group_by(ID, Name, Month) %>%
  summarise(se = mean(se), ds = mean(ds), .groups = "drop")

