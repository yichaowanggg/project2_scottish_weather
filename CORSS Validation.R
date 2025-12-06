###QUESTION2###
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


model0 <- stratified_cv(ghcnd_data,0)
model1 <- stratified_cv(ghcnd_data,1)
model2 <- stratified_cv(ghcnd_data,2)
model3 <- stratified_cv(ghcnd_data,3)
model4 <- stratified_cv(ghcnd_data,4)

