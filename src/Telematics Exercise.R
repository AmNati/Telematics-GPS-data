pacman::p_load(data.table,
               tidyverse,
               mlr, # for Model building
               lubridate,
               ggmap, # Plotting spatial data
               geosphere, # for Haversine distance
               kableExtra,
               factoextra
)


#' Setting Working directory path and loading data
path <- getwd() %>%
  setwd()

#' GPS data
gps <- fread("data/sample_trips.csv") %>%
  as.data.frame()


#' Coercing local_dtm to POSIXct time
gps <- gps %>%
  mutate(date.time = dmy_hms(local_dtm))

#' Checking for duplicates
# gps[duplicated.data.frame(gps)]

# Removing duplicates
gps <- gps %>%
  distinct(longitude, latitude, .keep_all = TRUE)
#' Structure of the data
glimpse(gps)


#' Setting Working directory path and loading data
path <- getwd() %>%
  setwd()

#' GPS data
gps <- fread("data/sample_trips.csv") %>%
  as.data.frame()


#' Coercing local_dtm to POSIXct time
gps <- gps %>%
  mutate(date.time = dmy_hms(local_dtm))

#' Checking for duplicates
# gps[duplicated.data.frame(gps)]

# Removing duplicates
gps <- gps %>%
  distinct(longitude, latitude, .keep_all = TRUE)
#' Structure of the data
glimpse(gps)



map_plot_func = function(t_n, df, x = "longitude", y = "latitude"){
  dat <- df %>% filter(trip_nb == t_n)
  lon <- median(as.matrix(dat[[x]]))
  lat <- median(as.matrix(dat[[y]]))
  map_g <- get_map(location = c(lon, lat), maptype = "terrain",
                   source = "google",
                   crop = TRUE,
                   scale = 1,
                   zoom = 14)
  
  ggmap(map_g, base_layer = ggplot(dat, aes(x = .data[[x]], y = .data[[y]]))) +
    geom_point() +
    ggtitle(paste("Trip", t_n)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

cowplot::plot_grid(map_plot_func(t_n = 1, gps), map_plot_func(t_n = 2, gps),     
                   map_plot_func(t_n = 3, gps), map_plot_func(t_n = 4, gps))

cowplot::plot_grid(map_plot_func(t_n = 1, gps), map_plot_func(t_n = 2, gps), ... =    
                     map_plot_func(t_n = 3, gps), map_plot_func(t_n = 4, gps), 
                   map_plot_func(t_n = 5, gps), map_plot_func(t_n = 6, gps),
                   map_plot_func(t_n = 7, gps), map_plot_func(t_n = 8, gps)
)

cowplot::plot_grid(map_plot_func(t_n = 9, gps), map_plot_func(t_n = 10, gps), ... =    
                     map_plot_func(t_n = 11, gps), map_plot_func(t_n = 12, gps), 
                   map_plot_func(t_n = 13, gps), map_plot_func(t_n = 14, gps),
                   map_plot_func(t_n = 15, gps), map_plot_func(t_n = 16, gps)
)

# Anomaly detection with DBSCAN
Dbscan_func = function(d.f, tn){
  dat <- d.f %>% filter(trip_nb == tn)
  db <- fpc::dbscan(as.matrix(dat$longitude, dat$latitude), 
                    eps = 0.15, MinPts = 5)
  
  fviz_cluster(db, data.frame(dat$longitude, dat$latitude),
               stand = FALSE, elipse = FALSE, geom = "point" ,palette = "jco") +
    ggtitle(paste("Trip", tn)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}



cowplot::plot_grid(Dbscan_func(gps, tn = 1), Dbscan_func(gps, tn = 2), 
                   Dbscan_func(gps, tn = 3), Dbscan_func(gps, tn = 4),
                   Dbscan_func(gps, tn = 5), Dbscan_func(gps, tn = 6),
                   Dbscan_func(gps, tn = 7), Dbscan_func(gps, tn = 8)
)


cowplot::plot_grid(Dbscan_func(gps, tn = 9), Dbscan_func(gps, tn = 10), 
                   Dbscan_func(gps, tn = 11), Dbscan_func(gps, tn = 12),
                   Dbscan_func(gps, tn = 13), Dbscan_func(gps, tn = 14),
                   Dbscan_func(gps, tn = 15), Dbscan_func(gps, tn = 16)
)



## Removing anomalous points from trip 1 and 8 and replacing with NA
#' Adding a row_index as column
gps <- gps %>%
  mutate(row_index = as.numeric(row.names(gps)), .before = trip_nb)

#' Function to obtain row indices of outliers
anomalous_indices_func <- function(tri.no, df){
  trip_df <- df %>% filter(trip_nb == tri.no)
  trip_dbscan <- fpc::dbscan(as.matrix(trip_df$longitude,
                                       trip_df$latitude), eps = 0.01, MinPts = 10)
  
  ind <- which(trip_dbscan$cluster == 0) # indices of outliers
  row.indices <- trip_df[ind, ]$row_index
  return(row.indices)
}

anomalous_points <- c(anomalous_indices_func(1, gps), anomalous_indices_func(8,gps))

#' Replacing anomalous points with NA
gps_df <- gps %>%
  mutate(longitude = replace(longitude, anomalous_points , NA),
         latitude = replace(latitude, anomalous_points, NA)
  )





# Interpolating missing instances to fill gaps
gps_1 <- gps_df %>%
  select(trip_nb, longitude, latitude, date.time)

interpo_func <- function(var,dt) approx(dt,var,xout=seq(min(dt),max(dt), by = 1))$y
gps_interpo <- setDT(gps_1)[,lapply(.SD,interpo_func,
                                    dt = date.time), by = list(trip_nb),
                            .SDcols = c("latitude","longitude","date.time")]


get.dist <- function(lon, lat) distHaversine(tail(cbind(lon,lat),-1), p1 =  
                                               head(cbind(lon,lat),-1))

# Calculating distances
gps_interpo_df <- gps_interpo %>%
  group_by(trip_nb) %>%
  mutate(date.time = as.POSIXct(date.time,
                                origin = "1970-01-0100:00:00", tz = "UTC"),
         time.diff = as.numeric(difftime(date.time,
                                         lag(date.time, default = first(date.time)),
                                         units = "secs")))

# Calculating distance, velocity and acceleration
gps_interpo_df_vel <- gps_interpo_df %>%
  group_by(trip_nb) %>%
  mutate(distance = c(0, get.dist(longitude,latitude)), velocity = 3.6*distance, 
         acceleration = c(0,diff(velocity,1)))


#' Determining Thresholds of Hard Events

#' Threshold for hard brakes
thresh_brake <- fivenum(gps_interpo_df_vel$acceleration)[2] - 
  1.5*IQR(gps_interpo_df_vel$acceleration)

#' Threshold for hard acceleration
thresh_acc <- fivenum(gps_interpo_df_vel$acceleration)[4] +          
  1.5*IQR(gps_interpo_df_vel$acceleration)



ggplot(gps_interpo_df_vel, aes(x = acceleration)) +
  geom_histogram(bins = 40, colour="black", fill="white") +
  labs(title = "Distribution of acceleration events", x  = " acceleration km/h.sq      ")+
  theme_bw()+
  xlim(-40,40) + 
  geom_vline(xintercept = c(-20,20), linetype="dotted", 
             color = "red", size=1.5)


Events_summary <- gps_interpo_df_vel %>%
  group_by(trip_nb) %>%
  summarise(Distance = round(sum(distance),0),
            HardAccelerations = sum(acceleration > 20),
            HardBrakes = sum(acceleration < -20),
            Idling = sum(acceleration == 0))
Events_summary %>%
  kable(caption = 'Summary of Events Per Trip',
        col.names = c("Trip No.","Distance (meters)","Hard Accelerations", "Hard Brakes", "Idlings"),
        booktabs = TRUE) %>%
  kable_classic_2(full_width = T)


## Modeling

#' Splitting the simulated Telematics data into Training and Test sets

telem <- fread("data/simulated_summary_total.csv") %>%
  as.data.frame()


train_index <- sample(1:nrow(telem), 0.8 * nrow(telem))
test_index <- setdiff(1:nrow(telem), train_index)

telem_train <- telem[train_index,-c(1,2)]
telem_test <- telem[test_index, -c(1,2)]


#' Exploratory data Analysis


plt_1 <- ggplot() + 
  geom_boxplot(aes(y = telem_train$Distance)) +
  labs(title = " Total distance in meters ", y = "Distance") +
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

plt_2 <- ggplot() + 
  geom_boxplot(aes(y = telem_train$NightTime_Pct)) +
  labs(title = "% of distance driven at Night", y = " Distance ") +
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

plt_3 <-ggplot() + 
  geom_boxplot(aes(y = telem_train$HardAccelerations)) +
  labs(title = "Hard Accelerations Event", y = "Hard Accelerations") +
  theme_bw()+
  ylim(0,2000)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

plt_4 <- ggplot() + 
  geom_boxplot(aes(y = telem_train$HardBrakes)) +
  labs(title = "Hard Brakes Event", y = "Hard Brakes ") +
  theme_bw()+
  ylim(0,2000)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


cowplot:: plot_grid(plt_1, plt_2, plt_3, plt_4)



#Contigency Table
tab <- telem%>% 
  group_by(VehicleType, Loss) %>%
  summarise(n = n()) %>%
  spread(Loss, n)


tab %>%
  kable(caption = 'Contigency Table of Vehicle Type and Class',
        col.names = c("Vehicle Type","No collision", "Collision"),
        booktabs = TRUE) %>%
  kable_classic_2(full_width = T) 


# Target Variable
# Categorical Variables
ggplot(data = telem_train, aes(as.factor(Loss), fill = as.factor(Loss)))+
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  labs(x = "Loss", y = "") +
  scale_y_continuous(labels=scales::percent) +
  guides(fill = guide_legend(title="Loss")) +
  scale_fill_discrete(
    labels = c("No collision", "Collision"))



#' quantile distribution of numeric variables
fivenum_sum_acc <- fivenum(telem_train$HardAccelerations)
fivenum_sum_Brake <- fivenum(telem_train$HardBrakes)
fivenum_sum_night_drive <- fivenum(telem_train$NightTime_Pct)

telem_train <- telem_train %>%               
  mutate(Hard_Acc_category = cut(HardAccelerations, 
                                 breaks = c(-Inf,fivenum_sum_acc[2],Inf), 
                                 labels = c("lowRisk","HighRisk")),
         .after= HardAccelerations) %>%
  
  mutate(Hard_Brake_category = cut(HardBrakes, 
                                   breaks = c(-Inf,fivenum_sum_Brake[2],Inf), 
                                   labels = c("lowRisk","HighRisk")),
         .after= HardBrakes) %>%
  
  mutate( Night_Time_category = cut(NightTime_Pct, 
                                    breaks =   c(-Inf,fivenum_sum_night_drive[4],Inf), 
                                    labels = c("lowRisk","HighRisk")),.after= NightTime_Pct)


#' 
telem_test <- telem_test %>% 
  mutate(Hard_Acc_category = cut(HardAccelerations,
                                 breaks = c(-Inf,fivenum_sum_acc[2], Inf), 
                                 labels = c("lowRisk","HighRisk")),
         .after= HardAccelerations) %>%
  
  mutate(Hard_Brake_category = cut(HardBrakes,
                                   breaks = c(-Inf,fivenum_sum_Brake[2], Inf),
                                   labels = c("lowRisk","HighRisk")),.after= HardBrakes) %>%
  
  mutate(Night_Time_category = cut(NightTime_Pct,
                                   breaks = c(-Inf,fivenum_sum_night_drive[4],Inf),
                                   labels = c("lowRisk","HighRisk")),.after= NightTime_Pct) 


train_df <- telem_train %>% 
  select(!c("HardAccelerations", "HardBrakes", "NightTime_Pct"))%>%
  mutate(VehicleType = as.factor(VehicleType), Loss = as.factor(Loss))


test_df <- telem_test %>% 
  select(!c("HardAccelerations", "HardBrakes", "NightTime_Pct")) %>%
  mutate(VehicleType = as.factor(VehicleType), Loss = as.factor(Loss))


#' Capping and normalizing Distance

cap <- fivenum(telem_train$Distance)[4] + 1.5*IQR(telem_train$Distance) 

train_df <- capLargeValues(train_df, target = "Loss", 
                           cols = c("Distance"),
                           threshold = cap)

test_df <- capLargeValues(test_df, target = "Loss",
                          cols = c("Distance"),
                          threshold = cap)

# Normalize Train set
train_df <- normalizeFeatures(train_df,
                              target = "Loss",
                              method = "scale",
                              cols = NULL,
                              on.constant = "quiet"
) 

# Normalize Test set
test_df <- normalizeFeatures(test_df,
                             target = "Loss",
                             method = "scale",
                             cols = NULL,
                             on.constant = "quiet"
) 


#' Creating Learning Task
#' Unbalanced class
trainTask <- makeClassifTask(data = train_df,
                             target = "Loss",
                             positive = "1")

testTask <- makeClassifTask(data = test_df,
                            target = "Loss", 
                            positive = "1")

#' Balancing class using different sampling technique
train_over <-  oversample(trainTask, rate = 6)
test_over  <-  oversample(testTask, rate = 6)

train_under <- undersample(trainTask, rate = 1/6)
test_under <- undersample(testTask, rate = 1/6)

train_smote <- smote(trainTask, rate = 6, nn = 8)
test_smote <- smote(testTask, rate = 6, nn = 8)




fv = generateFilterValuesData(train_smote, method = "FSelectorRcpp_information.gain")

plotFilterValues(fv) +
  labs(title = "Feature importance")




# Making a logistic Learner

logistic_learner <- makeLearner("classif.logreg", 
                                predict.type = "prob",
                                fix.factors.predict = TRUE,
                                predict.threshold = 0.40)

# cross validation (cv) accuracy
cv.logistic <- crossval(learner = logistic_learner,
                        task = train_over,iters = 15L, 
                        stratify = TRUE,
                        measures = list(acc ,tpr),
                        show.info = F)

#cv.logistic$measures.test

#training model
log_fit <- train(logistic_learner, train_over)
#getLearnerModel(log_model)

#predicting on test data
log_pred <- predict(log_fit, test_over)

perf_metrics <- performance(log_pred, measures = list(tpr,ppv, acc, auc,f1))

# Random Forest
param_set_ranf <- getParamSet("classif.randomForest")

# create a learner
ranf_learner <- makeLearner("classif.randomForest", 
                            predict.type = "prob", 
                            par.vals = list(ntree = 100, mtry = 2),
                            predict.threshold = 0.4)

ranf_learner$par.vals <- list(importance = TRUE)

#' setting  tunable parameters
#' grid search to find hyperparameters

ranf_param <- makeParamSet(
  makeIntegerParam("ntree",lower = 50, upper = 100),
  makeIntegerParam("mtry", lower = 2, upper = 10),
  makeIntegerParam("nodesize", lower = 10, upper = 50)
)

#' random search for 50 iterations
rand_search <- makeTuneControlRandom(maxit = 40L )

#' Cross-validation
cv_rand <- makeResampleDesc("CV",iter =3L , predict = "both")

#hypertuning
ranf_tune <- tuneParams(learner = ranf_learner, 
                        task = train_over, 
                        par.set = ranf_param, 
                        resampling = cv_rand, 
                        control = rand_search, 
                        measures = list(tpr, ppv, acc))


#cv accuracy
ranf_tune$y
ranf_tune$x #optimal parameters

#using hyperparameters for modeling
ranf_tree <- setHyperPars(ranf_learner, par.vals = ranf_tune$x)


rand_forest_fit <- train(ranf_tree, train_over)

rand_pred <- predict(rand_forest_fit, test_over)

#' Performance metrics
rand_perf <- performance(rand_pred, list(tpr,ppv, acc, auc,f1))


gb_param_set <- getParamSet("classif.gbm")

gb_learner <- makeLearner("classif.gbm",
                          predict.type = "prob",
                          predict.threshold = 0.4)

#specify tuning method
gb_rand_search <- makeTuneControlRandom(maxit = 50L)

#3 fold cross validation
cv_gb <- makeResampleDesc("CV",iters = 3L)

#parameters
gb_param <- makeParamSet(
  makeDiscreteParam("distribution", values = "bernoulli"),
  makeIntegerParam("n.trees", lower = 10, upper = 100), 
  makeIntegerParam("interaction.depth", lower = 2, upper = 5), 
  makeIntegerParam("n.minobsinnode", lower = 2, upper = 8),
  makeNumericParam("shrinkage",lower = 0.01, upper = 1)
)

#tune parameters
gb_tune <- tuneParams(learner = gb_learner, 
                      task = train_smote,
                      resampling = cv_gb,
                      measures = list(acc ,tpr),
                      par.set = gb_param,
                      control = gb_rand_search)



#check CV accuracy
gb_tune$y

#set parameters
gb_optim <- setHyperPars(learner = gb_learner, par.vals = gb_tune$x)

#train
gb_fit <- train(gb_optim, train_smote)

#test
gb_pred <- predict(gb_fit, test_smote)

performance(gb_pred, list(tpr,ppv, acc, auc,f1))


# Hot encoding
train_xg <- createDummyFeatures(train_smote)
test_xg <- createDummyFeatures(test_smote)

xgb_learner <- makeLearner("classif.xgboost", 
                           predict.type = "prob",
                           predict.threshold = 0.4)

xgb_learner$par.vals <- list(
  objective = "binary:logistic",
  eval_metric = "error",
  nrounds = 250 )

#define parameters for tuning
xgb_param <- makeParamSet(
  makeIntegerParam("nrounds",lower=200,upper=600),
  makeIntegerParam("max_depth",lower=3,upper=20),
  makeNumericParam("lambda",lower=0.55,upper=0.60),
  makeNumericParam("eta", lower = 0.001, upper = 0.5),
  makeNumericParam("subsample", lower = 0.10, upper = 0.80),
  makeNumericParam("min_child_weight",lower=1,upper=5),
  makeNumericParam("colsample_bytree",lower = 0.2,upper = 0.8)
)


#define search function
xgb_rand_search <- makeTuneControlRandom(maxit = 50L) #do 100 iterations

#3 fold cross validation
cv_xgb <- makeResampleDesc("CV",iters = 3L)

#tune parameters
xgb_tune <- tuneParams(learner = xgb_learner,
                       task = train_xg, 
                       resampling = cv_xgb,
                       measures = list(acc ,tpr),
                       par.set = xgb_param, 
                       control = xgb_rand_search)

#set parameters
xgb_optim <- setHyperPars(learner = xgb_learner, par.vals = xgb_tune$x)

#train model
xgb_fit <- train(xgb_optim, train_xg)

#test model
xgb_pred <- predict(xgb_fit, test_xg)
performance(xgb_pred, list(ppv, tpr, acc, auc, f1, gmean))
xgb_Thresh <- setThreshold(xgb_pred, threshold = c("1" = .4, "0" = .60 ))
performance(xgb_Thresh, measures = list(ppv, tpr, acc, auc, f1))


