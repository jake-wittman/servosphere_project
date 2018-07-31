# Libraries ---------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(doParallel)
pacman::p_load(ggplot2)
pacman::p_load(dplyr)
pacman::p_load(magrittr)
pacman::p_load(purrr)
pacman::p_load(stringr)
pacman::p_load(car)
pacman::p_load(zoo)
pacman::p_load(agricolae)
pacman::p_load(beepr)
pacman::p_load(nlme)
pacman::p_load(future.apply)
pacman::p_load(data.table)
pacman::p_load(glmmTMB)
pacman::p_load(multcomp)
#pacman::p_load(gridExtra)
set.seed(1)

# Data cleaning -----------------------------------------------------------
# Import all data
fullpath <- list.files(path = "./data/", pattern = "_servosphere.csv", full.names = T) # get full file names
# Uncomment below for a subset of fullpath for development
# fullpath <- fullpath[1:5]
plan(multiprocess, workers = 3)
dat_files <- future_lapply(fullpath, fread) # put all file names in a list for R
dat_files <- future_lapply(dat_files, as.data.frame)
trial_num <- fullpath %>% # extract trial numbers with this pipeline
  basename() %>%
  strsplit(split = "_")
trial_num <- unlist(trial_num)[c(T, F, F)] # remove extraneous text
names(dat_files) <- trial_num
for (i in 1:length(dat_files)) { # add trial number to each df
  dat_files[[i]]$id <- trial_num[i]
}

trial_record <- read.csv("data/food_trial_record.csv", na.strings = ".") # import record of trials
names(trial_record) <- c("date", "id", "food", "mass", "loops",
                         "include", "noted_dist", "notes")
trial_record <- trial_record[which(trial_record$id %in% as.numeric(trial_num)), ] # get rid of light trial records

# split food and starvation into separate variables in trial record, calculate actual starvation time
trial_record$food <- ifelse(trial_record$food == "fed", "diet_0",
                            ifelse(trial_record$food == "starved_24", "diet_24",
                                   ifelse(trial_record$food == "starved_48", "diet_48",
                                          as.character(trial_record$food))))
group_ids <- unlist(strsplit(as.character(trial_record$food), "_"))
group_ids <- ifelse(group_ids == "fed", 0, group_ids) # need to change "fed" to 0 for time since feeding
trial_record$food <- group_ids[seq(1, length(group_ids), by = 2)]
trial_record$starve <- as.numeric(group_ids[seq(2, length(group_ids), by = 2)]) # create starve variable
mtimes <- strftime(file.info(fullpath)$mtime, format = "%H:%M")
hm <- as.numeric(unlist(strsplit(mtimes, ":"))) # get hours and minutes
hm_df <- data.frame(h = hm[seq(1, length(hm), 2)], m = hm[seq(2, length(hm), 2)])
trial_record <- cbind(trial_record, hm_df)
trial_record$starve_time <- ifelse(trial_record$starve == 24, (24 + (trial_record$h - 9)),
                                   ifelse(trial_record$starve == 48, (48 + (trial_record$h - 9)),
                                          0))

#clean up a little more

trial_record <- trial_record[, c("date", "id", "food", "starve", "noted_dist", "mass", "loops", "include",
                                 "starve_time", "notes")]
trial_record$food <- as.factor(trial_record$food)
discard <- which(trial_record$include == 0)
keep_id <- as.character(trial_record$id[-discard])
dat_files <- dat_files[keep_id]
trial_record_include <- trial_record[trial_record$include == 1, ]

# Clean up files a bit
colnames <- c("stimulus", "dT", "dx", "dy", "enc1", "enc2", "enc3", "id") 
dat_files <- lapply(dat_files, setNames, colnames) #
Changetype <- function(df){
  # Change dx and dy to numeric variables 
  df$dx <- as.numeric(as.character(df$dx))
  df$dy <- as.numeric(as.character(df$dy))
  df
}
dat_files <- future_lapply(dat_files, Changetype)

# Fed vs Starved trials - drop stimulus 1
# For trial numbers < 397, drop stimulus 0 and 1 (1 - 81 in list)
# For trial numbers > 398 drop stimulus 0. Not sure how the others ended up how they did (> 81 in list)
# Droprows <- function(x){
#   x <- x[x$stimulus == 2, ]
# }
# dat_files <- lapply(dat_files[[1:81]], function(df) {
#   df <- df[df$stimulus == 2, ]
# })
for (i in 1:length(dat_files)) {
  if (i <= 61) {
    dat_files[[i]] <- dat_files[[i]][dat_files[[i]]$stimulus == 2, ]
  } else {
    dat_files[[i]] <- dat_files[[i]][dat_files[[i]]$stimulus == 1, ]
  }
}

# Check time for each stimulus (should be roughly the same)
SumdT <- function(x) {sum(x$dT)/1000/60} # function for use in sapply
sapply(dat_files, FUN = SumdT) #return a vector of times for each trial to confirm they're roughly correct

# Thin files (https://stackoverflow.com/questions/30359427/calculate-the-mean-of-every-13-rows-in-data-frame)
Thin <- function(dat, n){ # Function to use with thinning the data
  # For each selected column, sum all n consecutive rows (i.e. rows 1:10, 11:20, etc.) to condense the
  # dataset and reduce noise. [-1] removes unecessary first column from output
  # The second aggregate is to verify that each new row is a sum of n rows (or the remainder when you get to the end)
  x <- cbind(
    aggregate(dat[, c("dT", "dx", "dy")],
              list(rep(1:(nrow(dat) %/% n + 1), each = n, len = nrow(dat))),
              sum)[, -1],
    aggregate(dat[, c("dT")],
              list(rep(1:(nrow(dat) %/% n + 1), each = n, len = nrow(dat))),
              length)[, -1])
  colnames(x)[4] <- "length"
  x$id <- unique(dat$id)
  return(x)
}
dat_files_agg <- future_lapply(dat_files, Thin, n = 100) # thin the data
dat_files_agg <- future_lapply(dat_files_agg, cbind, x = NA, y = NA) # add empty columns for calculating x and y coords for path
dat_files_agg <- lapply(dat_files_agg, # Merge the trial_record information with the servosphere data
                        function(x, y) {merge(x, y)},
                        y = trial_record_include)

# Calculate (x, y) coords.

# Single core method
# for (j in 1:length(dat_files_agg)) {
#   for (i in 1:nrow(dat_files_agg[[j]])) {
#     if (i == 1){
#       dat_files_agg[[j]]$x[i] <- 0
#       dat_files_agg[[j]]$y[i] <- 0
#     } else {
#       dat_files_agg[[j]]$x[i] <- dat_files_agg[[j]]$x[i-1] + dat_files_agg[[j]]$dx[i]
#       dat_files_agg[[j]]$y[i] <- dat_files_agg[[j]]$y[i-1] + dat_files_agg[[j]]$dy[i]
#     }
#   }
# }

# (x,y) coords ------------------------------------------------------------

# Parallel method for calculating (x, y) coords.
registerDoParallel(cores = detectCores() - 1)
# Send each data frame to a different core for processing in foreach loop.
dat_files_agg <- foreach(j = 1:length(dat_files_agg)) %dopar% { 
  temp_dat <- dat_files_agg[[j]]
  for (i in 1:nrow(temp_dat)) {
    if (i == 1){
      temp_dat$x[i] <- 0 # set start position as (0,0)
      temp_dat$y[i] <- 0
    } else {
      temp_dat$x[i] <- temp_dat$x[i-1] + temp_dat$dx[i] # calculate new (x,y) 
      temp_dat$y[i] <- temp_dat$y[i-1] + temp_dat$dy[i]
    }
  }
  temp_dat # call temp_dat to get output
}
closeAllConnections()
#names(dat_files_agg) <- trial_num


# Plot all paths  ---------------------------------------------------------
# Set theme for plotting function
simple_theme <- theme(panel.background=element_rect(fill=NA), 
                      panel.border = element_rect(colour = "black",
                                                  fill=NA, size=1),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=12),
                      legend.position = "none",
                      strip.background = element_rect(fill = NA))
PlotPath <- function(df){ 
  # takes a servosphere data file and plots the path
  # title is a vector of titles to put on the plot to help identify the plot
  if (df$include[1] == 1){
    plot <- ggplot(dat = df, aes(x = x, y = y))
    plot + geom_path(size = 1.2) + ggtitle(paste(df$id[1]))
  }
}

plotlist <- lapply(dat_files_agg, PlotPath) # create plots
# save plots as one pdf to look at them
#ggsave("movementplots.pdf", marrangeGrob(grobs = plotlist, nrow = 2, ncol = 2))
#names(plotlist) <- trial_num


# Calculating individual derived variables -------------------------------------------

# Calculate bearing & ta
dat_files_agg <- lapply(dat_files_agg,
                        cbind,
                        ang = NA, ta = NA, turn_velocity = NA, vi = NA)
registerDoParallel(cores = detectCores() - 1)
dat_files_agg <- foreach(k = 1:length(dat_files_agg)) %dopar% {
  temp_dat <- dat_files_agg[[k]]
  for (i in 1:nrow(temp_dat)) {
    if (i == 1) {
      temp_dat$ang[i] <- atan2(temp_dat$x[i], temp_dat$y[i])
      temp_dat$turn_velocity[i] <- NA
    } else {
      temp_dat$ang[i] <-
        atan2(temp_dat$x[i] - temp_dat$x[i - 1], temp_dat$y[i] - temp_dat$y[i - 1])
      # calculate bearing in radians
    }
    
    if (temp_dat$ang[i] > 0) {
      temp_dat$bearing[i] <-
        temp_dat$ang[i] * (180 / pi) # converts radians to degrees
    } else if (temp_dat$ang[i] < 0) {
      temp_dat$bearing[i] <- (temp_dat$ang[i] + 2 * pi) * (180 / pi)
      # converts negative radians to positive and then to degrees
    } else {
      temp_dat$bearing[i] <- NA
    }
    
    if (i > 1) {
      # Calculate turn angle
      if (is.na(temp_dat$bearing[i])) {
        temp_dat$ta[i] <- NA
      } else if (!is.na(temp_dat$bearing[i]) & !is.na(temp_dat$bearing[i - 1])){
        temp_dat$ta[i] <- (temp_dat$bearing[i] - temp_dat$bearing[i - 1]) * (pi / 180)
      } else if (is.na(temp_dat$bearing[i - 1])) {
        # if last bearing was NA need to keep going back till you find a non-NA bearing
        j <- 1
        while (is.na(temp_dat$bearing[i - j])) {
          j <- j + 1
          if( j == i) {
            break
          }
        }
        if (j == i){
          temp_dat$ta[i] <- NA
        } else {
          temp_dat$ta[i] <- (temp_dat$bearing[i] - temp_dat$bearing[i - j]) * (pi / 180)
        }
        
        # convert to radians
        # Keep TA between -180 and 180
        if (is.na(temp_dat$ta[i])) {
          temp_dat$ta[i] <- NA
        } else if (temp_dat$ta[i] > pi) {
          temp_dat$ta[i] <- (temp_dat$ta[i] - (2 * pi)) * (180 / pi)
        } else if (temp_dat$ta[i] < -pi) {
          temp_dat$ta[i] <- (temp_dat$ta[i] + (2 * pi)) * (180 / pi)
        } else {
          temp_dat$ta[i] <- temp_dat$ta[i] * (180 / pi)
        }
      }
    }
  }
  temp_dat
}

closeAllConnections()

# Circular mean function to get mean and concentration, see next comment for more info
RollCMean <- function(x){
  m <- suppressWarnings(circular::mean.circular(circular::as.circular(x,
                                                                      type = "angles",
                                                                      units = "degrees",
                                                                      zero = pi/2,
                                                                      rotation = "clock"),
                                                na.rm = T)[[1]])
  
  if (m < 0 & !is.na(m)) { # change negative angles to positive
    m <- 360 + m
  }
  return(m)
}

# Not currently using this function
# Apply the RollCmean function as a rolling mean to reduce noise
# By using a rolling mean on the bearings, I reduce noise
# created by wobbles in the insect movement. Width of rolling mean
# set by width = n, where n = # of seconds to use
# Rollbear <- function(y){
#   bear <- y$bearing
#   res <- rollapply(bear, width = 1, FUN = RollCMean, fill = NA, by = 1)
#   y$roll_bear <- res
#   return(y)
# }
# 
# dat_files_agg <- lapply(dat_files_agg, Rollbear)


# The above rolling mean functions calculate the rolling mean for bearing
# This function calculates the circular average and concentration for bearing for each trial
CMean <- function(x){
  b <- x$bearing
  # Get average movement direction, m
  m <- suppressWarnings(circular::mean.circular(circular::as.circular(b,
                                                                      type = "angles",
                                                                      units = "degrees",
                                                                      zero = pi/2,
                                                                      rotation = "clock"),
                                                na.rm = T)[[1]])
  # Get rho for average movement direction
  r <- suppressWarnings(circular::rho.circular(circular::as.circular(b,
                                                                     type = "angles",
                                                                     units = "degrees",
                                                                     zero = pi/2,
                                                                     rotation = "clock"),
                                               na.rm = T)[[1]])
  
  if (m < 0 & !is.na(m)) { # change negative angles to positive
    m <- 360 + m
  }
  if (is.nan(r)) {
    r <- NA
  }
  return(c(m, r))
}

# Turn velocity
dat_files_agg <- lapply(dat_files_agg, function(df) {
  df$turn_velocity <- abs(df$ta / (df$dT / 1000))
  return(df)
})

# Calculate instantaneous velocity
InstVel <- function(df){
  df$vi <- ifelse(df$x == 0 & df$y == 0,
                  0,
                  (sqrt((df$dx ^ 2) + (df$dy ^ 2))) / (df$dT / 1000))
  return(df)
}
dat_files_agg <- lapply(dat_files_agg, InstVel)


# Whole trial variables ---------------------------------------------------

# Calculate total distance
TotalDistance <- function(df){ 
  # calculate total distance moved for each data frame (1 trial)
  td <- sum(sqrt((df$dx ^ 2) + (df$dy ^ 2)))
  return(c(df$id[1], td))
}
tot_dist <- sapply(dat_files_agg, TotalDistance) # return total distance moved for each data frame
summary_df <- as.data.frame(x = t(tot_dist)) # start putting together a summary data frame
names(summary_df) <- c("id", "tot_dist") # rename columns
summary_df$tot_dist <- as.numeric(as.character(summary_df$tot_dist))

# Calculate displacement and tortuosity
summary_df$displacement <- sapply(dat_files_agg, function(x){sqrt(x$x[nrow(x)] ^ 2 + x$y[nrow(x)] ^ 2)})
summary_df$tortuosity <- summary_df$displacement / summary_df$tot_dist


# Average individual trial bearing
avg_bearing <- sapply(dat_files_agg, CMean)
summary_df$avg_bearing <- avg_bearing[1, ]
summary_df$rho <- avg_bearing[2, ]
summary_df$avg_turnangle <- sapply(dat_files_agg, function(df){mean(df$ta, na.rm = T)})

# Velocity averages
summary_df$avg_vi <- sapply(dat_files_agg, function(x) mean(x$vi)) 
summary_df$avg_moving_vi <- sapply(dat_files_agg, function(x) mean(x$vi[x$vi != 0])) #average velocity when moving (velocities > 0)
summary_df$avg_turn_v <- sapply(dat_files_agg, function(x) mean(x$turn_velocity, na.rm = T))

# Time spent walking (in seconds)
summary_df$walk_time <- sapply(dat_files_agg, function(df) {
  sum(df$dT[df$dx != 0 | df$dy != 0]) / 1000
})

# Number of stops
numStops <- function(df) {
  stops <- ifelse(df$vi <= 0.1, 0, df$vi) # movements with a velocity less than 0.1 are not true movements
  n_s <- rle(stops) # get length of sequences of numbers in the above vector
  return(length(n_s$lengths[n_s$values == 0])) # how many sequences of 0 velocity are there
}
summary_df$stops <- sapply(dat_files_agg, numStops)

# length of stops
lStops <- function(df) {
  stops <- ifelse(df$vi <= 0.1, 0, df$vi) # movements with a velocity less than 0.1 are not true movements
  n_s <- rle(stops) # get length of sequences of numbers in the above vector
  return(mean(n_s$lengths[n_s$values == 0])) # how many sequences of 0 velocity are there
}
summary_df$l_stops <- sapply(dat_files_agg, lStops)

names(dat_files_agg) <- keep_id

# Group summaries ---------------------------------------------------------
n_total <- trial_record %>% 
  group_by(food, starve) %>% 
  summarize(count = n())
n_10 <- trial_record %>%
  filter(include == 1) %>% 
  group_by(food, starve) %>% 
  summarize(count = n())
n_total <- as.data.frame(n_total)
n_10 <- as.data.frame(n_10)
n_10$total <- n_total$count
n_10$percent <- n_10$count / n_total$count

loops <- trial_record %>% 
  filter(loops == 1) %>% 
  group_by(food, starve) %>% 
  summarize(count = n())
sum(loops$count) # total number of "loops" paths
sum(n_10$count) # total number of larvae that went > 10 cm
37 / 244 # % loops

summary_df <- merge(summary_df, trial_record) #add in the trial record information
summary_df <- summary_df[summary_df$include == 1, ] # Need to remove larvae that went < 10 cm
summary_table <- summary_df %>%
  group_by(food, starve) %>%
  summarize(count = n(), 
            starve_time = mean(starve_time, na.rm = T),
            disp = mean(displacement, na.rm = T), 
            dist = mean(tot_dist, na.rm = T),
            tortuosity = mean(tortuosity, na.rm = T),
            sd.tort = mean(tortuosity, na.rm = T) / sqrt(n()),
            bearing = circular::mean.circular(
              circular::as.circular(avg_bearing,
                                    type = "angles",
                                    units = "degrees",
                                    zero = pi /
                                      2,
                                    rotation = "clock"
              ),
              na.rm = T
            ),
            rho = mean(rho, na.rm = T),
            avg_vi = mean(avg_vi, na.rm = T),
            avg_moving_vi = mean(avg_moving_vi, na.rm = T),
            avg_turn_v = mean(avg_turn_v, na.rm = T),
            avg_stops = mean(stops, na.rm = T),
            avg_lstops = mean(l_stops, na.rm = T),
            avg_walk_time = mean(walk_time, na.rm = T),
            avg_ta = mean(avg_turnangle, na.rm = T)
  )
summary_table$bearing <- ifelse(summary_table$bearing < 0,
                                360 + summary_table$bearing,
                                summary_table$bearing)
as.data.frame(summary_table)

# Summarize results just by food treatment
food_summary <- summary_df %>%
  group_by(food) %>%
  summarize(count = n(), 
            dist = mean(tot_dist, na.rm = T),
            se.dist = sd(tot_dist, na.rm = T) / sqrt(n()),
            disp = mean(displacement, na.rm = T),
            se.disp = sd(displacement, na.rm = T) / sqrt(n()),
            avg.tort = mean(tortuosity, na.rm = T),
            se.tort = sd(tortuosity, na.rm = T) / sqrt(n()),
            bearing = circular::mean.circular(
              circular::as.circular(avg_bearing,
                                    type = "angles",
                                    units = "degrees",
                                    zero = pi /
                                      2,
                                    rotation = "clock"
              ),
              na.rm = T
            ),
            rho = mean(rho, na.rm = T),
            avg.vi = mean(avg_vi, na.rm = T),
            se.vi = sd(avg_vi, na.rm = T) / sqrt(n()),
            avg_moving_vi = mean(avg_moving_vi, na.rm = T),
            avg_turn_v = mean(avg_turn_v, na.rm = T),
            avg_stops = mean(stops, na.rm = T),
            se.stops = sd(stops, na.rm = T) / sqrt(n()),
            avg_walk_time = mean(walk_time, na.rm = T),
            avg_ta = mean(avg_turnangle, na.rm = T)
  )
food_summary$bearing <- ifelse(food_summary$bearing < 0,
                               360 + food_summary$bearing,
                               food_summary$bearing)
as.data.frame(food_summary)


# Visualize data ----------------------------------------------------------
hist(dat_files_agg[[1]])
# Boxplots of summary_df variables
ggplot(summary_df, aes(x = starve_time, y = tot_dist)) +
  geom_point() +
  facet_wrap( ~ food) +
  simple_theme 

ggplot(summary_df, aes(x = starve_time, y = displacement)) +
  geom_point() +
  facet_wrap( ~ food) +
  simple_theme 

ggplot(summary_df, aes(x = starve_time, y = tortuosity)) +
  geom_point() +
  facet_wrap( ~ food) +
  simple_theme 

ggplot(summary_df, aes(x = starve_time, y = avg_vi)) +
  geom_point() +
  facet_wrap( ~ food) +
  simple_theme 

ggplot(summary_df, aes(x = starve_time, y = avg_moving_vi)) +
  geom_point() +
  facet_wrap( ~ food) +
  simple_theme 

ggplot(summary_df, aes(x = starve_time, y = stops)) +
  geom_point() +
  facet_wrap( ~ food) +
  simple_theme

ggplot(summary_df, aes(x = starve_time, y = walk_time)) +
  geom_point() +
  facet_wrap( ~ food) +
  simple_theme

ggplot(summary_df, aes(x = starve_time, y = l_stops)) +
  geom_point() +
  facet_wrap( ~ food) + 
  simple_theme

ggplot(summary_df, aes(x = food, y = l_stops)) +
  geom_boxplot() +
  simple_theme

ggplot(summary_df, aes(x = starve_time, y = avg_turn_v)) +
  geom_point() +
  facet_wrap( ~ food) +
  simple_theme

ggplot(summary_df, aes(x = starve_time, y = avg_turnangle)) +
  geom_point() +
  facet_wrap( ~ food) +
  simple_theme

ggplot(trial_record, aes(x = starve_time, y = loops)) +
  geom_point() +
  facet_wrap( ~ food) +
  simple_theme

# Histograms for tortuosity
ggplot(summary_df, aes(x = tortuosity)) +
  geom_histogram() +
  facet_wrap( ~ food) +
  simple_theme

# Tortuosity vs Displacement
ggplot(summary_df, aes(x = tot_dist, y = tortuosity)) +
  geom_point() +
  facet_wrap( ~ food) +
  simple_theme

ggplot(summary_df, aes(x = avg_vi, y = tot_dist)) +
  geom_point() +
  facet_wrap( ~ food) +
  simple_theme

# Plotting all data
dat_all <- do.call("rbind", dat_files_agg)
dat_all$food <- as.factor(dat_all$food)
dat_all$id <- as.factor(dat_all$id)

ggplot(dat_all, aes(x = starve_time, y = vi)) +
  geom_point(aes(colour = id), alpha = 0.2, position = "jitter") +
  facet_wrap( ~ food) +
  simple_theme +
  theme(legend.position = "none")

dat_all %>% 
  dplyr::filter(vi > 0) %>% 
  ggplot(aes(x = starve_time, y = vi)) +
  geom_point(aes(colour = id), alpha = 0.2, position = "jitter") +
  facet_wrap( ~ food) +
  simple_theme +
  theme(legend.position = "none")

ggplot(dat_all, aes(x = bearing)) +
  geom_histogram() +
  facet_wrap( ~ food) +
  simple_theme +
  theme(legend.position = "none")

ggplot(dat_all, aes(x = ta)) +
  geom_histogram() +
  facet_wrap( ~ food) +
  simple_theme



# Plots of average bearing
ggplot(summary_df, aes(x = avg_bearing)) +
  geom_histogram(bins = 36) +
  facet_wrap( ~ food) +
  simple_theme

# rose diagrams for feeding trials
par(mfrow = c(1, 2))
pacman::p_load(circular)
# Not sure why this is commented out...keeping for now.
# for(i in unique(summary_df$food)) {
#   temp_dat <-
#     summary_df[summary_df$food == i, ]
#   temp_dat$avg_bearing <- as.circular(
#     temp_dat$avg_bearing,
#     type = "angles",
#     units = "degrees",
#     zero = pi /
#       2,
#     rotation = "clock"
#   )
#   for (j in unique(temp_dat$food)) {
#     
#     temp_dat2 <- temp_dat[temp_dat$food == j,]
#     plot(
#       temp_dat2$avg_bearing,
#       stack = T,
#       bin = 72,
#       shrink = 1.3
#     )
#     title(paste(toString(i), toString(j)), line = 0)
#     axis.circular(
#       at = circular(seq(pi / 4, 7 * pi / 4, pi / 2)),
#       zero = pi / 2,
#       rotation = "clock",
#       cex = 1.1
#     ) #add intercardinal direction labels
#     rose.diag(
#       temp_dat2$avg_bearing,
#       bins = 24,
#       col = "darkgrey",
#       cex = 1.0,
#       prop = 1.8,
#       add = TRUE
#     ) #add rose diagram
#   }
# }

par(mfrow = c(1, 2))
for(i in unique(summary_df$food)) {
  temp_dat <-
    summary_df[summary_df$food == i,]
  temp_dat$avg_bearing <- as.circular(
    temp_dat$avg_bearing,
    type = "angles",
    units = "degrees",
    zero = pi /
      2,
    rotation = "clock"
  )
  plot(
    temp_dat$avg_bearing,
    stack = T,
    bin = 72,
    shrink = 1.3
  )
  title(paste(toString(i)), line = 0)
  axis.circular(
    at = circular(seq(pi / 4, 7 * pi / 4, pi / 2)),
    zero = pi / 2,
    rotation = "clock",
    cex = 1.1
  ) #add intercardinal direction labels
  rose.diag(
    temp_dat$avg_bearing,
    bins = 24,
    col = "darkgrey",
    cex = 1.0,
    prop = 1.8,
    add = TRUE
  ) #add rose diagram
}

detach(package:circular)



# Statistical Analysis -----------------------------------------------
interceptCompare <- function(model, data) {
  results <- list() # empty list to store results
  s <- "food" #ADD THIS LINE
  
  for (i in unique(data$food)) {
    updatedfactor <- relevel(data$food, ref = i)
    model[["model"]][[s]] <- updatedfactor #CHANGE THE DATA IN THE MODEL
    model <- update(model,data=model[["model"]])# UPDATE THE MODEL
    results[[i]] <- summary(model)$coefficients[1:5, ]
  }
  results <- lapply(results, function(x) round(x, 4))
  return(results)
}


# Differences in mass
lm_mass <- lm(mass ~ starve_time * food, data = summary_df)
par(mfrow = c(2,2))
plot(lm_mass)
ggplot(summary_df, aes(x = starve_time, y = mass)) +
  geom_point() +
  facet_wrap( ~ food) +
  simple_theme
summary(lm_mass)
anova(lm_mass)
# No significant differences in mass

# Probability an insect will move for trial
lm_prob_move <- glm(include ~ food / starve_time,
                    data = trial_record,
                    family = binomial)

plot(lm_prob_move)
fit_prob <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 100), 4))
fit_prob$food <- rep(unique(summary_df$food), each = 100)
pred_prob <- predict.glm(lm_prob_move, newdata = fit_prob, type = "link", se.fit = T)
fit_prob$predicted <- plogis(pred_prob$fit)
fit_prob$low_ci <- plogis(pred_prob$fit - 1.96 * pred_prob$se.fit)
fit_prob$high_ci <- plogis(pred_prob$fit + 1.96 * pred_prob$se.fit)
ggplot(data = trial_record, aes(x = starve_time, y = include)) +
  geom_point() +
  geom_line(data = fit_prob,
            aes(x = starve_time, y = predicted),
            size = 1.7,
            inherit.aes = F) +
  geom_ribbon(data = fit_prob, 
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) +
  facet_wrap( ~ food) +
  simple_theme
summary(lm_prob_move)
Anova(lm_prob_move)
interceptCompare(lm_prob_move, trial_record)
summary(glht(lm_prob_move, linfct = mcp(food = "Tukey", interaction_average = T, covariate_average = T)))


# In all cases, the larvae were more likely to move as time since last feeding increased. Larch fed larvae were had a greater than 50% probability of moving and bur oak had less than 50% probability of moving. The interaction of diet and starve time was significant, which I think suggests that the effect of starvation on the probability of moving was not as great for larvae fed diet.

# Probability of loops
lm_prob_loops <- glm(loops ~ starve_time * food,
                     family = binomial,
                     data = trial_record_include)
fit_loops <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 100), 4))
fit_loops$food <- rep(unique(summary_df$food), each = 100)
pred_loops <- predict.glm(lm_prob_loops, newdata = fit_loops, type = "link", se.fit = T)
fit_loops$predicted <- plogis(pred_loops$fit)
fit_loops$low_ci <- plogis(pred_loops$fit - 1.96 * pred_loops$se.fit)
fit_loops$high_ci <- plogis(pred_loops$fit + 1.96 * pred_loops$se.fit)
ggplot(data = trial_record, aes(x = starve_time, y = loops)) +
  geom_point() +
  geom_line(data = fit_loops,
            aes(x = starve_time, y = predicted),
            size = 1.7,
            inherit.aes = F) +
  geom_ribbon(data = fit_loops, 
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) +
  facet_wrap( ~ food) +
  simple_theme
summary(lm_prob_loops)
Anova(lm_prob_loops)
interceptCompare(lm_prob_loops, trial_record_include)
summary(glht(lm_prob_loops, linfct = mcp(food = "Tukey", interaction_average = T, covariate_average = T)))

# ANOVA means comparisons without starvation ------------------------------
# We decided not to use these ANOVA models and instead do means comparisons with the glm models below after averaging over the interaction and covariate effects.
par(mfrow = c(2,2))
# mass
lm(mass ~ food, data = subset(summary_df, starve == 0)) %>% 
  plot()
anova_mass <- lm(mass ~ food, data = subset(summary_df, starve == 0))
Anova(anova_mass)
summary(glht(anova_mass, linfct = mcp(food = "Tukey"))) # post-hoc tests

# Probability of moving
anova_prob_move <- glm(include ~ food ,
                       data = subset(trial_record, starve == 0),
                       family = binomial)

plot(anova_prob_move)
Anova(anova_prob_move)
summary(anova_prob_move)
plogis(summary(anova_prob_move)$coefficients)
plogis(summary(anova_prob_move)$coefficients[, 1] + 1.96 * summary(anova_prob_move)$coefficients[, 2])
plogis(summary(anova_prob_move)$coefficients[, 1] - 1.96 * summary(anova_prob_move)$coefficients[, 2])
summary(glht(anova_prob_move, linfct = mcp(food = "Tukey"))) # post-hoc tests


# Total Distance
lm(tot_dist ~ food, data = subset(summary_df, starve == 0)) %>% 
  plot()
anova_glm_totdist <- glm(tot_dist ~ food,
                         data = subset(summary_df, starve == 0),
                         family = Gamma(link = "identity"))
plot(anova_glm_totdist)
Anova(anova_glm_totdist)
summary(anova_glm_totdist)
summary(glht(anova_glm_totdist, linfct = mcp(food = "Tukey"))) # post-hoc tests
ggplot(summary_df, aes(x = food, y = tot_dist)) +
  geom_boxplot() +
  simple_theme

# Displacement
lm(displacement ~ food, data = subset(summary_df, starve == 0)) %>% 
  plot()
anova_glm_disp <- glm(displacement ~ food,
                      data = subset(summary_df, starve == 0),
                      family = Gamma(link = "identity"))
plot(anova_glm_disp)
Anova(anova_glm_disp)
summary(anova_glm_disp)
summary(glht(anova_glm_disp, linfct = mcp(food = "Tukey"))) # post-hoc tests
ggplot(summary_df, aes(x = food, y = displacement)) +
  geom_boxplot() +
  simple_theme

# Tortuosity
# I'm not sure what to do about tortuosity, can't seem to find a model that doesn't have heterogeneous variance.
lm(tortuosity ~ food, data = subset(summary_df, starve == 0)) %>% 
  plot()
ggplot(summary_df, aes(x = food, y = tortuosity)) +
  geom_boxplot() +
  simple_theme
anova_tort <- lm(asin(sqrt(tortuosity)) ~ food, data = subset(summary_df, starve == 0))
Anova(anova_tort)
summary(glht(anova_tort, linfct = mcp(food = "Tukey"))) # post-hoc tests

# Velocity
lm(avg_vi  ~ food, data = subset(summary_df, starve == 0)) %>% 
  plot()
anova_glm_vi <- glm(avg_vi ~ food,
                    data = subset(summary_df, starve == 0),
                    family = Gamma(link = "identity"))
plot(anova_glm_vi)
Anova(anova_glm_vi)
summary(anova_glm_vi)
summary(glht(anova_glm_vi, linfct = mcp(food = "Tukey"))) # post-hoc tests
ggplot(summary_df, aes(x = food, y = avg_vi)) +
  geom_boxplot() +
  simple_theme

# Number of stops
lm(stops ~ food, data = subset(summary_df, starve == 0)) %>% 
  plot()
anova_glm_stops <- glm.nb(stops ~ food, data = subset(summary_df, starve == 0))
plot(anova_glm_stops) # better
Anova(anova_glm_stops)
summary(anova_glm_stops)
ggplot(summary_df, aes(x = food, y = stops)) +
  geom_boxplot() +
  simple_theme
summary(glht(anova_glm_stops, linfct = mcp(food = "Tukey"))) # post-hoc tests

# Length of stops
lm(l_stops ~ food, data = subset(summary_df, starve == 0)) %>% 
  plot()
anova_glm_l_stops <- glm(l_stops ~ food,
                         data = subset(summary_df, starve == 0),
                         family = Gamma(link = "identity"))
plot(anova_glm_l_stops) # pretty bad...
summary(glht(anova_glm_l_stops, linfct = mcp(food = "Tukey"))) # post-hoc tests

# Plot ANOVA
summary_df %>% 
  filter(starve == 0 ) %>% 
  ggplot(aes(x = food, y = stops)) +
    geom_boxplot() +
    simple_theme

# Linear models -----------------------------------------------------------
# Total distance and displacement show increasing variance in residuals for larger fitted values, so trying glm
# Tried fitting some no-interaction models - thye were inferior
lm_totdist <- glm(tot_dist ~ starve_time + food,
                  data = summary_df,
                  family = Gamma(link = "identity"))

lm_disp <- glm(displacement ~ starve_time + food,
               data = summary_df,
               family = "Gamma")

lm_tort2 <- glm(tortuosity ~ starve_time * food,
                data = summary_df,
                family = binomial)
lm_rho <- lm(asin(sqrt(rho)) ~ starve_time * food, data = summary_df)

lm_stops <- glm(stops ~ starve_time * food,
                data = summary_df,
                family = poisson) # not great, try a negative binomial
# negbin meets assumptions much better (see below)
lm_length_stops <- lm(l_stops ~ starve_time * food, data = summary_df)
glm_length_stops <- glm(l_stops ~ starve_time * food,
                        data = summary_df,
                        family = Gamma(link = "identity"))
lm_walktime <- lm(walk_time ^ 2 ~ starve_time * food,
                  data = summary_df)



par(mfrow = c(2, 2))
# Residuals look okay for the two three gamma models with identity link

plot(lm_rho) # also an interaction model


plot(lm_length_stops)
plot(glm_length_stops) # still not great.
plot(lm_walktime)


# Total distance
lm_totdist_int <- glm(tot_dist ~ food / starve_time, 
                      data = summary_df,
                      family = Gamma(link = "identity"))
plot(lm_totdist_int)
anova(lm_totdist, lm_totdist_int, test = "Chisq")
# The interactions model explains more deviance
# get fitted values for plotting
fit_dist <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_dist$food <- rep(unique(summary_df$food), each = 1000)
pred_dist <- predict.glm(lm_totdist_int, newdata = fit_dist, type = "link", se = T)
fit_dist$predicted <- pred_dist$fit
fit_dist$low_ci <- (fit_dist$predicted - 1.96 * pred_dist$se.fit)
fit_dist$high_ci <-(fit_dist$predicted + 1.96 * pred_dist$se.fit)
ggplot(summary_df, aes(x = starve_time, y = tot_dist)) +
  geom_point() +
  geom_line(data = fit_dist, aes(x = starve_time, y = predicted), size = 1.7) +
  geom_ribbon(data = fit_dist,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) + 
  facet_wrap( ~ food) +
  simple_theme 
summary(lm_totdist_int)
interceptCompare(lm_totdist_int, summary_df)
Anova(lm_totdist_int)
summary(glht(lm_totdist_int, linfct = mcp(food = "Tukey", interaction_average = T, covariate_average = T))) # post-hoc tests


# With the exception of larch, as starve time increased so did the total distance moved in 10 minutes. Starvation time seemed to decrease the amount of distance moved by larch. (Discussion: seems to be supported by probability of movement. Larch larvae had a higher probability of moving to begin with.) Larvae fed diet seem to have moved significantly less than those on other diets as well. Possibly seeing the effect of better nutrition overall from bur oak and diet increasing ability to move further when starved? In other words, larvae that feed on larch or norway maple are less likely to move further when starved possibly because they can't move as fast? Larvae that went further also went faster - upper limit on velocity when feeding on lower quality food?

# Displacement
lm_disp_int <- glm(displacement ~ food / starve_time,
                   data = summary_df,
                   family = Gamma(link = "identity"))
plot(lm_disp_int)
anova(lm_disp, lm_disp_int, test = "Chisq")
# The interactions model explains more  deviance
summary(lm_disp_int)
# Get fitted values for plotting
fit_displacement <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_displacement$food <- rep(unique(summary_df$food), each = 1000)
pred_disp <- predict.glm(lm_disp_int, newdata = fit_displacement, type = "link", se = T)
fit_displacement$predicted <- pred_disp$fit
fit_displacement$low_ci <- (fit_displacement$predicted - 1.96 * pred_disp$se.fit)
fit_displacement$high_ci <-(fit_displacement$predicted + 1.96 * pred_disp$se.fit)
ggplot(summary_df, aes(x = starve_time, y = displacement)) +
  geom_point() +
  geom_line(data = fit_displacement, aes(x = starve_time, y = predicted), size = 1.7) +
  geom_ribbon(data = fit_displacement,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) + 
  facet_wrap( ~ food) +
  simple_theme 
summary(lm_disp_int)
interceptCompare(lm_disp_int, summary_df)
Anova(lm_disp_int)
summary(glht(lm_disp_int, linfct = mcp(food = "Tukey", interaction_average = T, covariate_average = T))) # post-hoc tests
# Results are similar to total distance moved, although norway maple appears to have a significant increase in displacement as starvation time increases (paths become straighter?). Larch shows a negative trend in displacement with increasing starvation time. Diet on average displace less than the others. The negative trend for displacement in larch holds as in total distance. Discussion: possibly due to the nutrition available in larch? While they seem to develop as well as the Norway maple. I think there's a paper by Weseloh that showed that larvae on pines are more likely to disperse.

# Tortuosity
lm_tort <- lm(asin(sqrt(tortuosity)) ~ food / starve_time, data = summary_df)
plot(lm_tort) # also an interaction model
fit_tort <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_tort$food <- rep(unique(summary_df$food), each = 1000)
pred_tort <- predict(lm_tort, newdata = fit_tort, type = "response", se = T)
fit_tort$predicted <- sin(pred_tort$fit) ^ 2
fit_tort$low_ci <- sin(pred_tort$fit - (1.96 * pred_tort$se.fit)) ^ 2
fit_tort$high_ci <- sin(pred_tort$fit + (1.96 * pred_tort$se.fit)) ^ 2
ggplot(summary_df, aes(x = starve_time, y = tortuosity)) +
  geom_point() +
  geom_line(data = fit_tort, aes(x = starve_time, y = predicted), inherit.aes = F, size = 2) +
  geom_ribbon(data = fit_tort,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) +
  facet_wrap( ~ food) +
  simple_theme 
summary(lm_tort)
Anova(lm_tort)
interceptCompare(lm_tort, summary_df) # Doesn't work, so do it by hand
summary_df$food <- relevel(summary_df$food, ref = "diet")
lm_tort <- update(lm_tort)
round(summary(lm_tort)$coefficients[1:4, ], 4)
summary_df$food <- relevel(summary_df$food, ref = "larch")
lm_tort <- update(lm_tort)
round(summary(lm_tort)$coefficients[1:4, ], 4)
summary_df$food <- relevel(summary_df$food, ref = "norway")
lm_tort <- update(lm_tort)
round(summary(lm_tort)$coefficients[1:4, ], 4)
summary_df$food <- relevel(summary_df$food, ref = "buroak")
lm_tort <- update(lm_tort)
round(summary(lm_tort)$coefficients[1:4, ], 4)

summary(glht(lm_tort, linfct = mcp(food = "Tukey", interaction_average = T, covariate_average = T))) # post-hoc tests
# There's a lot of variation in this data but in general, bur oak and norway have higher tortuosity as starvation increases (paths are straighter) than diet and larch, which have negative slopes.

# Velocity
lm_vi_int <- glm(avg_vi ~ food / starve_time,
                 data = summary_df,
                 family = Gamma(link = "identity")) 
plot(lm_vi_int)
fit_vi <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_vi$food <- rep(unique(summary_df$food), each = 1000)
pred_vi <- predict.glm(lm_vi_int, newdata = fit_vi, type = "link", se = T)
fit_vi$predicted <- pred_vi$fit
fit_vi$low_ci <- (fit_vi$predicted - 1.96 * pred_vi$se.fit)
fit_vi$high_ci <-(fit_vi$predicted + 1.96 * pred_vi$se.fit)
ggplot(summary_df, aes(x = starve_time, y = avg_vi)) +
  geom_point() +
  geom_line(data = fit_vi, aes(x = starve_time, y = predicted), inherit.aes = F, size = 2) +
  geom_ribbon(data = fit_vi,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) +
  facet_wrap( ~ food) +
  simple_theme
summary(lm_vi_int)
Anova(lm_vi_int)
interceptCompare(lm_vi_int, summary_df)
summary(glht(lm_vi_int, linfct = mcp(food = "Tukey", interaction_average = T, covariate_average = T)))# post-hoc tests
# Velocity increases as starvation increases for bur oak and diet but stays flat for norway and has a slight decrease for larch. Supports total distance moved

# Turn velocity
ggplot(summary_df, aes(x = starve_time, y = avg_turn_v)) +
  geom_point() +
  facet_wrap( ~ food) +
  simple_theme

# Stops
negbin_stops <- glm.nb(stops ~ food / starve_time,
                       data = summary_df)
plot(negbin_stops)
fit_stops <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_stops$food <- rep(unique(summary_df$food), each = 1000)
pred_stops <- predict.glm(negbin_stops, newdata = fit_stops, type = "response", se = T)
fit_stops$predicted <- (pred_stops$fit)
fit_stops$low_ci <- (fit_stops$predicted - 1.96 * pred_stops$se.fit)
fit_stops$high_ci <- (fit_stops$predicted + 1.96 * pred_stops$se.fit)
ggplot(summary_df, aes(x = starve_time, y = stops)) +
  geom_point() +
  geom_line(data = fit_stops, aes(x = starve_time, y = predicted), inherit.aes = F, size = 2) +
  geom_ribbon(data = fit_stops,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) +
  facet_wrap( ~ food) +
  simple_theme
summary(negbin_stops)
Anova(negbin_stops)
interceptCompare(negbin_stops, summary_df)
cld(glht(negbin_stops, linfct = mcp(food = "Tukey", interaction_average = T, covariate_average = T)))

# Length of stops
fit_l_stops <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_l_stops$food <- rep(unique(summary_df$food), each = 1000)
pred_l_stops <- predict.glm(glm_length_stops, newdata = fit_l_stops, type = "response", se = T)
fit_l_stops$predicted <- (pred_l_stops$fit)
fit_l_stops$low_ci <- (fit_l_stops$predicted - 1.96 * pred_l_stops$se.fit)
fit_l_stops$high_ci <- (fit_l_stops$predicted + 1.96 * pred_l_stops$se.fit)
ggplot(summary_df, aes(x = starve_time, y = l_stops)) +
  geom_point() +
  geom_line(data = fit_l_stops, aes(x = starve_time, y = predicted), inherit.aes = F, size = 2) +
  geom_ribbon(data = fit_l_stops,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) +
  facet_wrap( ~ food) +
  simple_theme
summary(glm_length_stops)
Anova(glm_length_stops)
summary(glht(glm_length_stops, linfct = mcp(food = "Tukey")))

# Walk time
fit_walk <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_walk$food <- rep(unique(summary_df$food), each = 1000)
pred_walk <- predict(lm_walktime, newdata = fit_walk, type = "response", se = T)
fit_walk$predicted <- sqrt(pred_walk$fit)
fit_walk$low_ci <- sqrt(fit_walk$predicted) - 1.96 * sqrt(pred_walk$se.fit)
fit_walk$high_ci <- sqrt(fit_walk$predicted) + 1.96 * sqrt(pred_walk$se.fit)
ggplot(summary_df, aes(x = starve_time, y = walk_time)) +
  geom_point() +
  geom_line(data = fit_walk, aes(x = starve_time, y = predicted), inherit.aes = F, size = 2) +
  geom_ribbon(data = fit_walk,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) +
  facet_wrap( ~ food) +
  simple_theme
summary(lm_walktime)
summary(glht(lm_walktime, linfct = mcp(food = "Tukey")))

# Random effects models
lm_vi2 <- lme(vi ~ starve_time * food, data = dat_all, random = ~1 | id)
summary(lm_vi2)
plot(lm_vi2)

# GLM Mixed effects - can only use moving velocities (no 0s)

if (file.exists("data/glmm_vi.rds")) {
  glmm_vi <- readRDS("data/glmm_vi.rds")
} else {
  glmm_vi <- glmmTMB(vi ~ starve_time * food + (1|id), data = dat_all[dat_all$vi > 0, ], family = Gamma, verbose = T)
}

summary(glmm_vi)


# Publication plots: example paths -------------------------------------------------------
# Each "fitted" dataframe has a "significance" variable to change line type
# For labelling subplots, use this labeller function

# Movement paths
straight_path_plot <- ggplot(dat = dat_files_agg[["520"]], aes(x = x, y = y)) +
  geom_path(size = 0.7) +
  coord_cartesian(xlim = c(-65, 65), ylim = c(-65, 65)) +
  scale_x_continuous(name = "Distance (cm)", breaks = seq(-60, 60, 30), labels = seq(-60, 60, 30)) +
  scale_y_continuous(name = "Distance (cm)", breaks = seq(-60, 60, 30), labels = seq(-60, 60, 30)) +
  labs(tag = "a") +
  theme(plot.tag = element_text(size = 16, face = "bold")) + 
  simple_theme
straight_path_plot

loop_path_plot <- ggplot(dat_files_agg[["510"]], aes(x = x, y = y)) +
  geom_path(size = 0.7) +
  coord_cartesian(xlim = c(-65, 65), ylim = c(-65, 65)) +
  scale_x_continuous(name = "Distance (cm)", breaks = seq(-60, 60, 30), labels = seq(-60, 60, 30)) +
  scale_y_continuous(name = "Distance (cm)", breaks = seq(-60, 60, 30), labels = seq(-60, 60, 30)) +
  labs(tag = "b", size = 120) +
  theme(plot.tag = element_text(size = 16, face = "bold")) +
  simple_theme
loop_path_plot
# tiff("figures/example_paths.tiff",
#      width = 174,
#      height = 80,
#      units = "mm",
#      res = 1200)
# gridExtra::grid.arrange(straight_path_plot, loop_path_plot, ncol = 2)
# dev.off()


# Publication plots: ANOVA ------------------------------------------------
summary_df$food <- factor(summary_df$food, levels = c("buroak", "diet", "larch", "norway"))
# Probability of moving - do as a table

# Total distance
totdist_text <- data.frame(
  label = c("A", "B", "A", "A"),
  food = c("buroak", "diet", "larch", "norway"))
anova_totdist_plot <- summary_df %>% 
  filter(starve == 0) %>% 
  ggplot(aes(x = food, y = tot_dist)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("buroak" = "Bur Oak",
                              "norway" = "Nor. Maple",
                              "larch" = "Larch",
                              "diet" = "Artif. Diet")) +
  labs(x = "Food type", y = "Total distance moved (cm)", tag = "a") +
  geom_text(data = totdist_text,
            aes(x = food,
                y = c(600, 350, 600, 600),
                label = label),
            size = 4) +
  theme(plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0.18, 1.045),
        plot.margin = unit(c(0.8, 0.1, 0.1, 0.1), "cm")) +
  simple_theme
anova_totdist_plot

# Displacement
disp_text <- data.frame(
  label = c("AB", "B", "A", "AB"),
  food = c("buroak", "diet", "larch", "norway"))
anova_disp_plot <- summary_df %>% 
  filter(starve == 0) %>% 
  ggplot(aes(x = food, y = displacement)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("buroak" = "Bur Oak",
                              "norway" = "Nor. Maple",
                              "larch" = "Larch",
                              "diet" = "Artif. Diet")) +
  labs(x = "Food type", y = "Net displacement (cm)", tag = "b") +
  geom_text(data = disp_text,
            aes(x = food,
                y = c(300, 300, 550, 240),
                label = label),
            size = 4) +
  theme(plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0.18, 1.045),
        plot.margin = unit(c(0.8, 0.1, 0.1, 0.1), "cm")) +
  simple_theme

# Tortuosity
tort_text <- data.frame(
  label = c("C", "AB", "A", "BC"),
  food = c("buroak", "diet", "larch", "norway"))
anova_tort_plot <- summary_df %>% 
  filter(starve == 0) %>% 
  ggplot(aes(x = food, y = tortuosity)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("buroak" = "Bur Oak",
                              "norway" = "Nor. Maple",
                              "larch" = "Larch",
                              "diet" = "Artif. Diet")) +
  labs(x = "Food type", y = "Tortuosity", tag = "e") +
  geom_text(data = tort_text,
            aes(x = food,
                y = c(0.45, 0.58, 0.63, 0.5),
                label = label),
            size = 4) +
  theme(plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0.18, 1.045),
        plot.margin = unit(c(0.8, 0.1, 0.1, 0.1), "cm")) +
  simple_theme
anova_tort_plot

# Average velocity
avgvi_text <- data.frame(
  label = c("A", "B", "A", "A"),
  food = c("buroak", "diet", "larch", "norway")
)
anova_vi_plot <- summary_df %>% 
  filter(starve == 0) %>% 
  ggplot(aes(x = food, y = avg_vi)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("buroak" = "Bur Oak",
                              "norway" = "Nor. Maple",
                              "larch" = "Larch",
                              "diet" = "Artif. Diet")) +
  labs(x = "Food type", y = "Average velocity (cm/s)", tag = "c") +
  geom_text(data = avgvi_text,
            aes(x = food,
                y = c(1, 0.6, 1, 1),
                label = label),
            size = 4) +
  theme(plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0.18, 1.045),
        plot.margin = unit(c(0.8, 0.1, 0.1, 0.1), "cm")) +
  simple_theme
anova_vi_plot

# Number of stops
stops_text <- data.frame(
  label = c("AB", "A", "C", "B"),
  food = c("buroak", "diet", "larch", "norway"))
anova_stops_plot <- summary_df %>% 
  filter(starve == 0) %>% 
  ggplot(aes(x = food, y = stops)) +
    geom_boxplot() +
    scale_x_discrete(labels = c("buroak" = "Bur Oak",
                                "norway" = "Nor. Maple",
                                "larch" = "Larch",
                                "diet" = "Artif. Diet")) +
    scale_y_continuous(breaks = seq(0, 150, 50),
                       limits = c(0, 170)) +
    labs(x = "Food type", y = "Number of stops", tag = "d") +
    geom_text(data = stops_text,
              aes(x = food,
                  y = c(110, 160, 110, 110),
                  label = label),
              size = 4) +
  theme(plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0.18, 1.045),
        plot.margin = unit(c(0.8, 0.1, 0.1, 0.1), "cm")) +
    simple_theme
anova_stops_plot

# Grid of all plots
# tiff("figures/anova.tiff",
#      width = 174,
#      height = 234,
#      units = "mm",
#      res = 1200)
# gridExtra::grid.arrange(anova_totdist_plot,
#                         anova_disp_plot,
#                         anova_vi_plot,
#                         anova_stops_plot,
#                         anova_tort_plot,
#                         ncol = 2)
# 
# dev.off()

# Publication plots: regression ------------------------------------------
# use with facet_wrap labeller to label panels by tree
facet_labels <- c("buroak" = "Bur oak",
                  "diet" = "Diet",
                  "larch" = "Larch",
                  "norway" = "Norway maple") 

# Use with facet_wrap labeller to panel by letter
make_labelstring <- function(mypanels) {
  mylabels <- sapply(mypanels, 
                     function(x) {letters[which(mypanels == x)]})
  
  return(mylabels)
}
label_panels <- ggplot2::as_labeller(make_labelstring)

# Probability of moving
probmodel_eq1 <- substitute(logit(italic(y)) == a + b*italic(x),
                               list(a = -0.89,
                                    b = 0.07))
probmodel_eq1 <- as.character(as.expression(probmodel_eq1))
probmodel_eq2 <- substitute(logit(italic(y)) == a + b*italic(x),
                            list(a = 0.46,
                                 b = 0.03))
probmodel_eq2 <- as.character(as.expression(probmodel_eq2))
probmodel_eq3 <- substitute(logit(italic(y)) == a + b*italic(x),
                            list(a = -0.35,
                                 b = 0.10))
probmodel_eq3 <- as.character(as.expression(probmodel_eq3))
probmodel_text <- data.frame(
  label = c(probmodel_eq1,
            probmodel_eq2,
            probmodel_eq3),
  food = c("buroak","diet", "norway"))
fit_prob <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 100), 4))
fit_prob$food <- rep(c("buroak", "diet", "larch", "norway"), each = 100)
pred_prob <- predict.glm(lm_prob_move, newdata = fit_prob, type = "link", se.fit = T)
fit_prob$predicted <- plogis(pred_prob$fit)
fit_prob$low_ci <- plogis(pred_prob$fit - 1.96 * pred_prob$se.fit)
fit_prob$high_ci <- plogis(pred_prob$fit + 1.96 * pred_prob$se.fit)
fit_prob$sig <- ifelse(fit_prob$food == "larch", "no", "yes")
probmove_plot <- ggplot(data = trial_record, aes(x = starve_time, y = include)) +
  geom_point(position = position_jitter(w = 0.2, h = 0)) +
  labs(x = "Food deprivation period (hours)",
       y = "P(Moving > 10 cm)") +
  geom_line(data = fit_prob,
            aes(x = starve_time, y = predicted, linetype = sig),
            size = 0.8,
            colour = "black") +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_manual(values = c("black", "black")) +
  geom_ribbon(data = fit_prob,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) +
  facet_wrap( ~ food,
              labeller = as_labeller(label_panels)) +
  theme(strip.text = element_text(size = 16, face = "bold", hjust = -0)) +
  geom_text(data = probmodel_text,
            mapping = aes(x = 33, y = 0.30, label = label),
            parse = T) +
  simple_theme
probmove_plot
# ggsave("figures/probmove.tiff",
#        plot = probmove_plot,
#        width = 129,
#        height = 129,
#        units = "mm",
#        dpi = 1200)


# Total distance
totdistmodel_eq1 <- substitute(italic(y) == a + b*italic(x),
                            list(a = 277.81,
                                 b = 3.22))
totdistmodel_eq1 <- as.character(as.expression(totdistmodel_eq1))
totdistmodel_eq2 <- substitute(italic(y) == a + b*italic(x),
                               list(a = 79.41,
                                    b = 1.69))
totdistmodel_eq2 <- as.character(as.expression(totdistmodel_eq2))

totdistmodel_text1 <- data.frame(
  label = c(totdistmodel_eq1,
            totdistmodel_eq2),
  food = c("buroak", "diet")
)
fit_dist <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_dist$food <- rep(c("buroak", "diet", "larch", "norway"), each = 1000)
pred_dist <- predict.glm(lm_totdist_int, newdata = fit_dist, type = "link", se = T)
fit_dist$predicted <- pred_dist$fit
fit_dist$low_ci <- (fit_dist$predicted - 1.96 * pred_dist$se.fit)
fit_dist$high_ci <-(fit_dist$predicted + 1.96 * pred_dist$se.fit)
fit_dist$sig <- ifelse(fit_dist$food == "diet" | fit_dist$food == "buroak", "yes", "no")
totdist_plot <- ggplot(summary_df, aes(x = starve_time, y = tot_dist)) +
  geom_point() +
  labs(x = "Food deprivation period (hours)",
       y = "Total distance moved (cm)") +
  geom_line(data = fit_dist,
            aes(x = starve_time, y = predicted, linetype = sig),
            size = 0.8,
            colour = "black") +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_manual(values = c("black", "black")) +
  geom_ribbon(data = fit_dist,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) + 
  facet_wrap( ~ food,
              labeller = as_labeller(label_panels)) +
  theme(strip.text = element_text(size = 16, face = "bold", hjust = -0)) +
  geom_text(data = totdistmodel_text1,
            mapping = aes(x = 20, y = c(50, 450), label = label),
            parse = T) +
  simple_theme 
totdist_plot 
# ggsave("figures/totdist.tiff",
#        plot = totdist_plot,
#        width = 129,
#        height = 129,
#        units = "mm",
#        dpi = 1200)

# Displacement
dispmodel_eq1 <- substitute(italic(y) == a + b*italic(x),
                        list(a = 118.72,
                             b = 2.79))
dispmodel_eq1 <- as.character(as.expression(dispmodel_eq1))
dispmodel_eq2 <- substitute(italic(y) == a + b*italic(x),
                            list(a = 133.81,
                                 b = 2.45))
dispmodel_eq2 <- as.character(as.expression(dispmodel_eq2))
dispmodel_text1 <- data.frame(
  label = c(dispmodel_eq1,
            dispmodel_eq2),
  food = c("buroak", "norway")
)
fit_displacement <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_displacement$food <- rep(c("buroak", "diet", "larch", "norway"), each = 1000)
pred_disp <- predict.glm(lm_disp_int, newdata = fit_displacement, type = "link", se = T)
fit_displacement$predicted <- pred_disp$fit
fit_displacement$low_ci <- (fit_displacement$predicted - 1.96 * pred_disp$se.fit)
fit_displacement$high_ci <-(fit_displacement$predicted + 1.96 * pred_disp$se.fit)
fit_displacement$sig <- ifelse(fit_displacement$food == "buroak" | fit_displacement$food == "norway", "yes", "no")
disp_plot <- ggplot(summary_df, aes(x = starve_time, y = displacement)) +
  geom_point() +
  labs(x = "Food deprivation period (hours)",
       y = "Net displacement (cm)") +
  geom_line(data = fit_displacement,
            aes(x = starve_time, y = predicted, linetype = sig),
            size = 0.8) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_manual(values = c("black", "black")) +
  geom_ribbon(data = fit_displacement,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) + 
  facet_wrap( ~ food,
              labeller = as_labeller(label_panels)) + 
  theme(strip.text = element_text(size = 16, face = "bold", hjust = -0)) +
  geom_text(
    data = dispmodel_text1,
    mapping = aes(x = 25, y = 800, label = label),
    parse = T) +
  simple_theme 
disp_plot
# ggsave("figures/disp.tiff",
#        plot = disp_plot,
#        width = 129,
#        height = 129,
#        units = "mm",
#        dpi = 1200)

# Tortuosity
tortmodeq <- substitute(asin(sqrt(italic(y))) == a + b*italic(x),
                        list(a = 0.88,
                             b = 0.005))
tortmodeq <- as.character(as.expression(tortmodeq))
tortmodel_text <- data.frame(
  label = tortmodeq,
  food = "norway"
)
fit_tort <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_tort$food <- rep(c("buroak", "diet", "larch", "norway"), each = 1000)
pred_tort <- predict(lm_tort, newdata = fit_tort, type = "response", se = T)
fit_tort$predicted <- sin(pred_tort$fit) ^ 2
fit_tort$low_ci <- sin(pred_tort$fit - (1.96 * pred_tort$se.fit)) ^ 2
fit_tort$high_ci <- sin(pred_tort$fit + (1.96 * pred_tort$se.fit)) ^ 2
fit_tort$sig <- ifelse(fit_tort$food == "norway", "yes", "no")
tort_plot <- ggplot(summary_df, aes(x = starve_time, y = tortuosity)) +
  geom_point() +
  labs(x = "Food deprivation period (hours)",
       y = "Tortuosity") +
  geom_line(data = fit_tort,
            aes(x = starve_time, y = predicted, linetype = sig), 
            size = 0.8) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_manual(values = c("black", "black")) +
  geom_ribbon(data = fit_tort,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) + 
  facet_wrap( ~ food,
              labeller = as_labeller(label_panels)) +
  theme(strip.text = element_text(size = 16, face = "bold", hjust = -0)) +
  geom_text(data = tortmodel_text,
    mapping = aes(x = 25, y = 0.4, label = label),
    parse = T) +
  simple_theme 
tort_plot
# ggsave("figures/tort.tiff",
#        plot = tort_plot,
#        width = 129,
#        height = 129,
#        units = "mm",
#        dpi = 1200)




# Avg velocity
vimodel_eq1 <- substitute(italic(y) == a + b*italic(x),
                               list(a = 0.45,
                                    b = 0.005))
vimodel_eq1 <- as.character(as.expression(vimodel_eq1))
vimodel_eq2 <- substitute(italic(y) == a + b*italic(x),
                          list(a = 0.13,
                               b = 0.003))
vimodel_eq2 <- as.character(as.expression(vimodel_eq2))
vimodel_text1 <- data.frame(
  label = c(vimodel_eq1,
            vimodel_eq2),
  food = c("buroak", "diet")
)
fit_vi <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_vi$food <- rep(c("buroak", "diet", "larch", "norway"), each = 1000)
pred_vi <- predict.glm(lm_vi_int, newdata = fit_vi, type = "link", se = T)
fit_vi$predicted <- pred_vi$fit
fit_vi$low_ci <- (fit_vi$predicted - 1.96 * pred_vi$se.fit)
fit_vi$high_ci <-(fit_vi$predicted + 1.96 * pred_vi$se.fit)
fit_vi$sig <- ifelse(fit_vi$food == "buroak" | fit_vi$food == "diet", "yes", "no")
avgvi_plot <- ggplot(summary_df, aes(x = starve_time, y = avg_vi)) +
  geom_point() +
  labs(x = "Food deprivation period (hours)",
       y = "Average velocity (cm/s)") +
  geom_line(data = fit_vi,
            aes(x = starve_time, y = predicted, linetype = sig),
            size = 0.8) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_manual(values = c("black", "black")) +
  geom_ribbon(data = fit_vi,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) + 
  facet_wrap( ~ food,
              labeller = as_labeller(label_panels)) +
  theme(strip.text = element_text(size = 16, face = "bold", hjust = -0)) +
  geom_text(
    data = vimodel_text1,
    mapping = aes(x = 20, y = c(0.1, 0.9), label = label),
    parse = T
  ) +
  simple_theme 
avgvi_plot
# ggsave("figures/velocity.tiff",
#        plot = avgvi_plot,
#        width = 129,
#        height = 129,
#        units = "mm",
#        dpi = 1200)

# Number of stops
numstops_eq1 <- substitute(log(italic(y)) == a - b*italic(x),
                          list(a = 3.74,
                               b = 0.01))
numstops_eq1 <- as.character(as.expression(numstops_eq1))
numstops_eq2 <- substitute(log(italic(y)) == a - b*italic(x),
                           list(a = 4.19,
                                b = 0.01))
numstops_eq2 <- as.character(as.expression(numstops_eq2))
numstops_eq3 <- substitute(log(italic(y)) == a - b*italic(x),
                           list(a = 2.90,
                                b = 0.01))
numstops_eq3 <- as.character(as.expression(numstops_eq3))
numstops_text <- data.frame(
  label = c(numstops_eq1,
            numstops_eq2,
            numstops_eq3),
  food = c("buroak", "diet", "larch")
)
fit_stops <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_stops$food <- rep(c("buroak", "diet", "larch", "norway"), each = 1000)
pred_stops <- predict.glm(negbin_stops, newdata = fit_stops, type = "response", se = T)
fit_stops$predicted <- (pred_stops$fit)
fit_stops$low_ci <- (fit_stops$predicted - 1.96 * pred_stops$se.fit)
fit_stops$high_ci <- (fit_stops$predicted + 1.96 * pred_stops$se.fit)
fit_stops$sig <- ifelse(fit_dist$food == "norway", "no", "yes")
stops_plot <- ggplot(summary_df, aes(x = starve_time, y = stops)) +
  geom_point() +
  labs(x = "Food deprivation period (hours)",
       y = "Number of stops") +
  geom_line(data = fit_stops,
            aes(x = starve_time, y = predicted, linetype = sig),
            size = 0.8) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  geom_ribbon(data = fit_stops,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) + 
  facet_wrap( ~ food,
              labeller = as_labeller(label_panels)) +
  geom_text(
    data = numstops_text,
    mapping = aes(x = 23, y = c(130, 100, 100), label = label),
    parse = T
  ) +
  theme(strip.text = element_text(size = 16, face = "bold", hjust = -0)) +
  simple_theme 
stops_plot 
# ggsave("figures/stops.tiff",
#        plot = stops_plot,
#        width = 129,
#        height = 129,
#        units = "mm",
#        dpi = 1200)

# Presentation plots: example paths -------------------------------------------------------
# Each "fitted" dataframe has a "significance" variable to change line type
# For labelling subplots, use this labeller function

# Presentation theme
presentation_theme <- theme(panel.background=element_rect(fill=NA), 
                      panel.border = element_rect(colour = "black",
                                                  fill=NA, size=1),
                      axis.text = element_text(size = 18, colour="black"),
                      axis.title = element_text(size = 20),
                      legend.position = "none",
                      strip.background = element_rect(fill = NA))
# Movement paths
pres_straight_path_plot <- ggplot(dat = dat_files_agg[["520"]], aes(x = x, y = y)) +
  geom_path(size = 0.7) +
  coord_cartesian(xlim = c(-65, 65), ylim = c(-65, 65)) +
  scale_x_continuous(name = "Distance (cm)", breaks = seq(-60, 60, 30), labels = seq(-60, 60, 30)) +
  scale_y_continuous(name = "Distance (cm)", breaks = seq(-60, 60, 30), labels = seq(-60, 60, 30)) +
  theme(plot.tag = element_text(size = 16, face = "bold")) + 
  presentation_theme
pres_straight_path_plot

pres_loop_path_plot <- ggplot(dat_files_agg[["510"]], aes(x = x, y = y)) +
  geom_path(size = 0.7) +
  coord_cartesian(xlim = c(-65, 65), ylim = c(-65, 65)) +
  scale_x_continuous(name = "Distance (cm)", breaks = seq(-60, 60, 30), labels = seq(-60, 60, 30)) +
  scale_y_continuous(name = "Distance (cm)", breaks = seq(-60, 60, 30), labels = seq(-60, 60, 30)) +
  theme(plot.tag = element_text(size = 16, face = "bold")) +
  presentation_theme
pres_loop_path_plot
# tiff("figures/pres_example_paths.tiff",
#      width = 11.5,
#      height = 6,
#      units = "in",
#      res = 1200)
# gridExtra::grid.arrange(pres_straight_path_plot, pres_loop_path_plot, ncol = 2)
# dev.off()


# Presentation plots: ANOVA ------------------------------------------------
summary_df$food <- factor(summary_df$food, levels = c("buroak", "diet", "larch", "norway"))
# Probability of moving - do as a table

# Total distance
chisq_totdist_eq <- substitute(chi[3]^2 == a,
                            list(a = 35.140))
chisq_totdist_eq <- as.character(as.expression(chisq_totdist_eq))
p_totdist <- substitute(italic(p) < 0.001)
p_totdist <- as.character(as.expression(p_totdist))
chisq_totdist_text <- data.frame(
  label = c(chisq_totdist_eq,
            p_totdist))
totdist_text <- data.frame(
  label = c("A", "B", "A", "A"),
  food = c("buroak", "diet", "larch", "norway"))
pres_totdist <- summary_df %>% 
  filter(starve == 0) %>% 
  ggplot(aes(x = food, y = tot_dist)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("buroak" = "Bur Oak",
                              "norway" = "Norway Maple",
                              "larch" = "Larch",
                              "diet" = "Artificial Diet")) +
  labs(x = "Food type", y = "Total distance moved (cm)") +
  geom_text(data = totdist_text,
            aes(x = food,
                y = c(600, 350, 600, 600),
                label = label),
            size = 7) +
  theme(plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0.18, 1.045),
        plot.margin = unit(c(0.8, 0.1, 0.1, 0.1), "cm")) +
  geom_text(data = chisq_totdist_text,
            mapping = aes(x = 2, y = c(500, 440), label = label),
            size =6,
            parse = T) +
  presentation_theme
pres_totdist
# ggsave("figures/pres_totdist_anova.tiff",
#        plot = pres_totdist,
#        width = 9,
#        height = 6,
#        units = "in",
#        dpi = 1200)
       

# Displacement
disp_text <- data.frame(
  label = c("AB", "B", "A", "AB"),
  food = c("buroak", "diet", "larch", "norway"))
pres_disp <- summary_df %>% 
  filter(starve == 0) %>% 
  ggplot(aes(x = food, y = displacement)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("buroak" = "Bur Oak",
                              "norway" = "Nor. Maple",
                              "larch" = "Larch",
                              "diet" = "Artif. Diet")) +
  labs(x = "Food type", y = "Net displacement (cm)", tag = "b") +
  geom_text(data = disp_text,
            aes(x = food,
                y = c(300, 300, 550, 240),
                label = label),
            size = 4) +
  theme(plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0.18, 1.045),
        plot.margin = unit(c(0.8, 0.1, 0.1, 0.1), "cm")) +
  presentation_theme

# Tortuosity
tort_text <- data.frame(
  label = c("C", "AB", "A", "BC"),
  food = c("buroak", "diet", "larch", "norway"))
pres_tort <- summary_df %>% 
  filter(starve == 0) %>% 
  ggplot(aes(x = food, y = tortuosity)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("buroak" = "Bur Oak",
                              "norway" = "Nor. Maple",
                              "larch" = "Larch",
                              "diet" = "Artif. Diet")) +
  labs(x = "Food type", y = "Tortuosity", tag = "e") +
  geom_text(data = tort_text,
            aes(x = food,
                y = c(0.45, 0.58, 0.63, 0.5),
                label = label),
            size = 4) +
  theme(plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0.18, 1.045),
        plot.margin = unit(c(0.8, 0.1, 0.1, 0.1), "cm")) +
  presentation_theme
pres_tort

# Average velocity
chisq_vi_eq <- substitute(chi[3]^2 == a,
                               list(a = 35.263))
chisq_vi_eq <- as.character(as.expression(chisq_vi_eq))
p_vi <- substitute(italic(p) < 0.001)
p_vi <- as.character(as.expression(p_vi))
chisq_vi_text <- data.frame(
  label = c(chisq_vi_eq,
            p_vi))
avgvi_text <- data.frame(
  label = c("A", "B", "A", "A"),
  food = c("buroak", "diet", "larch", "norway")
)
pres_vi <- summary_df %>% 
  filter(starve == 0) %>% 
  ggplot(aes(x = food, y = avg_vi)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("buroak" = "Bur Oak",
                              "norway" = "Norway Maple",
                              "larch" = "Larch",
                              "diet" = "Artificial Diet")) +
  labs(x = "Food type", y = "Average velocity (cm/s)") +
  geom_text(data = avgvi_text,
            aes(x = food,
                y = c(1, 0.6, 1, 1),
                label = label),
            size = 7) +
  theme(plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0.18, 1.045),
        plot.margin = unit(c(0.8, 0.1, 0.1, 0.1), "cm")) +
  geom_text(data = chisq_vi_text,
            mapping = aes(x = c(2, 2), y = c(0.9, 0.8), label = label),
            size =6,
            parse = T) +
  presentation_theme
pres_vi
# ggsave("figures/pres_vi_anova.tiff",
#        plot = pres_vi,
#        width = 9,
#        height = 6,
#        units = "in",
#        dpi = 1200)

# Number of stops
chisq_stops_eq <- substitute(chi[3]^2 == a,
                          list(a = 38.810))
chisq_stops_eq <- as.character(as.expression(chisq_stops_eq))
p_stops <- substitute(italic(p) < 0.001)
p_stops <- as.character(as.expression(p_stops))
chisq_stops_text <- data.frame(
  label = c(chisq_stops_eq,
            p_stops))
stops_text <- data.frame(
  label = c("AB", "A", "C", "B"),
  food = c("buroak", "diet", "larch", "norway"))
pres_stops <- summary_df %>% 
  filter(starve == 0) %>% 
  ggplot(aes(x = food, y = stops)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("buroak" = "Bur Oak",
                              "norway" = "Norway Maple",
                              "larch" = "Larch",
                              "diet" = "Artificial Diet")) +
  scale_y_continuous(breaks = seq(0, 150, 50),
                     limits = c(0, 170)) +
  labs(x = "Food type", y = "Number of stops") +
  geom_text(data = stops_text,
            aes(x = food,
                y = c(110, 160, 110, 110),
                label = label),
            size = 7) +
  theme(plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0.18, 1.045),
        plot.margin = unit(c(0.8, 0.1, 0.1, 0.1), "cm")) +
  geom_text(data = chisq_vi_text,
            mapping = aes(x = c(3, 3), y = c(145, 130), label = label),
            size =6,
            parse = T) +
  presentation_theme
pres_stops

# ggsave("figures/pres_stops_anova.tiff",
#        plot = pres_stops,
#        width = 9,
#        height = 6,
#        units = "in",
#        dpi = 1200)

# Grid of all plots
# tiff("figures/anova.tiff",
#      width = 174,
#      height = 234,
#      units = "mm",
#      res = 1200)
# gridExtra::grid.arrange(pres_totdist,
#                         pres_disp,
#                         pres_vi,
#                         pres_stops,
#                         pres_tort,
#                         ncol = 2)
# 
# dev.off()

# Presentation plots: regression ------------------------------------------
# use with facet_wrap labeller to label panels by tree
facet_labels <- c("buroak" = "Bur oak",
                  "diet" = "Diet",
                  "larch" = "Larch",
                  "norway" = "Norway maple") 


# Probability of moving
probmodel_eq1 <- substitute(logit(italic(y)) == a + b*italic(x),
                            list(a = -0.89,
                                 b = 0.07))
probmodel_eq1 <- as.character(as.expression(probmodel_eq1))
probmodel_eq2 <- substitute(logit(italic(y)) == a + b*italic(x),
                            list(a = 0.46,
                                 b = 0.03))
probmodel_eq2 <- as.character(as.expression(probmodel_eq2))
probmodel_eq3 <- substitute(logit(italic(y)) == a + b*italic(x),
                            list(a = -0.35,
                                 b = 0.10))
probmodel_eq3 <- as.character(as.expression(probmodel_eq3))
probmodel_text <- data.frame(
  label = c(probmodel_eq1,
            probmodel_eq2,
            probmodel_eq3),
  food = c("buroak","diet", "norway"))
fit_prob <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 100), 4))
fit_prob$food <- rep(c("buroak", "diet", "larch", "norway"), each = 100)
pred_prob <- predict.glm(lm_prob_move, newdata = fit_prob, type = "link", se.fit = T)
fit_prob$predicted <- plogis(pred_prob$fit)
fit_prob$low_ci <- plogis(pred_prob$fit - 1.96 * pred_prob$se.fit)
fit_prob$high_ci <- plogis(pred_prob$fit + 1.96 * pred_prob$se.fit)
fit_prob$sig <- ifelse(fit_prob$food == "larch", "no", "yes")
pres_prob_reg <- ggplot(data = trial_record, aes(x = starve_time, y = include)) +
  geom_point(position = position_jitter(w = 0.2, h = 0)) +
  labs(x = "Food deprivation period (hours)",
       y = "P(Moving > 10 cm)") +
  geom_line(data = fit_prob,
            aes(x = starve_time, y = predicted, linetype = sig),
            size = 0.8,
            colour = "black") +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_manual(values = c("black", "black")) +
  geom_ribbon(data = fit_prob,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) +
  facet_wrap( ~ food,
              labeller = as_labeller(facet_labels)) +
  theme(strip.text = element_text(size = 16, face = "bold", hjust = -0)) +
  geom_text(data = probmodel_text,
            mapping = aes(x = 33, y = 0.30, label = label),
            size = 6,
            parse = T) +
  presentation_theme
pres_prob_reg
# ggsave("figures/pres_probmove_reg.tiff",
#        plot = pres_prob_reg,
#        width = 9,
#        height = 6.3,
#        units = "in",
#        dpi = 1200)


# Total distance
totdistmodel_eq1 <- substitute(italic(y) == a + b*italic(x),
                               list(a = 277.81,
                                    b = 3.22))
totdistmodel_eq1 <- as.character(as.expression(totdistmodel_eq1))
totdistmodel_eq2 <- substitute(italic(y) == a + b*italic(x),
                               list(a = 79.41,
                                    b = 1.69))
totdistmodel_eq2 <- as.character(as.expression(totdistmodel_eq2))

totdistmodel_text1 <- data.frame(
  label = c(totdistmodel_eq1,
            totdistmodel_eq2),
  food = c("buroak", "diet")
)
fit_dist <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_dist$food <- rep(c("buroak", "diet", "larch", "norway"), each = 1000)
pred_dist <- predict.glm(lm_totdist_int, newdata = fit_dist, type = "link", se = T)
fit_dist$predicted <- pred_dist$fit
fit_dist$low_ci <- (fit_dist$predicted - 1.96 * pred_dist$se.fit)
fit_dist$high_ci <-(fit_dist$predicted + 1.96 * pred_dist$se.fit)
fit_dist$sig <- ifelse(fit_dist$food == "diet" | fit_dist$food == "buroak", "yes", "no")
pres_totdist_reg <- ggplot(summary_df, aes(x = starve_time, y = tot_dist)) +
  geom_point() +
  labs(x = "Food deprivation period (hours)",
       y = "Total distance moved (cm)") +
  geom_line(data = fit_dist,
            aes(x = starve_time, y = predicted, linetype = sig),
            size = 0.8,
            colour = "black") +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_manual(values = c("black", "black")) +
  geom_ribbon(data = fit_dist,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) + 
  facet_wrap( ~ food,
              labeller = as_labeller(facet_labels)) +
  theme(strip.text = element_text(size = 16, face = "bold", hjust = -0)) +
  geom_text(data = totdistmodel_text1,
            mapping = aes(x = 20, y = c(50, 450), label = label),
            size = 6,
            parse = T) +
  presentation_theme 
pres_totdist_reg 

# ggsave("figures/pres_totdist_reg.tiff",
#        plot = pres_totdist_reg,
#        width = 9,
#        height = 6.3,
#        units = "in",
#        dpi = 1200)

# Displacement
dispmodel_eq1 <- substitute(italic(y) == a + b*italic(x),
                            list(a = 118.72,
                                 b = 2.79))
dispmodel_eq1 <- as.character(as.expression(dispmodel_eq1))
dispmodel_eq2 <- substitute(italic(y) == a + b*italic(x),
                            list(a = 133.81,
                                 b = 2.45))
dispmodel_eq2 <- as.character(as.expression(dispmodel_eq2))
dispmodel_text1 <- data.frame(
  label = c(dispmodel_eq1,
            dispmodel_eq2),
  food = c("buroak", "norway")
)
fit_displacement <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_displacement$food <- rep(c("buroak", "diet", "larch", "norway"), each = 1000)
pred_disp <- predict.glm(lm_disp_int, newdata = fit_displacement, type = "link", se = T)
fit_displacement$predicted <- pred_disp$fit
fit_displacement$low_ci <- (fit_displacement$predicted - 1.96 * pred_disp$se.fit)
fit_displacement$high_ci <-(fit_displacement$predicted + 1.96 * pred_disp$se.fit)
fit_displacement$sig <- ifelse(fit_displacement$food == "buroak" | fit_displacement$food == "norway", "yes", "no")
pres_disp_reg <- ggplot(summary_df, aes(x = starve_time, y = displacement)) +
  geom_point() +
  labs(x = "Food deprivation period (hours)",
       y = "Net displacement (cm)") +
  geom_line(data = fit_displacement,
            aes(x = starve_time, y = predicted, linetype = sig),
            size = 0.8) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_manual(values = c("black", "black")) +
  geom_ribbon(data = fit_displacement,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) + 
  facet_wrap( ~ food,
              labeller = as_labeller(label_panels)) + 
  theme(strip.text = element_text(size = 16, face = "bold", hjust = -0)) +
  geom_text(
    data = dispmodel_text1,
    mapping = aes(x = 25, y = 800, label = label),
    parse = T) +
  presentation_theme 
pres_disp_reg
# ggsave("figures/pres_disp_reg.tiff",
#        plot = pres_disp_reg,
#        width = 129,
#        height = 129,
#        units = "mm",
#        dpi = 1200)

# Tortuosity
tortmodeq <- substitute(asin(sqrt(italic(y))) == a + b*italic(x),
                        list(a = 0.88,
                             b = 0.005))
tortmodeq <- as.character(as.expression(tortmodeq))
tortmodel_text <- data.frame(
  label = tortmodeq,
  food = "norway"
)
fit_tort <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_tort$food <- rep(c("buroak", "diet", "larch", "norway"), each = 1000)
pred_tort <- predict(lm_tort, newdata = fit_tort, type = "response", se = T)
fit_tort$predicted <- sin(pred_tort$fit) ^ 2
fit_tort$low_ci <- sin(pred_tort$fit - (1.96 * pred_tort$se.fit)) ^ 2
fit_tort$high_ci <- sin(pred_tort$fit + (1.96 * pred_tort$se.fit)) ^ 2
fit_tort$sig <- ifelse(fit_tort$food == "norway", "yes", "no")
pres_tort_reg <- ggplot(summary_df, aes(x = starve_time, y = tortuosity)) +
  geom_point() +
  labs(x = "Food deprivation period (hours)",
       y = "Tortuosity") +
  geom_line(data = fit_tort,
            aes(x = starve_time, y = predicted, linetype = sig), 
            size = 0.8) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_manual(values = c("black", "black")) +
  geom_ribbon(data = fit_tort,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) + 
  facet_wrap( ~ food,
              labeller = as_labeller(facet_labels)) +
  theme(strip.text = element_text(size = 16, face = "bold", hjust = -0)) +
  geom_text(data = tortmodel_text,
            mapping = aes(x = 25, y = 0.4, label = label),
            size = 6,
            parse = T) +
  presentation_theme 
pres_tort_reg
# ggsave("figures/pres_tort_reg.tiff",
#        plot = pres_tort_reg,
#        width = 9,
#        height = 6.3,
#        units = "in",
#        dpi = 1200)




# Avg velocity
vimodel_eq1 <- substitute(italic(y) == a + b*italic(x),
                          list(a = 0.45,
                               b = 0.005))
vimodel_eq1 <- as.character(as.expression(vimodel_eq1))
vimodel_eq2 <- substitute(italic(y) == a + b*italic(x),
                          list(a = 0.13,
                               b = 0.003))
vimodel_eq2 <- as.character(as.expression(vimodel_eq2))
vimodel_text1 <- data.frame(
  label = c(vimodel_eq1,
            vimodel_eq2),
  food = c("buroak", "diet")
)
fit_vi <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_vi$food <- rep(c("buroak", "diet", "larch", "norway"), each = 1000)
pred_vi <- predict.glm(lm_vi_int, newdata = fit_vi, type = "link", se = T)
fit_vi$predicted <- pred_vi$fit
fit_vi$low_ci <- (fit_vi$predicted - 1.96 * pred_vi$se.fit)
fit_vi$high_ci <-(fit_vi$predicted + 1.96 * pred_vi$se.fit)
fit_vi$sig <- ifelse(fit_vi$food == "buroak" | fit_vi$food == "diet", "yes", "no")
pres_vi_reg <- ggplot(summary_df, aes(x = starve_time, y = avg_vi)) +
  geom_point() +
  labs(x = "Food deprivation period (hours)",
       y = "Average velocity (cm/s)") +
  geom_line(data = fit_vi,
            aes(x = starve_time, y = predicted, linetype = sig),
            size = 0.8) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_colour_manual(values = c("black", "black")) +
  geom_ribbon(data = fit_vi,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) + 
  facet_wrap( ~ food,
              labeller = as_labeller(facet_labels)) +
  theme(strip.text = element_text(size = 16, face = "bold", hjust = -0)) +
  geom_text(
    data = vimodel_text1,
    mapping = aes(x = 20, y = c(0.1, 0.9), label = label),
    size = 6,
    parse = T
  ) +
  presentation_theme 
pres_vi_reg
# ggsave("figures/pres_velocity_reg.tiff",
#        plot = pres_vi_reg,
#        width = 9,
#        height = 6.3,
#        units = "in",
#        dpi = 1200)

# Number of stops
numstops_eq1 <- substitute(log(italic(y)) == a - b*italic(x),
                           list(a = 3.74,
                                b = 0.01))
numstops_eq1 <- as.character(as.expression(numstops_eq1))
numstops_eq2 <- substitute(log(italic(y)) == a - b*italic(x),
                           list(a = 4.19,
                                b = 0.01))
numstops_eq2 <- as.character(as.expression(numstops_eq2))
numstops_eq3 <- substitute(log(italic(y)) == a - b*italic(x),
                           list(a = 2.90,
                                b = 0.01))
numstops_eq3 <- as.character(as.expression(numstops_eq3))
numstops_text <- data.frame(
  label = c(numstops_eq1,
            numstops_eq2,
            numstops_eq3),
  food = c("buroak", "diet", "larch")
)
fit_stops <- data.frame(starve_time = rep(seq(0, max(summary_df$starve_time), length.out = 1000), 4))
fit_stops$food <- rep(c("buroak", "diet", "larch", "norway"), each = 1000)
pred_stops <- predict.glm(negbin_stops, newdata = fit_stops, type = "response", se = T)
fit_stops$predicted <- (pred_stops$fit)
fit_stops$low_ci <- (fit_stops$predicted - 1.96 * pred_stops$se.fit)
fit_stops$high_ci <- (fit_stops$predicted + 1.96 * pred_stops$se.fit)
fit_stops$sig <- ifelse(fit_dist$food == "norway", "no", "yes")
pres_stops_reg <- ggplot(summary_df, aes(x = starve_time, y = stops)) +
  geom_point() +
  labs(x = "Food deprivation period (hours)",
       y = "Number of stops") +
  geom_line(data = fit_stops,
            aes(x = starve_time, y = predicted, linetype = sig),
            size = 0.8) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  geom_ribbon(data = fit_stops,
              aes(x = starve_time, ymin = low_ci, ymax = high_ci),
              inherit.aes = F,
              fill = "grey",
              alpha = 0.4) + 
  facet_wrap( ~ food,
              labeller = as_labeller(facet_labels)) +
  geom_text(
    data = numstops_text,
    mapping = aes(x = 23, y = c(130, 100, 100), label = label),
    size = 6,
    parse = T
  ) +
  theme(strip.text = element_text(size = 16, face = "bold", hjust = -0)) +
  presentation_theme 
pres_stops_reg 
# ggsave("figures/pres_stops_reg.tiff",
#        plot = pres_stops_reg,
#        width = 9,
#        height = 6.3,
#        units = "in",
#        dpi = 1200)

# Example of a tidyeval function for plotting
testPlotFunction <- function(df, xvar, yvar, facet){
  xvar <- enexpr(xvar)
  yvar <- enexpr(yvar)
  facet <- enquo(facet)
  ggplot(df, aes(!!xvar, !!yvar)) +
    geom_point() +
    facet_wrap(facet,
               labeller = as_labeller(facet_labels)) +
    presentation_theme
}
testPlotFunction(summary_df, starve_time, stops, food)
