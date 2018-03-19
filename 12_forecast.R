require(MASS)

#### Actually, I should be doing each year, but I don't have time

head(plot_df)
plot_df$trigger <- factor(plot_df$trigger, c("None", "Moderate", "Severe", "Extreme"))

### Number of samples
n_samples <- dim(plot_df)[1]
n_train <- round(n_samples*0.9)

### Random sample (10-fold, 90% for training, 10% for testing)
train <- sample(1:n_samples, n_train)
table(plot_df$trigger[train])

### Priors based on trigger values * mean
priors <- c(0,0.5*0.7 - 0.5*0.5,0.5*0.5 - 0.5*0.25,0.5*0.25)
priors[1] <- 1- sum(priors[2:4])
priors

### Linear Discriminant Function Analysis
#z <- lda(trigger ~ ., plot_df, prior = priors, subset = train)
z <- lda(trigger ~ min_perc, plot_df, prior = priors, subset = train)

(z1 <- update(z, . ~ . + dura_months))


### Check results
predict(z, plot_df[-train, ])$class
plot_df[-train, ]$trigger
table(predict(z, plot_df[-train, ])$class, plot_df[-train, ]$trigger)
##  [1] s s s s s s s s s s s s s s s s s s s s s s s s s s s c c c
## [31] c c c c c c c v c c c c v c c c c c c c c c c c c v v v v v
## [61] v v v v v v v v v v v v v v v






require(caret)


inTraining <- createDataPartition(plot_df$trigger, p = .75, list = FALSE)
training <- plot_df[ inTraining,]
testing  <- plot_df[-inTraining,]






to_train <- plot_df[,c(seq(3,16))]
to_train$trigger <- plot_df$trigger
to_train <- to_train[complete.cases(to_train),]

inTraining <- createDataPartition(to_train$trigger, p = .75, list = FALSE)
training <- to_train[ inTraining,]
testing  <- to_train[-inTraining,]


### Parameter Tuning
fitControl <- trainControl(## 10-fold CV
                           method = "repeatedcv",
                           number = 10,
                           ## repeated ten times
                           repeats = 10)
                           

   
Fit1 <- train( trigger ~ min_perc, data = to_train, 
                 method = "lda", 
                 trControl = fitControl)
Fit1
          
     
          
Fit2 <- train( trigger ~ ., data = training, 
                 method = "stepLDA", maxvar = 5,
                 trControl = fitControl)
                 
 
 
 Fit2 <- train( x=to_train[, !names(to_train) %in% c("trigger")], y=to_train$trigger, 
                 method = "stepLDA", maxvar = 5,
                 trControl = fitControl)
                 
                 
                 
                                 
                 , maxvar=5, direction="both")



sc_obj <- stepclass(trigger ~ ., data= to_train, "lda", start.vars = "dura_months")
sc_obj
plot(sc_obj)


sc_obj <- stepclass(trigger ~ ., data= to_train, method = "qda", 
    start.vars = "dura_months", maxvars=5, criterion = "AS")  # same as above 
sc_obj



                           
fds <- train(trigger ~ ., data = training,
              method = "stepLDA",
              trControl = trainControl(method = "cv"))


