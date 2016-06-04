#Particle Swarm Optimization for SVR
setwd("~/Backrest Upholstery Durability Development/Publishing/Git")
data <- read.csv("Fabric1 ALT.csv")
library(kernlab)
library(plyr)
library(dplyr)

t9 <- Sys.time()
#Define Training and Test Set
training.data <- data[,c(1:4)] #Removing Gender Classifications of Male/Female
data.count <- count(training.data) #Count Observations
training.count <- as.integer(round(0.85*data.count)) #Specifying training count at 85%
test.count <- round(data.count - training.count) #Test Count at 15%
set.seed(training.count) #Set seed to ensure repeatability
index <- 1:nrow(training.data)
training.sample <- sample(index,training.count) #Sampling 85% of Data
train <- training.data[training.sample,]#Separating Training Data
test <- training.data[-training.sample,] #separating/Indentifying Test Data 15% of Data
training.sort <- arrange(train, Cycles) #Ordering Training Data by Cycles
test.sort <- arrange(test, Cycles) #Ordering Test Data by Cycles

#Define Validation Set
val.count <- as.integer(data.count - training.count) #Validation count
val.sample <- sample(index, val.count) #sampling validation set 15% of data
val <- training.data[val.sample,]
val.sort <- arrange(val, Cycles) #Sorting validation set by Cycles
rmse <- function(error){
      sqrt(mean(error^2))
      }

#Define number of particles set to 30
particle <- 30
Chi <- 0.72984 #Constriction Factor

#Define Particle's Search Bounds
set.seed(231365384)
n <- 5 # Number of iterations
#cost
C.min <- 0.1
C.max <- 2000
vel.C.max <- 0.1*(C.max-C.min)
C <- runif(particle, C.min,C.max)
#epsilon
eps.min <- 0.00085358
eps.max <- 0.25
vel.eps.max <- 0.1*(eps.max-eps.min)
eps <- runif(particle, eps.min,eps.max)
#gamma
gam.min <- 0.000001
gam.max <- 50
vel.gam.max <- 0.1*(gam.max-gam.min)
gam <- runif(particle, gam.min,gam.max)
library(kernlab)
#sigma estimation
srange <- sigest(Tensile.Strength ~ ., training.sort)
s <- runif(particle, srange[1], srange[3])
vel.s.max <- 0.1*(srange[3] - srange[1])

#Initialize Particle's Position
s.i.j <- mapply(c, C, eps, gam, s)
pbest <- s.i.j
#Define Particle's Neighborhood
gbest <- s.i.j

#Initialize Particle's Velocities
v.i.j.init <- c(vel.C.max, vel.eps.max, vel.gam.max, vel.s.max)
v.i.j <- matrix(v.i.j.init,4,particle)
#Initialize Particle Best Matrix for comparison
value <- 1
RMSE.best <- matrix(value, 1, particle)
t10 <- Sys.time()
Data.Loading.Time <- difftime(t10,t9)
#Run Particle Swarm Optimization
set.seed(13846598)
result = vector("list")
#rbf <-rbfdot(sigma = srange[2])
#train.matrix <- data.matrix(training.sort)
#test.matrix <- data.matrix(test.sort)
#Kern.Matrix <- kernelMatrix(rbf, train.matrix)
#Kern.Matrix.Test <- kernelMatrix(rbf, test.matrix)
t1 <- Sys.time()
for(i in 1:n){
      #SVR RMSE Model
      t3 <- Sys.time()
      for(j in 1:particle){
            #Train SVM
            t11 <- Sys.time()
            model <- ksvm(Tensile.Strength ~ ., training.sort, type = "eps-svr", kernel = "rbfdot", 
                        kpar = list(sigma = s.i.j[4,j]), C = s.i.j[1,j], epsilon = s.i.j[2,j], gamma = s.i.j[3,j], cross=0)
            t12 <- Sys.time()
            SVR.Training.Time <- difftime(t12, t11)
            #Predict SVM
            predictedY <- predict(model, test.sort)
            t13 <- Sys.time()
            SVR.Prediction.Time <- difftime(t13, t12)
            #Calculate NRMSE
            error <- test.sort$Tensile.Strength - predictedY
            particleRMSE <- rmse(error)
            #Store iterative Results of RMSE
            result[[j]] <- particleRMSE
            rmse.matrix <- matrix(unlist(result), ncol = particle, byrow = FALSE)
      }
      t4 <- Sys.time()
      SVR.Particle.Train.Time <- difftime(t4, t3)
      if(particleRMSE < RMSE.best[,j]){
            RMSE.best[,j] = particleRMSE
            #Update Particle's Best Position pBest
            pbest[,j] <- s.i.j[,j]
      }
      #Update Global Best and particle's neighbors
      glob.best <- rbind(pbest, rmse.matrix)      
      for(k in 1:particle){
            if (k - 1 == 0){
                  if (glob.best[5,k] < glob.best[5,k+1] & glob.best[5,particle]){
                        gbest[1:4,k] = glob.best[1:4,k]      
                  }  else if (glob.best[5,k+1] < glob.best[5,particle]){
                        gbest[1:34,k] = glob.best[1:4,k+1]
                  }  else {
                        gbest[1:4,k] = glob.best[1:4,particle]
                  }
            } else if (k + 1 == particle + 1){
                  if(glob.best[5,k] < glob.best[5,1] & glob.best[5,k-1]){
                        gbest[1:4,k] = glob.best[1:4,k]
                  }  else if (glob.best[5,1] < glob.best[5,k-1]){
                        gbest[1:4,k] = glob.best[1:4,1]
                  }  else {
                        gbest[1:4,k] = glob.best[1:4,k-1]
                  }  
            } else {
                  if(glob.best[5,k] < glob.best[5,k+1] & glob.best[5,k-1]){
                        gbest[1:4,k] = glob.best[1:4,k]
                  }  else if (glob.best[5,k+1] < glob.best[5,k-1]){
                        gbest[1:4,k] = glob.best[1:4,k + 1] 
                  }  else {
                        gbest[1:4,k] = glob.best[1:4,k - 1]
                  }  
            }
      }
      #Update Particle's Velocities and Positions
      t7 <- Sys.time()
      u1 <- runif(1)
      u2 <- runif(1)
      c1 <- 1.49 #Acceleration Constant
      c2 <- 1.49 #Acceleration Constant
      v.i.j <- Chi * (v.i.j +c1*u1*(pbest - s.i.j) +c2*u2*(gbest - s.i.j))
      s.i.j <- s.i.j + v.i.j
      remove(glob.best)
      #Is Particle Feasible?
      for(l in 1:particle){
            if(s.i.j[1,l] < C.min){
                  s.i.j[1,l] = C.min
            }
            if(s.i.j[1,l] > C.max){
                  s.i.j[1,l] = C.max
            }
            if(s.i.j[2,l] < eps.min){
                  s.i.j[2,l] = eps.min
            }
            if(s.i.j[2,l] > eps.max){
                  s.i.j[2,l] = eps.max
            }
            if(s.i.j[3,l] < gam.min){
                  s.i.j[3,l] = gam.min
            }
            if(s.i.j[3,l] > gam.max){
                  s.i.j[3,l] = gam.max
            }
            if(s.i.j[4,l] < srange[1]){
                  s.i.j[4,l] = srange[1]
            }
            if(s.i.j[4,l] > srange[3]){
                  s.i.j[4,l] = srange[3]
            }
      }
      #Stop Criterion Met?
      
}
t2 <- Sys.time()
PSO.Iteration.Time <- difftime(t2, t1)
SVR.Training.Time
SVR.Particle.Train.Time
SVR.Prediction.Time
PSO.Iteration.Time
Data.Loading.Time
#BEST SVM
#best.model <- ksvm(Tensile.Strength ~ . , training.sort, type = "eps-svr", kernel = "rbfdot", 
                   #kpar = list(sigma = s.i.j[4,1]), C = s.i.j[1,1], epsilon = s.i.j[2,1], gamma = s.i.j[3,1], cross=0)
best.model <- ksvm(Tensile.Strength ~ . , training.sort, type = "eps-svr", kernel = "rbfdot", 
                   kpar = list(sigma = 0.68886841), C = 554.25474130, epsilon = 0.06986307, gamma = 17.31763199, cross=0)
#Predicting the best model
best.predictedY <- predict(best.model, val.sort)
error <- val.sort$Tensile.Strength - best.predictedY
best.RMSE <- rmse.matrix[,1]
plot(training.sort$Cycles, training.sort$Tensile.Strength, pch = 16,
            xlab = "Cycles", ylab = "Remaining Tensile Strength")
      lines(val.sort$Cycles, best.predictedY, col = "red", pch = 4)
      legend("topright",c('ALT Data', 'Best SVR Model'),
            col=c(1,2),
            pch=c(16,4), 
            text.col=c(1,2))
      title(main="ALT SVR Model")

#Calculate Best NRMSE
#s.i.j[,1]
best.RMSE
#Bootstrapped Confidence Intervals
#Re-Center Residuals

set.seed(9000000) #set seed to ensure repeatability
mean.Res.Center <- best.model@fitted/training.count
Res.Center <- best.model@fitted - mean.Res.Center
B1 <- 99
B2 <- 499
B3 <- 999
Boot.list1 <- list()
Boot.list2 <- list()
Boot.list3 <- list()

#Bootstrap "for" loop.
for(b in 1:B3){
      #Generate Rademacher Variables
      rademacher <- sample(c(-1,1), training.count, replace = TRUE, prob = c(0.5,0.5))
      #Calculate
      Rad.Res <- rademacher * Res.Center
      #Generate i.b
      
      #Set variable range update training.sort
      D <- training.sort$Tensile.Strength + Rad.Res
      Boot <- cbind(training.sort[,1:3], D)
      colnames(Boot) <- c("Cycles", "Process", "Weight", "Tensile.Strength")
      #Train SVR
      Boot.SVR <- ksvm(Tensile.Strength ~ ., Boot, type = "eps-svr", kernel = "rbfdot", 
                       kpar = list(sigma = s.i.j[4,1]),
                       cost = s.i.j[1,1], epsilon = s.i.j[2,1], gamma = s.i.j[3,1])
      Boot.pred <- predict(Boot.SVR, test.sort)
      Boot.list1[[b]] <- Boot.pred
}
df <- do.call("rbind", Boot.list1)
Boot.Data <- stack(as.data.frame(Boot.list1))
Boot.DF <- t(mapply(c, test.sort$Cycles, t(Boot.Data$values)))
colnames(Boot.DF) <- c("Cycles", "Tensile.Strength")
Bag <- rowSums(as.data.frame(Boot.list1))
ybag <- t(mapply(c, test.sort$Cycles, (test.sort$Tensile.Strength + Bag)/(B3 + 1)))
colnames(ybag) <- c("Cycles", "Tensile.Strength")
CI.SVR <- ksvm(Tensile.Strength ~ ., ybag, type = "eps-svr", kernel = "rbfdot",
               kpar = list(sigma = s.i.j[4,1]),
               cost = s.i.j[1,1], epsilon = s.i.j[2,1], gamma = s.i.j[3,1])
set.seed(500000)
#bootstrapping the validation set.
Boot.Test <- t(mapply(c, val.sort$Cycles, val.sort$Tensile.Strength))
colnames(Boot.Test) <- c("Cycles","Tensile.Strength")
CI <- predict(CI.SVR, val.sort)
CI.Line <- data.frame(t(mapply(c, val.sort$Cycles, CI)))
colnames(CI.Line) <- c("Cycles","Tensile.Strength")
CI.sort <- arrange(CI.Line, Cycles)

#Viewing Predicted Model by Sampling Cycles
set.seed(500000)
Cycles <- sample(0:120000, 5000, replace=T)
#Process <- matrix(2, nrow = 1, ncol = 5000)
Process <- sample(1:3, 5000, replace=T)
#Weight <- runif(5000, 0.49, 1.5)
Weight <- sample(c(0.598,0.651,0.492,0.827,0.862,.9375,1.3125), 5000, replace = TRUE)
TTF.sampling <- t(mapply(c, Cycles, Process, Weight))
colnames(TTF.sampling) <- c("Cycles", "Process", "Weight")
dist.TTF <- predict(best.model, TTF.sampling)
TTF <- data.frame(t(mapply(c, Cycles, dist.TTF)))
colnames(TTF) <- c("Cycles", "Tensile.Strength")
#pairs(~Cycles+Tensile.Strength+Weight+Process,data=TTF,
      #main="Simple Scatterplot Matrix")
TTF.sort <- arrange(TTF, Cycles)
#Preciting the Confidence inter vals.
CI.Pred <- predict(CI.SVR, Cycles)
CI.TTF <- data.frame(t(mapply(c, Cycles, CI.Pred)))
colnames(CI.TTF) <- c("Cycles", "Tensile.Strength")
CI.TTF.sort <- arrange(CI.TTF, Cycles)

#Plotting the data
pdf(file="plot.pdf",family="Times", pointsize=16, width=16,height=10)
plot(training.sort$Cycles, training.sort$Tensile.Strength, pch = 1,
     xlab = "Cycles", ylab = "Remaining Tensile Strength", xlim=c(0, 120000), ylim = c(0,11))
      #lines(val.sort$Cycles, best.predictedY, col = "red", pch = 4)
      #points(TTF.sort$Cycles, TTF.sort$Tensile.Strength, col = "purple", pch = 4, cex = 0.3)
      #lines(, lty=2, col="purple")
      #lines(CI.sort$Cycles, CI.sort$Tensile.Strength, col = "green", pch = 4)
      lines(CI.TTF.sort, col = "black", lty = 1, cex = 0.4)
      #points(ybag, col = "green", pch = 4, cex = 0.6)
      #legend("topright",c('Combined ALT and Field Data'
                          #, 'Predicted Values'
                          #, 'Boostrapped SVR Model'
                          #, 'ybag points'
                          #, 'Warranty Time Period'),
           # col=c(1
                  #,"purple"
                  #, "blue"
                  #, "green"
                  #, "red"),
           # pch=c(16
                  #,4
                  #,4,4), 
           # text.col=c(1
                       #,"purple"
                       #, "blue"
                       #, "green"
                       #, "red"))
      #title(main="Suspension Fabric Raw Data")
      #abline(v = 73618, col = "black", lty = 5)

library(ggvis)
library(shiny)
training.data %>%
      ggvis(~Cycles, ~Tensile.Strength,
            fill := "blue", stroke := "black") %>%
      Process = input_select(
            c("Process 1" = "1",
              "Process 2" = "2",
              "Process 3" = "3"),
            label = "Process") %>%
      layer_points(fill = ~factor(Process)) %>%
      scale_numeric("x", domain = c(0, 120000), nice = FALSE) %>%
      scale_numeric("y", domain = c(0, 10), nice = FALSE)