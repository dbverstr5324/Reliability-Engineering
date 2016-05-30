TOL <- 0.00001; 
nf <- 0; 
minimize <- TOL; #Tolerance limits
m <- 2000000; #Number of Monte Carlo Simulations
#Failure Terminated Data
t1 = 392; 
n=6; #Test Failure Count
result = vector("list")
time1 <- Sys.time()
for (i in 1:m){
      #Beta Random Generator
      beta <- runif(1, min = 0, max = 10)
      #Failure Times Random Generator
      t2 <- sample(327:624,1);
      t3 <- sample(327:624,1);
      t4 <- sample(327:624,1);
      t5 <- sample(327:624,1);
      t6 <- sample(327:624,1);
      lnt1 <- log(t1);
      lnt2 <- log(t2);
      lnt3 <- log(t3);
      lnt4 <- log(t4);
      lnt5 <- log(t5);
      lnt6 <- log(t6);
      sumlnti <- lnt1+lnt2+lnt3+lnt4+lnt5+lnt6;
      sumtilnti <- ((t1^beta)*lnt1)+((t2^beta)*lnt2)+((t3^beta)*lnt3)+((t4^beta)*lnt4)+((t5^beta)*lnt5)+((t6^beta)*lnt6);

      sumti = t1^beta+t2^beta+t3^beta+t4^beta+t5^beta+t6^beta;
      minsum = sumtilnti/sumti-1/beta-sumlnti/n
      theta = (sumti/n)^(1/beta);
      optimumtheta = 327/((log(1/(1-0.0085)))^(1/beta));
      if (theta>optimumtheta){
            if (minsum>0 & minsum<TOL){
                  if (minsum<minimize){
                  nf=nf+1;
                  minimize=minsum;
                  minbeta = beta;
                  simt1=t1;
                  simt2=t2;
                  simt3=t3;
                  simt4=t4;
                  simt5=t5;
                  simt6=t6;
                  mintheta=(sumti/n)^(1/minbeta);
                  times <- mapply(c, simt1, simt2, simt3, simt4, simt5, simt6)
                  result[[m]] <- times
                  }
            }
      }
}
time2 <- Sys.time()
Simulation.Time <- difftime(time2, time1)
Simulation.Time
final.result <- result[!sapply(result, is.null)]
x <- c(simt1, simt2, simt3, simt4, simt5, simt6)
y <- dweibull(x, scale = mintheta, shape = minbeta)
z <- t(mapply(c, x, y))
a <- sample(327:800, 1000, replace = TRUE)
f <- dweibull(a, scale = mintheta, shape = minbeta)
g <- t(mapply(c, a, f))
dist <- rweibull(1000, scale = mintheta, shape = minbeta)
plot(z, ylim=c(0, 0.008), xlim= c(327, 800), pch = 16, col = "blue")
      points(g, col = "red", pch = 4, cex = 0.3)