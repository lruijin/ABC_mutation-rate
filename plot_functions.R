library(R.matlab)
library(hetGP)
init_data <- readMat('init1.mat')
par(mfrow = c(1,2))
plot(init_data$theta.list[,1], apply(init_data$X,c(1,2),mean), pch = 20,
     xlab = "Delta", ylab = "Mean")
abline(v = log10(0.8),col = "red")
plot(init_data$theta.list[,2], apply(init_data$X,c(1,2),mean), pch = 20,
     xlab = "P", ylab = "Mean")
abline(v = -6, col = "red")

par(mfrow = c(1,2))
plot(init_data$theta.list[,1], apply(init_data$X,c(1,2),sd), pch = 20,
     xlab = "Delta", ylab = "Standard deviation")
abline(v = log10(0.8),col = "red")
plot(init_data$theta.list[,2], apply(init_data$X,c(1,2),sd), pch = 20,
     xlab = "P", ylab = "Standard deviation")
abline(v = -6, col = "red")
## ----------------- Mean ---------------------
y = apply(init_data$X,c(1,2),mean)
X = init_data$theta.list.s

gpfit <- mleHomGP(X,y,lower =c(0.1,0.1), upper = c(4,4), init = list(theta_init = c(1,1)))
par(mfrow = c(1,2))
X.fit <- as.matrix(expand.grid(seq(0,1,len=10),seq(0,1,len=10)), ncol =2)
y.fit <- predict(gpfit,X.fit)
image(seq(0,1,len=10),seq(0,1,len=10), matrix(y.fit$mean,ncol=10), xlab = "Delta",
      ylab = "P", main = "Mean")
points(1+log10(0.8), 0.5, col="green", pch = 4)
## ----------------- sd ---------------------
y = apply(init_data$X,c(1,2),sd)
X = init_data$theta.list.s

gpfit <- mleHomGP(X,y,lower =c(0.1,0.1), upper = c(4,4), init = list(theta_init = c(1,1)))
X.fit <- as.matrix(expand.grid(seq(0,1,len=10),seq(0,1,len=10)), ncol =2)
y.fit <- predict(gpfit,X.fit)
image(seq(0,1,len=10),seq(0,1,len=10), matrix(y.fit$mean,ncol=10), xlab = "Delta",
      ylab = "P", main = "Standard deviation")
points(1+log10(0.8), 0.5, col="green", pch = 4)

#------------------- Sample ---------------------
library(ggplot2)
library(ggpubr)
sample <- readMat("sample_981408.mat")
plot(sample$Sample[,1], pch = 20, type = "l", xlab = "Steps", ylab = "Delta")
plot(sample$Sample[,2], pch = 20, type = "l", xlab = "Steps", ylab = "P")

sum(diff(sample$Sample[5001:10000,1]) !=0)/5000
res <- sample$Sample[5001:1000,]
res <- data.frame(res)
colMeans(res)
apply(res,2,quantile, c(0.025, 0.5, 0.975))
names(res) <- c("x","y")
head(res)
apply(res,2,quantile, c(0.025,0.975))

p_2d = ggplot(res, aes(x=10^x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_point(aes(x= 10^x, y = y), col = "green", size = 0.8, shape = 1, alpha = 0.5)+
  geom_point(x = 0.8, y = -6, col = "red", shape = "cross")+
  theme(
    legend.position='none'
  )+ 
  xlab("Delta")+
  ylab("P")

p_delta <- ggplot(res, aes(x = 10^x)) +
  geom_histogram(aes(y = stat(density)), bins = 15)+
  geom_vline(xintercept = 0.8, col = "red") +
  theme_classic()+xlab(NULL)
p_p <- ggplot(res, aes(y = y)) +
  geom_histogram(aes(x = stat(density)), bins = 15)+
  geom_hline(yintercept = -6, col = "red") +
  theme_classic()+ylab(NULL)

ggarrange(p_delta,NULL,p_2d, p_p, ncol=2,nrow = 2,widths = c(3,1),heights= c(1,3),
          align = "hv")
#----------------------- Mean ----------------------------
library(R.matlab)
library(hetGP)
init_data <- readMat('init2.mat')
par(mfrow=c(2,2))
xNames = c("delta","p1","p2","tau")
for(i in 1:4){
  plot(init_data$theta.list.s[,i], apply(init_data$X,c(1,2),mean), pch = 20,
       xlab = xNames[i],ylab = "Mean")
}
for(i in 1:4){
  plot(init_data$theta.list.s[,i], apply(init_data$X,c(1,2),sd), pch = 20,
       xlab = xNames[i],ylab = "Sd")
}

plot(init_data$theta.list[,1], apply(init_data$X,c(1,2),mean), pch = 20)
## ----------------- Mean and sd ---------------------
y = apply(init_data$X,c(1,2),sd)
X = init_data$theta.list.s

gpfit <- mleHomGP(X,y,lower =c(0.1,0.1,0.1,0.1), upper = c(4,4,4,4), init = list(theta_init = c(1,1,1,1)))
X.fit <- as.matrix(cbind(log10(.8) + 1,expand.grid(seq(0,1,len=10),seq(0,1,len=10)),(4-0.1)/13.2), ncol =4)
y.fit <- predict(gpfit,X.fit)
image(seq(0,1,len=10),seq(0,1,len=10), matrix(y.fit$mean,ncol=10), xlab = "p1",
      ylab = "p2")
points((2+log10(5))/6,0.5,col = "green", pch = 4)

X.fit <- as.matrix(expand.grid(seq(0,1,len=10),seq(0,1,len=10)), ncol =2)
X.fit <- cbind(X.fit[,1], (2+log10(5))/6, X.fit[,2],0.9/13.2)
y.fit <- predict(gpfit,X.fit)
image(seq(0,1,len=10),seq(0,1,len=10), matrix(y.fit$mean,ncol=10), xlab = "delta",
      ylab = "p2")
points(log10(.8) + 1, 0.5,col = "green", pch = 4)

X.fit <- as.matrix(expand.grid(seq(0,1,len=10),seq(0,1,len=10)), ncol =2)
X.fit <- cbind(log10(.8) + 1, (2+log10(5))/6, X.fit[,2], X.fit[,1])
y.fit <- predict(gpfit,X.fit)
image(seq(0,1,len=10),seq(0,1,len=10), matrix(y.fit$mean,ncol=10), xlab = "tau",
      ylab = "p2")
points(3.9/13.2, 0.5,col = "green", pch = 4)

X.fit <- as.matrix(expand.grid(seq(0,1,len=10),seq(0,1,len=10)), ncol =2)
X.fit <- cbind(X.fit[,1], (2+log10(5))/6, 0.5, X.fit[,2])
y.fit <- predict(gpfit,X.fit)
image(seq(0,1,len=10),seq(0,1,len=10), matrix(y.fit$mean,ncol=10), xlab = "delta",
      ylab = "tau")
points(log10(.8) + 1,3.9/13.2,col = "green", pch = 4)

#------------------- Sample ---------------------


#### ---------------- For real data --------------------------
setwd("../Result/werngren_v1")
bacName = c("H37Rv", "E865", "E729", "E740", "E1221", "E1449", "Harlingen",
            "E26_95", "E80", "E55", "E26_94", "E3942", "E47")
burnin = 20000

for(i in bacName[12]){
  sample <- readMat(paste0("sample_a10_",i,"_1.mat"))
  sample1 <- readMat(paste0("sample_a10_",i,"_1_m1.mat"))
  sample1a <- readMat(paste0("sample_a10_",i,"_1_m1a.mat"))
  
  data <- readMat(paste0(i,".mat"))
  assign(paste0("MOM_",i),data$MOM)
  assign(paste0("p_init_",i),data$p.init)
  assign(paste0("MLE_",i),data$MLE)
  
  if(i == "H37Rv"){
    sample$sample = sample$Sample
  }
  assign(paste0("acc_",i),sum(diff(sample$sample[,1])!=0)/nrow(sample$sample))
  assign(paste0("acc_",i,"_m1"),sum(diff(sample1$sample)!=0)/length(sample1$sample))
  assign(paste0("acc_",i,"_m1a"),sum(diff(sample1a$sample[,1])!=0)/nrow(sample1a$sample))
  
  res <- sample$sample[-(1:burnin),]
  res[,1:3] <- 10^res[,1:3]
  res1 <- sample1$sample[-(1:burnin),]
  res1 <- 10^res1
  res1a <- sample1a$sample[-(1:burnin),]
  res1a <- 10^res1a
  
  assign(paste0("med_",i), apply(res,2, quantile, 0.5))
  assign(paste0("med_",i,"_m1"), quantile(res1,0.5))
  assign(paste0("med_",i,"_m1a"), apply(res1a,2, quantile, 0.5))
  
  assign(paste0("Q1_",i), apply(res,2, quantile, 0.025))
  assign(paste0("Q1_",i,"_m1"), quantile(res1,0.025))
  assign(paste0("Q1_",i,"_m1a"), apply(res1a,2, quantile, 0.025))
  
  assign(paste0("Q3_",i), apply(res,2, quantile, 0.975))
  assign(paste0("Q3_",i,"_m1"), quantile(res1,0.975))
  assign(paste0("Q3_",i,"_m1a"), apply(res1a,2, quantile, 0.975))
  
  assign(paste0("M_",i), apply(res,2, mean))
  assign(paste0("M_",i,"_m1"), mean(res1))
  assign(paste0("M_",i,"_m1a"), apply(res1a,2, mean))
  
  assign(paste0("Std_",i), apply(res,2, sd))
  assign(paste0("Std_",i,"_m1"), sd(res1))
  assign(paste0("Std_",i,"_m1a"), apply(res1a,2, sd))
  
  assign(paste0("mode_",i), apply(res,2, Modes))
  assign(paste0("mode_",i,"_m1"), Modes(res1))
  assign(paste0("mode_",i,"_m1a"), apply(res1a,2, Modes))
  
  assign(paste0("res_",i), res)
  assign(paste0("res_",i,"_m1"), res1)
  assign(paste0("res_",i,"_m1a"), res1a)
}

for(i in bacName[12]){
  cat(i,":\n")
  cat("Mutation rate:", get(paste0("p_init_",i)),"\n")
  cat("MOM:", signif(get(paste0("MOM_",i)),3),"\n")
  cat("Model 1: \n")
  cat("Mean: ",signif(get(paste0("M_",i,"_m1")),3), "\n")
  cat("Median: ",signif(get(paste0("med_",i,"_m1")),3), "\n")
  cat("Modes: ",signif(get(paste0("mode_",i,"_m1"))$modes[1],3), "\n")
  cat("Standard deviation: ", signif(get(paste0("Std_",i,"_m1")),3), "\n")
  cat("Confidence interval: (", signif(get(paste0("Q1_",i,"_m1")),3),",",
      signif(get(paste0("Q3_",i,"_m1")),3), ")", "\n")
}

for(i in bacName[12]){
  res <- get(paste0("res_",i,"_m1a"))
  res <- log10(res)
  kern.fit <- kde(res,xmin = c(log10(0.5),-11), 
                  xmax = c(log10(2),-7), 
                  gridsize = c(50,50))
  top = kern.fit$estimate[order(kern.fit$estimate,decreasing = T)[1]]
  idx = which(kern.fit$estimate == top,arr.ind = T)
  joint.mode_m1a <- c(10^kern.fit$eval.points[[1]][idx[1]],10^kern.fit$eval.points[[2]][idx[2]])
  
  assign(paste0("joint.mode_",i,"_m1a"), joint.mode_m1a)
}

for(i in bacName[12]){
  cat(i,":\n")
  
  cat("Model 1a: \n")
  cat("Mean: ",signif(get(paste0("M_",i,"_m1a")),3), "\n")
  cat("Median: ",signif(get(paste0("med_",i,"_m1a")),3), "\n")
  cat("Modes of delta: ",signif(get(paste0("mode_",i,"_m1a"))[[1]]$modes,3), "\n")
  cat("Modes of p: ",signif(get(paste0("mode_",i,"_m1a"))[[2]]$modes,3), "\n")
  cat("Joint mode:", signif(get(paste0("joint.mode_",i,"_m1a")),3), "\n")
  cat("Standard deviation: ", signif(get(paste0("Std_",i,"_m1a")),3), "\n")
  cat("Confidence interval of delta: (", signif(get(paste0("Q1_",i,"_m1a"))[1],3),",",
      signif(get(paste0("Q3_",i,"_m1a"))[1],3), ")", "\n")
  cat("Confidence interval of p: (", signif(get(paste0("Q1_",i,"_m1a"))[2],3),",",
      signif(get(paste0("Q3_",i,"_m1a"))[2],3), ")", "\n")
}
for(i in bacName[12]){
  res <- get(paste0("res_",i))
  res[,1:3] <- log10(res[,1:3])
  kern.fit <- kde(res,xmin = c(log10(0.5),-11,-11, 4), 
                  xmax = c(log10(2),-7,-7, 18), 
                  gridsize = c(4,15,15,7))
  top = kern.fit$estimate[order(kern.fit$estimate,decreasing = T)[1]]
  idx = which(kern.fit$estimate == top,arr.ind = T)
  joint.mode <- c(10^kern.fit$eval.points[[1]][idx[1]],10^kern.fit$eval.points[[2]][idx[2]],
                      10^kern.fit$eval.points[[3]][idx[3]],kern.fit$eval.points[[4]][idx[4]])
  
  assign(paste0("joint.mode_",i), joint.mode)
}

for(i in bacName[12]){
  cat(i,":\n")
  
  cat("Two-stage model: \n")
  cat("Mean: ",signif(get(paste0("M_",i)),3), "\n")
  cat("Median: ",signif(get(paste0("med_",i)),3), "\n")
  cat("Modes of delta: ",signif(get(paste0("mode_",i))[[1]]$modes,3), "\n")
  cat("Modes of p1: ",signif(get(paste0("mode_",i))[[2]]$modes,3), "\n")
  cat("Modes of p2: ",signif(get(paste0("mode_",i))[[3]]$modes,3), "\n")
  cat("Modes of tau: ",signif(get(paste0("mode_",i))[[4]]$modes,3), "\n")
  #cat("Joint mode:", signif(get(paste0("joint.mode_",i)),3), "\n")
  cat("Standard deviation: ", signif(get(paste0("Std_",i)),3), "\n")
  cat("Confidence interval of delta: (", signif(get(paste0("Q1_",i))[1],3),",",
      signif(get(paste0("Q3_",i))[1],3), ")", "\n")
  cat("Confidence interval of p1: (", signif(get(paste0("Q1_",i))[2],3),",",
      signif(get(paste0("Q3_",i))[2],3),")", "\n")
  cat("Confidence interval of p2: (", signif(get(paste0("Q1_",i))[3],3),",",
      signif(get(paste0("Q3_",i))[3],3),")", "\n")
  cat("Confidence interval of tau: (", signif(get(paste0("Q1_",i))[4],3),",",
      signif(get(paste0("Q3_",i))[4],3),")", "\n")
}

setwd("../Result/werngren")
setwd("../simulation/simu2")
bacName = c("H37Rv", "E865", "E729", "E740", "E1221", "E1449", "Harlingen",
            "E26_95", "E80", "E55", "E26_94", "E3942", "E47")
sample <- readMat("sample_post_H37Rv.mat")
sample <- readMat("sample1_1.mat")
varNames <- c("Delta","P1","P2","Tau")
par(mfrow = c(2,2))
for(i in 1:4){
  plot(sample$Sample[,i], pch = 20, type = "l", xlab = "Steps", ylab = varNames[i])
}

burnin = 5000

res <- sample$Sample[-(1:burnin),]
res <- data.frame(res)
res[,1] <- 10^res[,1]
#apply(res,2,quantile, c(0.025, 0.5, 0.975))
names(res) <- c("delta","p1","p2","tau")
#head(res)
#apply(res,2,quantile, c(0.025,0.975))
#apply(res,2,Modes)

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
p_2d_p2_p1 = ggplot(res, aes(x=p1, y=p2) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_gradientn(colours = jet.colors(7)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  #geom_point(aes(x= p1, y = p2), col = "green", size = 0.8, shape = 1, alpha = 0.5)+
  #geom_point(y = median(res$p2), x = median(res$p1), col = "red", shape = "cross")+
  theme(
    legend.position='right'
  )+ 
  ggtitle(expression(paste("(A) Joint posterior of ", log[10](p[1]), " and ", log[10](p[2]))))+
  xlab(expression(log[10](P[1])))+
  ylab(expression(log[10](P[2])))

p_2d_tau_p1 = ggplot(res, aes(x=p1, y=tau) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_gradientn(colours = jet.colors(7)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  #geom_point(aes(x= p1, y = p2), col = "green", size = 0.8, shape = 1, alpha = 0.5)+
  #geom_point(y = median(res$p2), x = median(res$p1), col = "red", shape = "cross")+
  theme(
    legend.position='right'
  )+ 
  ggtitle(expression(paste("(B) Joint posterior of ", log[10](p[1]), " and ", tau)))+
  xlab(expression(log[10](P[1])))+
  ylab(expression(tau))

p_2d_p2_tau = ggplot(res, aes(x=tau, y=p2) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_gradientn(colours = jet.colors(7)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  #geom_point(aes(x= p1, y = p2), col = "green", size = 0.8, shape = 1, alpha = 0.5)+
  #geom_point(y = median(res$p2), x = median(res$p1), col = "red", shape = "cross")+
  theme(
    legend.position='right'
  )+ 
  ggtitle(expression(paste("(C) Joint posterior of ", tau, " and ", log[10](p[2]))))+
  xlab(expression(tau))+
  ylab(expression(log[10](P[2])))

p_p1 <- ggplot(res, aes(x = p1)) +
  geom_histogram(aes(y = stat(density)), bins = 15,
                 color = "black", fill = "lightblue")+
  geom_density(lwd = 1.2, colour = 2)+
  #geom_vline(xintercept = -9, col = "red") +
  scale_y_continuous(expand = c(0,0))+
  theme(panel.grid = element_blank(),
        axis.line.x = element_line(size = 1),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank())+xlab(NULL)
p_p2_v <- ggplot(res, aes(y = p2)) +
  geom_histogram(aes(x = stat(density)), bins = 15,
                 color = "black", fill = "lightblue")+
  geom_density(lwd = 1.2, colour = 2)+
  #geom_hline(yintercept = -8, col = "red") +
  scale_x_continuous(expand = c(0,0))+
  theme(panel.grid = element_blank(),
        axis.text.y.right = element_line(size = 1),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank())+ylab(NULL)
p_tau <- ggplot(res, aes(x = tau)) +
  geom_histogram(aes(y = stat(density)), bins = 15,
                 color = "black", fill = "lightblue")+
  geom_density(lwd = 1.2, colour = 2)+
  #geom_vline(xintercept = 14, col = "red") +
  theme_classic()+ylab(NULL)
p_tau_v <- ggplot(res, aes(y = tau)) +
  geom_histogram(aes(x = stat(density)), bins = 15,
                 color = "black", fill = "lightblue")+
  geom_density(lwd = 1.2, colour = 2)+
  #geom_hline(yintercept = 14, col = "red") +
  theme_classic()+ylab(NULL)




p_2d_delta_p1 = ggplot(res, aes(x=delta, y=p1) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_point(aes(x= delta, y = p1), col = "green", size = 0.8, shape = 1, alpha = 0.5)+
  geom_point(x = 0.8, y = -9, col = "red", shape = "cross")+
  theme(
    legend.position='none'
  )+ 
  xlab(expression(delta))+
  ylab(expression(log[10](P[1])))


p_2d_delta_tau = ggplot(res, aes(x=delta, y=tau) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_point(aes(x= delta, y = tau), col = "green", size = 0.8, shape = 1, alpha = 0.5)+
  geom_point(x = 0.8, y = 14, col = "red", shape = "cross")+
  theme(
    legend.position='none'
  )+ 
  xlab(expression(delta))+
  ylab(expression(tau))



p_2d_delta_p2 = ggplot(res, aes(x=delta, y=p2) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_point(aes(x= delta, y = p2), col = "green", size = 0.8, shape = 1, alpha = 0.5)+
  geom_point(x = 0.8, y = -8, col = "red", shape = "cross")+
  theme(
    legend.position='none'
  )+ 
  xlab(expression(delta))+
  ylab(expression(log[10](P[2])))

p_delta <- ggplot(res, aes(x = delta)) +
  geom_histogram(aes(y = stat(density)), bins = 15)+
  geom_vline(xintercept = 0.8, col = "red") +
  theme_classic()+xlab(NULL)
p_delta_v <- ggplot(res, aes(y = delta)) +
  geom_histogram(aes(x = stat(density)), bins = 15)+
  geom_hline(yintercept = 0.8, col = "red") +
  theme_classic()+xlab(NULL)
p_p1 <- ggplot(res, aes(x = p1)) +
  geom_histogram(aes(y = stat(density)), bins = 15)+
  geom_vline(xintercept = -9, col = "red") +
  theme_classic()+ylab(NULL)
p_p1_v <- ggplot(res, aes(y = p1)) +
  geom_histogram(aes(x = stat(density)), bins = 15)+
  geom_hline(yintercept = -9, col = "red") +
  theme_classic()+ylab(NULL)
p_p2 <- ggplot(res, aes(x = p2)) +
  geom_histogram(aes(y = stat(density)), bins = 15)+
  geom_vline(xintercept = -8, col = "red") +
  theme_classic()+ylab(NULL)
p_p2_v <- ggplot(res, aes(y = p2)) +
  geom_histogram(aes(x = stat(density)), bins = 15)+
  geom_hline(yintercept = -8, col = "red") +
  theme_classic()+ylab(NULL)
p_tau <- ggplot(res, aes(x = tau)) +
  geom_histogram(aes(y = stat(density)), bins = 15)+
  geom_vline(xintercept = 14, col = "red") +
  theme_classic()+ylab(NULL)
p_tau_v <- ggplot(res, aes(y = tau)) +
  geom_histogram(aes(x = stat(density)), bins = 15)+
  geom_hline(yintercept = 14, col = "red") +
  theme_classic()+ylab(NULL)


ggarrange(p_delta,NULL,p_2d_delta_p2, p_p2_v, ncol=2,nrow = 2,widths = c(3,1),heights= c(1,3),
          align = "hv")

ggarrange(p_p1,NULL,p_2d_p1_p2, p_p2_v, ncol=2,nrow = 2,widths = c(3,1),heights= c(1,3),
          align = "hv")
ggarrange(p_p2,NULL,p_2d_p2_tau, p_tau_v, ncol=2,nrow = 2,widths = c(3,1),heights= c(1,3),
          align = "hv")
ggarrange(p_delta,NULL,p_2d_delta_p1, p_p1_v, ncol=2,nrow = 2,widths = c(3,1),heights= c(1,3),
          align = "hv")
ggarrange(p_p1,NULL,p_2d_p1_tau, p_tau_v, ncol=2,nrow = 2,widths = c(3,1),heights= c(1,3),
          align = "hv")
ggarrange(p_delta,NULL,p_2d_delta_tau, p_tau_v, ncol=2,nrow = 2,widths = c(3,1),heights= c(1,3),
          align = "hv")
## -------------- Find the MAPs for the joint dist.--------------------
library(MASS)
library(ks)
kern.fit <- kde(res,
                xmin = c(log10(0.5),-11,-11, 4), xmax = c(log10(2),-7,-7, 18), gridsize = c(50,50,50,50))
#plot(kern.fit, drawpoints = T)
top = kern.fit$estimate[order(kern.fit$estimate,decreasing = T)[1]]
top
idx = which(kern.fit$estimate == top,arr.ind = T)
kern.fit$eval.points[[1]][idx[1]]
kern.fit$eval.points[[2]][idx[2]]
kern.fit$eval.points[[3]][idx[3]]
top5 = kern.fit$estimate[order(kern.fit$estimate,decreasing = T)[1:100]]
idx = matrix(NA,100,3)
for(i in 1:100){
  idx[i,] <- which(kern.fit$estimate==top5[i], arr.ind = T)
}

kern.pred <- predict(kern.fit, x = c(p1_high[1],p2_high[1],tau_high[1]))

#kern.fit$estimate[idx[1],idx[2],idx[3]]
p1_high <- kern.fit$eval.points[[1]][idx[,1]]
p2_high <- kern.fit$eval.points[[2]][idx[,2]]
tau_high <- kern.fit$eval.points[[3]][idx[,3]]
rank = 100:1

##### ------------------ Real Data -----------------------------
p_2d_delta_p1 = ggplot(res, aes(x=delta, y=p1) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_point(aes(x= delta, y = p1), col = "green", size = 0.8, shape = 1, alpha = 0.5)+
  geom_point(x = post.mean[1], y = post.mean[2], col = "red", shape = "cross")+
  theme(
    legend.position='none'
  )+ 
  xlab(expression(delta))+
  ylab(expression(log[10](P[1])))

p_2d_p1_p2 = ggplot(res, aes(x=p1, y=p2) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_point(aes(x= p1, y = p2), col = "green", size = 0.8, shape = 1, alpha = 0.5)+
  geom_point(x = post.mean[2], y = post.mean[3], col = "red", shape = "cross")+
  theme(
    legend.position='none'
  )+ 
  xlab(expression(log[10](P[1])))+
  ylab(expression(log[10](P[2])))

p_2d_p2_tau = ggplot(res, aes(x=p2, y=tau) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_point(aes(x= p2, y = tau), col = "green", size = 0.8, shape = 1, alpha = 0.5)+
  geom_point(x = post.mean[3], y = post.mean[4], col = "red", shape = "cross")+
  theme(
    legend.position='none'
  )+ 
  xlab(expression(log[10](P[2])))+
  ylab(expression(tau))

p_2d_delta_tau = ggplot(res, aes(x=delta, y=tau) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_point(aes(x= delta, y = tau), col = "green", size = 0.8, shape = 1, alpha = 0.5)+
  geom_point(x = post.mean[1], y = post.mean[4], col = "red", shape = "cross")+
  theme(
    legend.position='none'
  )+ 
  xlab(expression(delta))+
  ylab(expression(tau))

p_2d_p1_tau = ggplot(res, aes(x=p1, y=tau) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_point(aes(x= p1, y = tau), col = "green", size = 0.8, shape = 1, alpha = 0.5)+
  geom_point(x = post.mean[2], y = post.mean[4], col = "red", shape = "cross")+
  theme(
    legend.position='none'
  )+ 
  xlab(expression(log[10](P[1])))+
  ylab(expression(tau))

p_2d_delta_p2 = ggplot(res, aes(x=delta, y=p2) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_point(aes(x= delta, y = p2), col = "green", size = 0.8, shape = 1, alpha = 0.5)+
  geom_point(x = post.mean[1], y = post.mean[3], col = "red", shape = "cross")+
  theme(
    legend.position='none'
  )+ 
  xlab(expression(delta))+
  ylab(expression(log[10](P[2])))

p_delta <- ggplot(res, aes(x = delta)) +
  geom_histogram(aes(y = stat(density)), bins = 15)+
  geom_vline(xintercept = post.mean[1], col = "red") +
  theme_classic()+xlab(NULL)
p_delta_v <- ggplot(res, aes(y = delta)) +
  geom_histogram(aes(x = stat(density)), bins = 15)+
  geom_hline(yintercept = post.mean[1], col = "red") +
  theme_classic()+xlab(NULL)
p_p1 <- ggplot(res, aes(x = p1)) +
  geom_histogram(aes(y = stat(density)), bins = 15)+
  geom_vline(xintercept = post.mean[2], col = "red") +
  theme_classic()+ylab(NULL)
p_p1_v <- ggplot(res, aes(y = p1)) +
  geom_histogram(aes(x = stat(density)), bins = 15)+
  geom_hline(yintercept = post.mean[2], col = "red") +
  theme_classic()+ylab(NULL)
p_p2 <- ggplot(res, aes(x = p2)) +
  geom_histogram(aes(y = stat(density)), bins = 15)+
  geom_vline(xintercept = post.mean[3], col = "red") +
  theme_classic()+ylab(NULL)
p_p2_v <- ggplot(res, aes(y = p2)) +
  geom_histogram(aes(x = stat(density)), bins = 15)+
  geom_hline(yintercept = post.mean[3], col = "red") +
  theme_classic()+ylab(NULL)
p_tau <- ggplot(res, aes(x = tau)) +
  geom_histogram(aes(y = stat(density)), bins = 15)+
  geom_vline(xintercept = post.mean[4], col = "red") +
  theme_classic()+ylab(NULL)
p_tau_v <- ggplot(res, aes(y = tau)) +
  geom_histogram(aes(x = stat(density)), bins = 15)+
  geom_hline(yintercept = post.mean[4], col = "red") +
  theme_classic()+ylab(NULL)


ggarrange(p_delta,NULL,p_2d_delta_p2, p_p2_v, ncol=2,nrow = 2,widths = c(3,1),heights= c(1,3),
          align = "hv")

ggarrange(p_p1,NULL,p_2d_p1_p2, p_p2_v, ncol=2,nrow = 2,widths = c(3,1),heights= c(1,3),
          align = "hv")
ggarrange(p_p2,NULL,p_2d_p2_tau, p_tau_v, ncol=2,nrow = 2,widths = c(3,1),heights= c(1,3),
          align = "hv")
ggarrange(p_delta,NULL,p_2d_delta_p1, p_p1_v, ncol=2,nrow = 2,widths = c(3,1),heights= c(1,3),
          align = "hv")
ggarrange(p_p1,NULL,p_2d_p1_tau, p_tau_v, ncol=2,nrow = 2,widths = c(3,1),heights= c(1,3),
          align = "hv")
ggarrange(p_delta,NULL,p_2d_delta_tau, p_tau_v, ncol=2,nrow = 2,widths = c(3,1),heights= c(1,3),
          align = "hv")
          