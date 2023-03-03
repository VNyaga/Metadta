
#MADA
#install.packages("mada", dependencies = TRUE)
library(mada)

#============================================================================================
# 5.1
telo <- read.csv('https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/telomerase.csv', sep=',')

fit = reitsma(telo)

ss = summary(fit)
options(digits=2) #set 2 dp
c(ss$coefficients[3,1], ss$coefficients[3,5], ss$coefficients[3,6]) #Sens
c(1-ss$coefficients[4,1], 1-ss$coefficients[4,6], 1-ss$coefficients[4,5]) #Spec

ss$Psi

par(mfcol=c(1,2), cex=0.75)
forest(madad(telo), type = "sens")
forest(madad(telo), type = "spec")

plot(fit, sroclwd = 2, asp=0.5)
points(fpr(telo), sens(telo), pch = 2)
legend("bottomright", c("data", "summary estimate"), pch = c(2,1))
legend("bottomleft", c("SROC", "conf. region"), lwd = c(2,1))
#============================================================================================
#5.2
dem <- read.csv('https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/dementia.csv', sep=',')
library(stringr)
names(dem)[2:5] <- str_to_upper(names(dem)[2:5])
fit = reitsma(dem)
summary(fit)

ss <- (summary(fit))

#Table 3, column 7
options(digits=2) #set 2 dp
c(ss$coefficients[3,1], ss$coefficients[3,5], ss$coefficients[3,6]) #Sens
c(1-ss$coefficients[4,1], 1-ss$coefficients[4,6], 1-ss$coefficients[4,5]) #Spec
ss$Psi

plot(fit, sroclwd = 2, asp=0.5)
points(fpr(telo), sens(telo), pch = 2)
legend("bottomright", c("data", "summary estimate"), pch = c(2,1))
legend("bottomleft", c("SROC", "conf. region"), lwd = c(2,1))
#============================================================================================
#5.3
ascus <- read.csv('https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/ascus.csv', sep=',')
(fit <- reitsma(ascus, formula = cbind(tsens, tfpr) ~ Test))
ss <- (summary(fit))

#Inverse function of log
expit <- function(x){exp(x)/(1 + exp(x))}

#Table 4, column 7
options(digits=2) #set 2 dp
c(expit(ss$coefficients[1,1]), expit(ss$coefficients[1,5]), expit(ss$coefficients[1,6])) #Sens HC2
expit(ss$coefficients[2,1]+ ss$coefficients[1,1]) #Sens RepC
c(1-expit(ss$coefficients[3,1]), 1-expit(ss$coefficients[3,6]), 1-expit(ss$coefficients[3,5])) #Spec HC2
1-expit(ss$coefficients[4,1] + ss$coefficients[3,1]) #Spec RepC

ss$Psi