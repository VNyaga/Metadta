# Comparison and validation of metadta for meta-analysis of diagnostic test accuracy studies.

Install the package

```{r mada}
#install.packages("mada", dependencies = TRUE)
library(mada)
```

## Section 5.1

Read in data

```{r telomerase}
telo <- read.csv('https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/telomerase.csv', sep=',')
```

Fit the normal-normal model

```{r telofit}
fit = reitsma(telo)
```

Transform to 0-1 scale

```{r telotransform}
ss = summary(fit)

#Table 2, column 7
options(digits=2) #set 2 dp
c(ss$coefficients[3,1], ss$coefficients[3,5], ss$coefficients[3,6]) #Sens
c(1-ss$coefficients[4,1], 1-ss$coefficients[4,6], 1-ss$coefficients[4,5]) #Spec
```

Obtain the variances table 2 column 7

```{r telovar}
ss$Psi
```

Forest plot in figure 2

```{r figure2}
par(mfcol=c(1,2), cex=0.75)
forest(madad(telo), type = "sens")
forest(madad(telo), type = "spec")
```

SROC plot in figure 4 (bottom middle)

```{r figure4}
plot(fit, sroclwd = 2, asp=0.5)
points(fpr(telo), sens(telo), pch = 2)
legend("bottomright", c("data", "summary estimate"), pch = c(2,1))
legend("bottomleft", c("SROC", "conf. region"), lwd = c(2,1))

```

## Section 5.2

```{r dementia}
dem <- read.csv('https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/dementia.csv', sep=',')

#Change names to caps
library(stringr)
names(dem)[2:5] <- str_to_upper(names(dem)[2:5])
```

Fit the normal-normal model

```{r demfit}
fit = reitsma(dem)

#Table 3, column 7
summary(fit)
```

Transform to 0-1 scale

```{r demtransform}
ss <- (summary(fit))

#Table 3, column 7
options(digits=2) #set 2 dp
c(ss$coefficients[3,1], ss$coefficients[3,5], ss$coefficients[3,6]) #Sens
c(1-ss$coefficients[4,1], 1-ss$coefficients[4,6], 1-ss$coefficients[4,5]) #Spec
```

Obtain the variances in table 3, column 7

```{r demvar}
ss$Psi
```

SROC plot in figure 5 (top middle)

```{r figure5}
plot(fit, sroclwd = 2, asp=0.5)
points(fpr(telo), sens(telo), pch = 2)
legend("bottomright", c("data", "summary estimate"), pch = c(2,1))
legend("bottomleft", c("SROC", "conf. region"), lwd = c(2,1))

```

## Section 5.3

Read in data

```{r ascus}
ascus <- read.csv('https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/ascus.csv', sep=',')
```

Fit the normal-normal model

```{r ascusfit}
(fit <- reitsma(ascus, formula = cbind(tsens, tfpr) ~ Test))
#Table 4, column 7
summary(fit)
```

Transform to 0-1 scale

```{r ascsustransform}
ss <- (summary(fit))

#Inverse function of log
expit <- function(x){exp(x)/(1 + exp(x))}

#Table 4, column 7
options(digits=2) #set 2 dp
c(expit(ss$coefficients[1,1]), expit(ss$coefficients[1,5]), expit(ss$coefficients[1,6])) #Sens HC2
expit(ss$coefficients[2,1]+ ss$coefficients[1,1]) #Sens RepC
c(1-expit(ss$coefficients[3,1]), 1-expit(ss$coefficients[3,6]), 1-expit(ss$coefficients[3,5])) #Spec HC2
1-expit(ss$coefficients[4,1] + ss$coefficients[3,1]) #Spec RepC
```

Obtain the variances in table 4 column 7

```{r ascusvar}
ss$Psi
```
