library(sp)
library(CASdatasets)
library(ChainLadder)
library(stats)
library(fitdistrplus)
library(MASS)
library(copula)

rm(list=ls())
data(ustri2GL)

data1<- matrix(unlist(ustri2GL$comauto), ncol = 10)
data1<- as.triangle(data1,
                    origin="origin",
                    dev="dev",
                    value="value")
data2<- matrix(unlist(ustri2GL$home), ncol = 10)
data2<- as.triangle(data2,
                    origin="origin",
                    dev="dev",
                    value="value")
data3<- matrix(unlist(ustri2GL$work), ncol = 10)
data3<- as.triangle(data3,
                    origin="origin",
                    dev="dev",
                    value="value")

data1incr <- cum2incr(data1)
data2incr <- cum2incr(data2)
data3incr <- cum2incr(data3)

fit1incr <- glmReserve(data1incr, var.power = 2, cum = FALSE)
summary(fit1incr)
summary(fit1incr, type = "model") 
tri1 <- as.data.frame(fit1incr[["FullTriangle"]])

fit2incr <- glmReserve(data2incr, var.power = 2, cum = FALSE)
summary(fit2incr)
summary(fit2incr, type = "model") 
tri2 <- as.data.frame(fit2incr[["FullTriangle"]])

fit3incr <- glmReserve(data3incr, var.power = 2, cum = FALSE)
summary(fit3incr)
summary(fit3incr, type = "model") 
tri3 <- as.data.frame(fit3incr[["FullTriangle"]])

data1incr <- as.data.frame(data1incr)
data2incr <- as.data.frame(data2incr)
data3incr <- as.data.frame(data3incr)

data1incr <- na.omit(data1incr)
data2incr <- na.omit(data2incr)
data3incr <- na.omit(data3incr)

data1incr$origin <- as.factor(data1incr$origin)
data2incr$origin <- as.factor(data2incr$origin)
data3incr$origin <- as.factor(data3incr$origin)

data1incr$dev <- as.factor(data1incr$dev)
data2incr$dev <- as.factor(data2incr$dev)
data3incr$dev <- as.factor(data3incr$dev)

gamma_glm1 <- glm(value ~ origin + dev, family = "Gamma"(link=log), data = data1incr)
alpha1 = gamma.shape(gamma_glm1)
tri1$origin <- as.numeric(tri1$origin)
tri1down <- tri1[tri1$origin + tri1$dev > 11, ]

gamma_glm2 <- glm(value ~ origin + dev, family = "Gamma"(link=log), data = data2incr)
alpha2 = gamma.shape(gamma_glm2)
tri2$origin <- as.numeric(tri2$origin)
tri2down <- tri2[tri2$origin + tri2$dev > 11, ]

gamma_glm3 <- glm(value ~ origin + dev, family = "Gamma"(link=log), data = data3incr)
alpha3 = gamma.shape(gamma_glm3)
tri3$origin <- as.numeric(tri3$origin)
tri3down <- tri3[tri3$origin + tri3$dev > 11, ]

#probability
dataa <- pobs(data1incr$value)
datab <- pobs(data2incr$value)
datac <- pobs(data3incr$value)
dataabc <- cbind(dataa, datab, datac)
fit1copula <- fitCopula(normalCopula(dim=3, dispstr = "un"), data = dataabc, method = "ml")
fit2copula <- fitCopula(claytonCopula(dim=3), data = dataabc, method = "ml")
fit3copula <- fitCopula(gumbelCopula(dim=3), data = dataabc, method = "ml")
fit4copula <- fitCopula(frankCopula(dim=3), data = dataabc, method = "ml")
fit5copula <- fitCopula(joeCopula(dim=3), data = dataabc, method = "ml")

#aic bic copula terbaik#########################################################
aic_fit1 <- -2 * fit1copula@loglik + 2 * length(fit1copula@estimate)
aic_fit2 <- -2 * fit2copula@loglik + 2 * length(fit2copula@estimate)
aic_fit3 <- -2 * fit3copula@loglik + 2 * length(fit3copula@estimate)
aic_fit4 <- -2 * fit4copula@loglik + 2 * length(fit4copula@estimate)
aic_fit5 <- -2 * fit5copula@loglik + 2 * length(fit5copula@estimate)
aic_fit1
aic_fit2
aic_fit3
aic_fit4
aic_fit5

bic_fit1 <- -2 * fit1copula@loglik + length(fit1copula@estimate) * log(length(dataabc))
bic_fit2 <- -2 * fit2copula@loglik + length(fit2copula@estimate) * log(length(dataabc))
bic_fit3 <- -2 * fit3copula@loglik + length(fit3copula@estimate) * log(length(dataabc))
bic_fit4 <- -2 * fit4copula@loglik + length(fit4copula@estimate) * log(length(dataabc))
bic_fit5 <- -2 * fit5copula@loglik + length(fit5copula@estimate) * log(length(dataabc))
bic_fit1
bic_fit2
bic_fit3
bic_fit4
bic_fit5

#simulasi 1000 kali#############################################################

random <- list()
prob_1 <- list()
prob_2 <- list()
prob_3 <- list()
beta_1 <- list()
beta_2 <- list()
beta_3 <- list()
tri1_fix <- list()
tri2_fix <- list()
tri3_fix <- list()
set.seed(100)
for (i in 1:1000) {
  random <- rCopula(45, fit2copula@copula)
  random <- as.data.frame(random)
  prob_a <- random$V1
  prob_b <- random$V2
  prob_c <- random$V3
  
  beta_a <- tri1down$value / alpha1$alpha
  tri1fix <- qgamma(prob_a, shape = alpha1$alpha, scale = beta_a)
  
  beta_b <- tri2down$value / alpha2$alpha
  tri2fix <- qgamma(prob_b, shape = alpha2$alpha, scale = beta_b)
  
  beta_c <- tri3down$value / alpha3$alpha
  tri3fix <- qgamma(prob_c, shape = alpha3$alpha, scale = beta_c)
  
  # Save the results in lists with different indices
  random[[i]] <- random
  prob_1[[i]] <- prob_a
  prob_2[[i]] <- prob_b
  prob_3[[i]] <- prob_c
  beta_1[[i]] <- beta_a
  beta_2[[i]] <- beta_b
  beta_3[[i]] <- beta_c
  tri1_fix[[i]] <- tri1fix
  tri2_fix[[i]] <- tri2fix
  tri3_fix[[i]] <- tri3fix
}

tri1_fix <- matrix(unlist(tri1_fix), ncol = 1000)
tri1_fix <- as.data.frame(tri1_fix)
tri1_down <- rowMeans(tri1_fix)
tri1down$value <- tri1_down
tri1_final <- tri1
tri1_final[tri1_final$origin+tri1_final$dev > 11,] <- tri1down
tri1_final <- as.triangle(tri1_final,
                          origin="origin",
                          dev="dev",
                          value="value")
tri1incr_final <- tri1_final
tri1_final <- incr2cum(tri1_final)
tri1_final <- as.data.frame(tri1_final)
tri1_final$origin <- as.numeric(tri1_final$origin)
paidtodate1 <- subset(tri1_final, origin + dev == 11, select = value)
paidtodate1$value <- rev(paidtodate1$value)
ultimateclaims1 <- subset(tri1_final, dev == 10, select = value)
IBNR1 <- ultimateclaims1-paidtodate1
hasil1 <- as.data.frame(c(ultimateclaims1, paidtodate1, IBNR1))
colnames(hasil1) <- c("ultimateclaims1", "paidtodate1", "IBNR1")
hasil1

tri2_fix <- matrix(unlist(tri2_fix), ncol = 1000)
tri2_fix <- as.data.frame(tri2_fix)
tri2_down <- rowMeans(tri2_fix)
tri2down$value <- tri2_down
tri2_final <- tri2
tri2_final[tri2_final$origin+tri2_final$dev > 11,] <- tri2down
tri2_final <- as.triangle(tri2_final,
                          origin="origin",
                          dev="dev",
                          value="value")
tri2incr_final <- tri2_final
tri2_final <- incr2cum(tri2_final)
tri2_final <- as.data.frame(tri2_final)
tri2_final$origin <- as.numeric(tri2_final$origin)
paidtodate2 <- subset(tri2_final, origin + dev == 11, select = value)
paidtodate2$value <- rev(paidtodate2$value)
ultimateclaims2 <- subset(tri2_final, dev == 10, select = value)
IBNR2 <- ultimateclaims2-paidtodate2
hasil2 <- as.data.frame(c(ultimateclaims2, paidtodate2, IBNR2))
colnames(hasil2) <- c("ultimateclaims2", "paidtodate2", "IBNR2")
hasil2

tri3_fix <- matrix(unlist(tri3_fix), ncol = 1000)
tri3_fix <- as.data.frame(tri3_fix)
tri3_down <- rowMeans(tri3_fix)
tri3down$value <- tri3_down
tri3_final <- tri3
tri3_final[tri3_final$origin+tri3_final$dev > 11,] <- tri3down
tri3_final <- as.triangle(tri3_final,
                          origin="origin",
                          dev="dev",
                          value="value")
tri3incr_final <- tri3_final
tri3_final <- incr2cum(tri3_final)
tri3_final <- as.data.frame(tri3_final)
tri3_final$origin <- as.numeric(tri3_final$origin)
paidtodate3 <- subset(tri3_final, origin + dev == 11, select = value)
paidtodate3$value <- rev(paidtodate3$value)
ultimateclaims3 <- subset(tri3_final, dev == 10, select = value)
IBNR3 <- ultimateclaims3-paidtodate3
hasil3 <- as.data.frame(c(ultimateclaims3, paidtodate3, IBNR3))
colnames(hasil3) <- c("ultimateclaims3", "paidtodate3", "IBNR3")
hasil3

#standar deviasi################################################################
tri1_final_1 <- list()
ultimateclaims_1 <- list()
IBNR_1 <- list()
tri1_fix <- matrix(unlist(tri1_fix), ncol = 1000)
tri1_fix <- as.data.frame(tri1_fix)
paidtodate1 <- subset(tri1_final, origin + dev == 11, select = value)
paidtodate1$value <- rev(paidtodate1$value)

for (i in 1:1000) {
  tri1_down_a <- tri1_fix[,i]
  tri1down$value <- tri1_down_a
  tri1_final_a <- tri1
  tri1_final_a[tri1_final_a$origin+tri1_final_a$dev > 11,] <- tri1down
  tri1_final_a <- as.triangle(tri1_final_a,
                              origin="origin",
                              dev="dev",
                              value="value")
  tri1_final_a <- incr2cum(tri1_final_a)
  tri1_final_a <- as.data.frame(tri1_final_a)
  tri1_final_a$origin <- as.numeric(tri1_final_a$origin)
  ultimateclaims1_a <- subset(tri1_final_a, dev == 10, select = value)
  IBNR1_a <- ultimateclaims1_a-paidtodate1
  
  tri1_final_1[[i]] <- tri1_final_a
  ultimateclaims_1[[i]] <- ultimateclaims1_a
  IBNR_1[[i]] <- IBNR1_a
}

IBNR_1
IBNR_1 <- as.data.frame(matrix(unlist(IBNR_1), ncol = 1000))
sd_1 <- apply(IBNR_1,1,sd)
mean_1 <- rowMeans(IBNR_1)

#standar deviasi tri2
tri2_final_2 <- list()
ultimateclaims_2 <- list()
IBNR_2 <- list()
tri2_fix <- matrix(unlist(tri2_fix), ncol = 1000)
tri2_fix <- as.data.frame(tri2_fix)
paidtodate2 <- subset(tri2_final, origin + dev == 11, select = value)
paidtodate2$value <- rev(paidtodate2$value)

for (i in 1:1000) {
  tri2_down_b <- tri2_fix[,i]
  tri2down$value <- tri2_down_b
  tri2_final_b <- tri2
  tri2_final_b[tri2_final_b$origin+tri2_final_b$dev > 11,] <- tri2down
  tri2_final_b <- as.triangle(tri2_final_b,
                              origin="origin",
                              dev="dev",
                              value="value")
  tri2_final_b <- incr2cum(tri2_final_b)
  tri2_final_b <- as.data.frame(tri2_final_b)
  tri2_final_b$origin <- as.numeric(tri2_final_b$origin)
  ultimateclaims2_b <- subset(tri2_final_b, dev == 10, select = value)
  IBNR2_b <- ultimateclaims2_b-paidtodate2
  
  tri2_final_2[[i]] <- tri2_final_b
  ultimateclaims_2[[i]] <- ultimateclaims2_b
  IBNR_2[[i]] <- IBNR2_b
}

IBNR_2 <- as.data.frame(matrix(unlist(IBNR_2), ncol = 1000))
sd_2 <- apply(IBNR_2,1,sd)
mean_2 <- rowMeans(IBNR_2)

#standar deviasi tri3
tri3_final_3 <- list()
ultimateclaims_3 <- list()
IBNR_3 <- list()
tri3_fix <- matrix(unlist(tri3_fix), ncol = 1000)
tri3_fix <- as.data.frame(tri3_fix)
paidtodate3 <- subset(tri3_final, origin + dev == 11, select = value)
paidtodate3$value <- rev(paidtodate3$value)

for (i in 1:1000) {
  tri3_down_c <- tri3_fix[,i]
  tri3down$value <- tri3_down_c
  tri3_final_c <- tri3
  tri3_final_c[tri3_final_c$origin+tri3_final_c$dev > 11,] <- tri3down
  tri3_final_c <- as.triangle(tri3_final_c,
                              origin="origin",
                              dev="dev",
                              value="value")
  tri3_final_c <- incr2cum(tri3_final_c)
  tri3_final_c <- as.data.frame(tri3_final_c)
  tri3_final_c$origin <- as.numeric(tri3_final_c$origin)
  ultimateclaims3_c <- subset(tri3_final_c, dev == 10, select = value)
  IBNR3_c <- ultimateclaims3_c-paidtodate3
  
  tri3_final_3[[i]] <- tri3_final_c
  ultimateclaims_3[[i]] <- ultimateclaims3_c
  IBNR_3[[i]] <- IBNR3_c
}

IBNR_3 <- as.data.frame(matrix(unlist(IBNR_3), ncol = 1000))
sd_3 <- apply(IBNR_3,1,sd)
mean_3 <- rowMeans(IBNR_3)
