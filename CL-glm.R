library(sp)
library(CASdatasets)
library(ChainLadder)
library(stats)
library(fitdistrplus)
library(MASS)

#Chain Ladder###################################################################
#transformasi data untuk chain ladder
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

mack <- MackChainLadder(data1, est.sigma="Mack")
mack$f #development factor
mack$FullTriangle
mack_smmry <- summary(mack)
mack_smmry$ByOrigin[4]
mack_smmry$ByOrigin[5]

mack2 <- MackChainLadder(data2, est.sigma="Mack")
mack2$f
mack2$FullTriangle
mack_smmry2 <- summary(mack2)
mack_smmry2$ByOrigin[4]
mack_smmry2$ByOrigin[5]

mack3 <- MackChainLadder(data3, est.sigma="Mack")
mack3$f
mack3$FullTriangle
mack_smmry3 <- summary(mack3)
mack_smmry3$ByOrigin[4]
mack_smmry3$ByOrigin[5]

#GLM Method#####################################################################
data1incr <- cum2incr(data1)
data2incr <- cum2incr(data2)
data3incr <- cum2incr(data3)

fit1incr <- glmReserve(data1incr, var.power = 2, cum = FALSE)
summary(fit1incr)
summary(fit1incr, type = "model")
fit1incr[["FullTriangle"]]

fit2incr <- glmReserve(data2incr, var.power = 2, cum = FALSE)
summary(fit2incr)
summary(fit2incr, type = "model")
fit2incr[["FullTriangle"]]

fit3incr <- glmReserve(data3incr, var.power = 2, cum = FALSE)
summary(fit3incr)
summary(fit3incr, type = "model") 
fit3incr[["FullTriangle"]]