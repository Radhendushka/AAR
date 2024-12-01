rm(list=ls())

library(fdrtool)
library(DescTools)
library(Iso)
library(Rfast)
library(gplm)
library(matrixStats)
## source of all estimation functions 
source("Estimation_functions.R")
## data consists of log AAR (y), temeprature (x) and age (z)
data=read.csv("Vostok_data.csv")

y=data$y
x=data$x
z=data$z


#Smoothing parameter
h=28

#Segmentwise Estimation

ind1=which(z<125)
y1=y[ind1]
x1=x[ind1]
z1=z[ind1]
S_Estimates1=Estimates(h,y1,x1,z1)
Bootstrap1=Estimates_Bootstrap(x1,z1,S_Estimates1$tilde_gamma,S_Estimates1$tilde_log_g,S_Estimates1$tilde_innov,S_Estimates1$tilde_phi)




ind2=which(z>125 & z<410)
y2=y[ind2]
x2=x[ind2]
z2=z[ind2]
S_Estimates2=Estimates(h,y2,x2,z2)
Bootstrap2=Estimates_Bootstrap(x2,z2,S_Estimates2$tilde_gamma,S_Estimates2$tilde_log_g,S_Estimates2$tilde_innov,S_Estimates2$tilde_phi)



##Combining g estimates of all three segments 

tilde_g1=exp(S_Estimates1$tilde_log_g)
tilde_g2=exp(S_Estimates2$tilde_log_g)
g=c(tilde_g1,tilde_g2)

##Adjusted and fitted AAR Scalled at g[1]
b1_adjusted=exp(y1)*g[1]/g[ind1]
b2_adjusted=exp(y2)*g[1]/g[ind2]


fit_b1_adjusted=(1+S_Estimates1$tilde_gamma*x1)*g[1]
fit_b2_adjusted=(1+S_Estimates2$tilde_gamma*x2)*g[1]


r1=cor(b1_adjusted,fit_b1_adjusted)
r2=cor(b2_adjusted,fit_b2_adjusted)

##Combined Adjusted AAR vector
b_adjusted=exp(y)*g[1]/g
fit_b_adjusted=c(fit_b1_adjusted,fit_b2_adjusted)

##predicted interval of fitted adjusted AAR scalled at g[1]
lower_limit_fit_b1_adjusted=Bootstrap1$predicted_Adjusted_AAR_quantile_025_star
lower_limit_fit_b2_adjusted=Bootstrap2$predicted_Adjusted_AAR_quantile_025_star
lower_limit_fit_b_adjusted=c(lower_limit_fit_b1_adjusted,lower_limit_fit_b2_adjusted)*g[1]

upper_limit_fit_b1_adjusted=Bootstrap1$predicted_Adjusted_AAR_quantile_975_star
upper_limit_fit_b2_adjusted=Bootstrap2$predicted_Adjusted_AAR_quantile_975_star
upper_limit_fit_b_adjusted=c(upper_limit_fit_b1_adjusted,upper_limit_fit_b2_adjusted)*g[1]


## Fitted y
fit_y1=log(1+S_Estimates1$tilde_gamma*x1)+S_Estimates1$tilde_log_g
fit_y2=log(1+S_Estimates2$tilde_gamma*x2)+S_Estimates2$tilde_log_g
fit_y=c(fit_y1,fit_y2)



##predicted interval of fitted log AAR
lower_limit_fit_y1=Bootstrap1$predicted_log_AAR_quantile_025_star
lower_limit_fit_y2=Bootstrap2$predicted_log_AAR_quantile_025_star
lower_limit_fit_y=c(lower_limit_fit_y1,lower_limit_fit_y2)

upper_limit_fit_y1=Bootstrap1$predicted_log_AAR_quantile_975_star
upper_limit_fit_y2=Bootstrap2$predicted_log_AAR_quantile_975_star
upper_limit_fit_y=c(upper_limit_fit_y1,upper_limit_fit_y2)



mean(b_adjusted>lower_limit_fit_b_adjusted & b_adjusted<upper_limit_fit_b_adjusted)
mean(b_adjusted[ind1]>lower_limit_fit_b_adjusted[ind1] & b_adjusted[ind1]<upper_limit_fit_b_adjusted[ind1])
mean(b_adjusted[ind2]>lower_limit_fit_b_adjusted[ind2] & b_adjusted[ind2]<upper_limit_fit_b_adjusted[ind2])


##Plots of figure 3
par(mar = c(6, 6, 4, 6))
plot(  z,b_adjusted, type ="l",lwd=2, 
       ylab = "",
       main = " ", 
       xlab = "",
       col = "#00BFFF",
       #yaxt="n",
       xaxt="n",
       xlim = c(-20,810),
       log="y",
       ylim = c(4,90),
       cex.axis=1.5
)


lines(z,lower_limit_fit_b_adjusted,
      ylab = "",
      main = " ", 
      xlab = "",
      col = "#ffc2ff",
      yaxt="n",
      xaxt="n",
      xlim = c(0,810),
      #xlim = c(10,150),
      lwd=0.1)
lines(z,upper_limit_fit_b_adjusted,
      ylab = "",
      main = " ", 
      xlab = "",
      col = "#ffc2ff",
      yaxt="n",
      xaxt="n",
      xlim = c(0,810),
      #xlim = c(10,150),
      lwd=0.1)

polygon(c(z, rev(z)), c(upper_limit_fit_b_adjusted, rev(lower_limit_fit_b_adjusted)),
        col = "#ffc2ff", lty = 0)


lines(  z,b_adjusted, type ="l",lwd=2, 
        ylab = "",
        main = " ", 
        xlab = "",
        col = "#00BFFF",
        #yaxt="n",
        xaxt="n",
        xlim = c(-20,810),
        log="y",
        ylim = c(4,100),
        cex.axis=1.5
)
lines(z,fit_b_adjusted,
      ylab = "",
      main = " ", 
      xlab = "",
      col = "#8B008B",
      yaxt="n",
      xaxt="n",
      xlim = c(0,810),
      #xlim = c(10,150),
      lwd=1.5)


abline(v=125)
abline(v=410)

axis(side=1,at = seq(0,810,100),cex.axis=1.5)
mtext("Age (kyr)", side = 1, line = 3,cex=1.8)
mtext(" Adjusted AAR (meter/kyr)", side = 2, line = 3,cex=1.8)
legend("bottomright",
       c("Observed"," Fitted"),
       col=c("#00BFFF","#8B008B"),
       lty = c(1,1),
       lwd=c(2,1.5),
       bty="n",
       cex=1.2
)
text(25,75,"ϱ=0.768",cex=1.2)
text(325,75,"ϱ=0.965",cex=1.2)


##Plots of Figure 4
par(mar = c(6, 6, 4, 6))
par(mfrow=c(1,3))
plot( fit_b_adjusted[ind1],b_adjusted[ind1],lwd=2, 
      ylab = "",
      main = " ", 
      xlab = "",
      col = "#00BFFF",
      #yaxt="n",
      #xaxt="n",
      xlim = c(5,50),
      #log="y",
      ylim = c(5,70),
      cex.axis=1.5,
      cex=0.2,pch=1
)
lines(10:50,10:50)

mtext("Adjusted AAR (meter/kyr)", side = 2, line = 3,cex=1.5)

text(10,70,"γ=0.059",cex=1.5)
text(10,65,"ϱ=0.768",cex=1.5)
plot( fit_b_adjusted[ind2],b_adjusted[ind2],lwd=2, 
      ylab = "",
      main = " ", 
      xlab = "",
      col = "#00BFFF",
      #yaxt="n",
      #xaxt="n",
      xlim = c(5,50),
      #log="y",
      ylim = c(5,70),
      cex.axis=1.5,
      cex=0.2,pch=1
      
)
lines(10:50,10:50)
text(10,70,"γ=0.056",cex=1.5)
text(10,65,"ϱ=0.965",cex=1.5)


