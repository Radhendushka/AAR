rm(list=ls())

##Loading Downloaded data calibrated with AICC age
data_raw=read.csv("EPICA_downloaded.csv")

##Loading data file (log(AAR), Temperature(x) and age(z)) derived from downloaded data  
data=read.csv("EPICA_data.csv")

par(mar = c(6, 6, 4, 6))
plot(   data_raw$AgeAICC/1000,data_raw$Temperature,type ="l",lwd=1.5, 
        ylab = "",
        main = " ", 
        xlab = "",
        col = "#CC0033",
        yaxt="n",
        xaxt="n",
        ylim = c(-12,12),
        xlim = c(0,810),
        cex.axis=1.5
)
axis(side=1,at = seq(0,810,100),cex.axis=1.5)
mtext("Age (kyr)", side = 1, line = 3,cex=1.8)

axis(side=4,at = seq(-12,12,4),cex.axis=1.5)
mtext("Temperature (degree Celsius)", side = 4, line = 3,cex=1.8)

par(new = TRUE)
plot(data$z,exp(data$y), lwd=2,
     type = "l", 
     xaxt = "n", 
     #yaxt = "n",
     ylab = "", 
     xlab = "", 
     col = "#00BFFF", lty = 1,
     ylim = c(0.2,60),
     xlim = c(0,810),
     log="y",
     cex.axis=1.5
)
#axis(side=4,at = seq(-1,4.5,1))
mtext("AAR (m/kyr)", side = 2, line = 3,
      cex=1.8)
legend("topright",
       c("Temperature","AAR "),
       col=c("#CC0033","#00BFFF"),
       lty = c(1,1),
       lwd=c(1.5,2),
       bty="n",
       cex=1.2
)

