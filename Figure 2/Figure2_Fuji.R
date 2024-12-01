rm(list=ls())

##Loading Downloaded data
data_raw=read.csv("Dome_Fuji_downloaded.csv")
## Age and T_source coulumns of the data corresponds to age and temperature 

##Loading data file (log(AAR), Temperature(x) and age(z)) derived from downloaded data  
data=read.csv("Dome_Fuji_data.csv")
#View(data)
data_raw$Age=sapply(data_raw$Age,as.numeric)
par(mar = c(6, 6, 4, 6))
plot(   data_raw$Age/1000,data_raw$T_source,type ="l",lwd=1.5, 
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
plot(data$z,exp(data$y), lwd=1.5,
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
mtext("AAR (m/Kyr)", side = 2, line = 3,
      cex=1.8)
legend("topright",
       c("Temperature","AAR "),
       col=c("#CC0033","#00BFFF"),
       lty = c(2,1),
       lwd=c(1.5,1.5),
       bty="n",
       cex=1.2
)

