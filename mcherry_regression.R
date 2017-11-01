#reading data
table <- read_excel("C:/Users/Ali/Desktop/igem/Wiki Files/table2.xlsx")

#Vectorizing Data

time <- table$`RFU/OD600`[c(3:15)]

time <- as.numeric(time)


x1 <-table$X__1[c(3:15)]

x1 <- as.numeric(x1)


x2 <- table$X__2[c(3:15)]
x2 <- as.numeric(x2)


x3 <- table$X__3[c(3:15)]
x3 <- as.numeric(x3)

x <- c(x1,x2,x3)

time_ <- c(time,time,time)
#plotting data vs time
#plot(c(time,time,time), c(x1,x2,x3), xlab = 'Time', ylab = 'RFU/OD600')

#Transforming variable

log_x = log(x)



plot(c(time,time,time), log_x, xlab = 'Time', ylab = 'log(RFu/OD600)')

#regression model

fit <- lm(log(x) ~ c(time,time,time))

#regression information
summary(fit)


#graphing best fit line
abline(fit, col='red')


#orginal data points

plot(c(time,time,time), x, xlab='Time', ylab='RFu/OD600')

#transformed prediction line

time_val <- seq(min(time),max(time), by = 13/38)
  
#prediction
lm2 <- exp(predict(fit,list(time=time_val)))


#plotting prediction

lines(time_val, lm2[c(1:39)], col="red")
  
  

