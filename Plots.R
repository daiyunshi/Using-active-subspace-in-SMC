par(mfcol=c(1,1))
plot(mean_value_M10[1:50000],type = 's',ylim = c(-1,1),col='yellow green',xlab = 'Runs',ylab = 'Mean Value')
par(new=T)
plot(mean_M10[1:50000],type = 's',ylim = c(-1,1),col='navy blue',xlab = 'Runs',ylab = 'Mean Value')


ASMH_eASMH <- data.frame(M=factor(rep(c("ASMH","eASMH"),each=10000)),Y=c(Y_ASMH_M10[90001:100000],RYM10[90001:100000]))
head(ASMH_eASMH)
ggplot(ASMH_eASMH, aes(x=Y, color=M)) + 
  geom_density(alpha=.5)

##²âÊÔMµÄÓ°Ïì
######
#ASMH#
######
par(mfcol=c(1,2))
plot(mean_value_M10[1:100000],type = 's',ylim = c(-1,1),col='seagreen',xlab = 'Runs',ylab = 'Mean Value')
plot(mean_value_M1[1:100000],type = 's',ylim = c(-1,1),col='skyblue4',xlab = 'Runs',ylab = 'Mean Value')

hist(Y_ASMH_M10[90000:100000],breaks = 100,col = 'seagreen' ,probability = T,xlab = 'y',ylim = c(0,0.8))
hist(Y_ASMH_M1[90000:100000],breaks = 100,col = 'skyblue4',probability = T ,xlab = 'y',ylim = c(0,0.8))
par(mfcol=c(1,2))
acf3 <- acf(Y_ASMH_M10_2.38[90000:100000],lag.max = 25)
acf4 <- acf(Y_ASMH_M1_2.38[90000:100000],lag.max = 25)

plot(Y_ASMH_M10[60000:70000],type = 's',col='seagreen')
plot(Y_ASMH_M1[60000:70000],type = 's',col='skyblue4')
#######
#eASMH#
#######

hist(RYM10[10000:100000],breaks = 100,col = 'seagreen' ,probability = T,xlab = 'y',ylim = c(0,0.8))
hist(RYM1[10000:100000],breaks = 100,col = 'skyblue4',probability = T ,xlab = 'y',ylim = c(0,0.8))

plot(RYM10[65000:70000],type = 's',col='seagreen')
plot(RYM1[65000:70000],type = 's',col='skyblue4')

plot(mean_M10[1:50000],type = 's',ylim = c(-1,1),col='seagreen',xlab = 'Runs',ylab = 'Mean Value')
plot(mean_M1[1:50000],type = 's',ylim = c(-1,1),col='skyblue4',xlab = 'Runs',ylab = 'Mean Value')

par(mfcol=c(2,1))
plot(RYM10_2.38[1,c(600000:610000)],type = 's',col="salmon",ylim = c(-4,4),xlab = 'Iteration',ylab = 'Markov chain states')
par(new=T)
plot(RYM10_2.38[2,c(600000:610000)],type = 's',col="skyblue",ylim = c(-4,4),xlab = 'Iteration',ylab = 'Markov chain states')

plot(RYM1_2.38[1,c(60000:70000)],type = 's',col="salmon",ylim = c(-4,4),xlab = 'Iteration',ylab = 'Markov chain states')
par(new=T)
plot(RYM1_2.38[2,c(60000:70000)],type = 's',col="skyblue",ylim = c(-4,4),xlab = 'Iteration',ylab = 'Markov chain states')

par(mfcol=c(1,3))
acf5 <- acf(RYM100_2.38[60000:100000],lag.max = 25)
acf7 <- acf(RYM10_2.38[60000:100000],lag.max = 25)
acf6 <- acf(RYM1_2.38[60000:100000],lag.max = 25)
par(mfcol=c(1,1))
plot(acf5,type='s')
par(new=T)
plot(acf6,type='s')
par(new=T)
plot(acf7,type='s')




REX2.38_mean1 <- c()
for(i in 1:100000){
  REX2.38_mean1[i] <- mean(REX2.38[1,1:i])
}
REX2.38_mean2 <- c()
for(i in 1:100000){
  REX2.38_mean2[i] <- mean(REX2.38[2,1:i])
}
plot(REX2.38_mean1,type='s',col='salmon',ylim = c(-1,1))
par(new=T)
plot(REX2.38_mean2,type='s',col='skyblue',ylim = c(-1,1))
REXEASMH1 <- c()
for(i in 1:100000){
  REXEASMH1[i] <- mean(REX2.38[1,1:i])
}
REXEASMH2 <- c()
for(i in 1:100000){
  REXEASMH2[i] <- mean(REX2.38[2,1:i])
}
plot(REXEASMH1,type='s',col='salmon',ylim = c(-1,1))
par(new=T)
plot(REXEASMH2,type='s',col='skyblue',ylim = c(-1,1))

######
#eq2##
######
plot(mean_eq2_M10,type='s',col='salmon',ylim = c(-1,1),ylab = 'Mean value',xlab = 'Runs')
par(new=T)
plot(mean_value_M10,type='s',col='skyblue',ylim = c(-1,1),ylab = 'Mean value',xlab = 'Runs')

par(mfcol=c(1,2))
plot(eq2_M10[20000:30000],type = 's')
plot(Y_ASMH_M10[20000:30000],type = 's')

acfeq2ASMH <- acf(eq2_M10,lag.max = 100)
acfeq1ASMH <- acf(Y_ASMH_M10,lag.max = 100)
