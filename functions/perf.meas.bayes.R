perf.bayes<-function(estim,se.estim,ul,ll,true.delta,n.sims){
delta<-se.delta<-   ul.delta<-  ll.delta<-c()

for (i in 1:n.sims){
  delta[i]    <- r.mu.cace[i]
  se.delta[i] <- r.se.cace[i]
  ll.delta[i] <- r.lb.cace[i]
  ul.delta[i] <- r.ub.cace[i]
}

z.test         <- qnorm(0.975)

coverage.delta <- true.delta >= t(ll.delta) & true.delta <= t(ul.delta)

bias.delta <-mean(delta)-true.delta
se.bias.delta<-mean(sd(delta)/sqrt(n.sims))
rmse.delta <- sqrt(mean((delta-true.delta)^2))
coverage.delta <- mean(coverage.delta)
monte.carlo.error=sqrt(var(delta))/sqrt(n.sims)
mean.ci.delta<- mean(ul.delta-ll.delta)
median.ci.delta<- median(ul.delta-ll.delta)


perf<-c(coverage.delta,bias.delta,se.bias.delta,rmse.delta, median.ci.delta)
names(perf)<-c("coverage","mean_bias","MC error bias","rmse","median ci width")

return(perf)

}
