#######################
# 获取数据
######################

library(tidyverse)
library(lubridate)

#remotes::install_github("GuangchuangYu/nCov2019")
library(nCov2019)
x <- load_nCov2019()
plot(x)
#saveRDS(x,"x.rds")

str(x)
chinadf<-x$province[which(x$province$country=='China'),]
head(chinadf)
beijingdf<-chinadf[which(chinadf$province=='beijing'),]
head(beijingdf)

#install.packages(c("shiny","shinycssloaders","shinydashboard","plotly"))
open_dashboard()


##############################
# 北京疫情模拟
################################

# 每日新增
beijingdf$diff<-c(0,diff(beijingdf$cases))
plot(diff~date,data=beijingdf)

bj202204<-beijingdf[which(beijingdf$date>=as.Date("2022-04-01")),]
head(bj202204)
bj202204$cum<-cumsum(bj202204$diff) # 累计
plot(diff~date,data=bj202204,type='b')
plot(cum~date,data=bj202204,type='b')


##########
#基本公式模型
##########


#日接触率，治愈率=1/14
sir1<-calc4(rio=0.15,mu=1/14,times=0:150,S=30*1000*1000,I=8,R=0)
sir1$time<-as.Date("2022-04-01")+0:150
head(sir1)
sir2<-merge(sir1,bj202204[,c("date","cum")],by.x="time",by.y="date",all = TRUE)

sdf4<-melt(sir2,id.vars = c("time"))
names(sdf4)<-c("time","type","value")

g<-ggplot(data = sdf4, mapping = aes(x=time,y=value,colour=type))
g<-g+geom_line() + geom_point()
g<-g+scale_y_log10()
g<-g+ggtitle("SIR模型")+xlab("时间周期")+ylab("感染人数")
g


##########################
# beijing
#########################
Dat <- bj202204$cum
N <- 30*1000*1000
StartDay <- 6

# 有效数据
Infected <- Dat[StartDay:length(Dat)]
# 有效数据天数
Days <- 1:(length(Infected))
# 以多少天为窗口计算R0
WatchingWindow <- 12

# 设置SIR模型的默认值
init <- c(S = N - Infected[1], I = Infected[1], R = 0)

sir_equations<-function(time,vars,params){
  with(as.list(c(vars,params)),{
    N<-S+I+R
    dS<- -rio*I*S/N
    dI<-rio*I*S/N - mu*I
    dR<-mu*I
    return(list(c(dS,dI,dR)))
  })
}

#定义一些损失函数，由于后续采用的是最小二乘，所以都是适用与回归的损失函数
MSE <- function(Infected, fit) {
  sum((Infected - fit)^2) / length(Infected)
}
MAE <- function(Infected, fit) {
  sum(abs(Infected - fit)) / length(Infected)
}

R2 <- function(Infected, fit) {
  1 - MSE(Infected, fit)/var(Infected)
}

R0 <- function(days){
  
  # 取窗口期内的数据拟合
  WindowInfected <- Infected[(days - WatchingWindow):(days - 1)]
  Days <- 1:(length(WindowInfected))
  init <- c(S = N - WindowInfected[1], I = WindowInfected[1], R = 0)
  
  LOSS <- function(parameters) {
    names(parameters) <- c("rio", "mu")
    out <- ode(y = init, times = Days, func = sir_equations, parms = parameters)
    fit <- out[ , 3]
    -R2(WindowInfected,fit)
    # MSE(WindowInfected, fit)
  }
  
  Opt = optim(c(0.5, 0.5), LOSS, method = "L-BFGS-B", lower = c(0, 0), upper = c(Inf, 1)) 
  Opt_par = setNames(Opt$par, c("beta", "gamma"))
  as_tibble(list(days = days, 
                 R0 = Opt_par["beta"] / Opt_par["gamma"],
                 beta = Opt_par["beta"], 
                 gamma = Opt_par["gamma"], 
                 value = Opt$value,
                 date = ymd("2022-04-01") + days(StartDay + days - 1) - 1)) # StartDay - 1 | days - 1
}

RODF <- map_df((WatchingWindow + 1):length(Infected), R0)
RODF

plot(R0~date,data=RODF)

par(mfrow = c(1, 2))
plot(diff~date,data=bj202204,type='b')
plot(R0~date,data=RODF,type='b')
