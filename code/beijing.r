#######################
# 获取数据
######################

setwd("C:/work/R/covid19")
#remotes::install_github("YuLab-SMU/nCov2019")

library(nCov2019)
x <- query()
#x<-readRDS("x.rds")

library(ggplot2)
library(magrittr)
library(reshape2)

# 画图函数
draw<-function(df){
  dat<-melt(df,id.vars = c("date"))
  g<-ggplot(data = dat, mapping = aes(x=date,y=value,colour=variable))
  g<-g+geom_line() + geom_point()   
  g<-g+scale_y_log10()
  g<-g+ggtitle("北京疫情统计")+xlab("日期")+ylab("Log(人数)")
  g
}

##############################
# 北京疫情模拟
################################

# 提取数据
hist<-x$historical
beijing<-hist["China","beijing"]
head(beijing)

# 计算每日新增
beijing$daily<-c(0,diff(beijing$cases))
draw(beijing[,c("date", "cases" ,"deaths", "recovered","daily")])

# 本次疫情从4月1日到7月4日数据
bj202204<-beijing[which(beijing$date>=as.Date("2022-04-01")),]
head(bj202204)
bj202204$cum<-cumsum(bj202204$daily) # 累计
draw(bj202204[,c("date", "cum","daily")])


##########
#基本公式模型
##########

# SIR
library(EpiModel)
param <- param.dcm(inf.prob = 0.5, act.rate = 0.5, rec.rate = 1/7)
init <- init.dcm(s.num = 30*1000*1000, i.num = 8, r.num = 0)
control <- control.dcm(type = "SIR", nsteps = 200, dt = 0.5)
mod <- dcm(param, init, control)
mod
plot(mod)

moddf<-as.data.frame(mod)
moddf$date<-as.Date("2022-04-01")+0:(nrow(moddf)-1)

d1<-data.frame(date=moddf$date,i=moddf$i.num,variable="sir")
d2<-data.frame(date=bj202204$date,i=bj202204$cum,variable="real")
mdf<-rbind(d1,d2)

g<-ggplot(data = mdf, mapping = aes(x=date,y=i,colour=variable))
g<-g+geom_line() + geom_point()   
g<-g+scale_y_log10()
g<-g+ggtitle("北京疫情统计")+xlab("日期")+ylab("Log(人数)")
g


##########
# 模型拟合
##########
library(lubridate)
library(plyr)

cases <- bj202204$cum     # 每日累计
N <- 30*1000*1000         # 总人数
startDay <- 6             # 开始累计时间

# 每日累计确诊
infected <- cases[startDay:length(cases)]
infected

# 以多少天为窗口计算R0
window <- 12


SIR<-function(time,vars,params){
  with(as.list(c(vars,params)),{ #beta=rio,gamma=mu
    num <- s.num + i.num + r.num
    dS<- -beta * i.num * s.num/num
    dI<- beta * i.num * s.num/num - gamma * i.num
    dR<- gamma * i.num
    return(list(c(dS,dI,dR)))
  })
}

#定义一些损失函数，由于后续采用的是最小二乘，所以都是适用与回归的损失函数
mse <- function(infected, fit) {
  sum((infected - fit)^2) / length(infected)
}

r2 <- function(infected, fit) {
  1 - mse(infected, fit)/var(infected)
}

R0 <- function(days){

  # 取窗口期内的数据拟合
  WindowInfected <- infected[(days - window):(days - 1)]
  # 持续天数
  Days <- 1:(length(WindowInfected))
  
  # 设置SIR模型的默认值，s.num易感人群，i.num已感人群，r.num已治愈人群
  init <- c(s.num = N - WindowInfected[1], i.num = WindowInfected[1], r.num = 0)

  # 定义损失函数
  LOSS <- function(parameters) {
    names(parameters) <- c("beta", "gamma")
    out <- ode(y = init, times = Days, func = SIR, parms = parameters)
    fit <- out[ , 3]
    -r2(WindowInfected,fit)
  }

  # 拟合计算
  Opt = optim(c(0.5, 0.5), LOSS, method = "L-BFGS-B", lower = c(0, 0), upper = c(Inf, 1))
  Opt_par = setNames(Opt$par, c("beta", "gamma"))
  data.frame(days = days,
             R0 = Opt_par["beta"] / Opt_par["gamma"],
             beta = Opt_par["beta"],
             gamma = Opt_par["gamma"],
             value = Opt$value,
             date = ymd("2022-04-01") + days(startDay + days - 1) - 1) # StartDay - 1 | days - 1
}

# 找到需要计算的日期
days<-(window+1):(length(infected));days

# 计算R0的计算
RODF<-lapply(days,R0) %>% ldply
head(RODF)

plot(R0~date,data=RODF)

par(mfrow = c(3, 1))
plot(cum~date,data=bj202204,type='b',col="blue")
plot(daily~date,data=bj202204,type='b',col="red")

d2<-data.frame(date=RODF$date,R0=RODF$R0)
d1<-data.frame(date=as.Date('2022-04-01')+0:16,R0=NA)
d2<-rbind(d1,d2)
plot(R0~date,data=d2,type='b',col="orange")

