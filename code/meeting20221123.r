#################
# R语言大会：20221123
# 主题：用R语言解读传染病模型
####################


setwd("C:/workspace/doc/meeting/ms-build-2022/infect/code")

library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(magrittr)
library(lubridate)
library(ggrepel)

# 画图函数
draw<-function(df,ylog=TRUE){
  dat<-melt(df,id.vars = c("date"))
  g<-ggplot(data = dat, mapping = aes(x=date,y=value,colour=variable))
  g<-g+geom_line() + geom_point()   
  if(ylog) g<-g+scale_y_log10()
  g<-g+ggtitle("北京疫情统计")+xlab("日期")+ylab("人数")
  g
}
###################################

#remotes::install_github("YuLab-SMU/nCov2019")
library(nCov2019)
x <- query()
saveRDS(x,"x20221123.rds")
#x<-readRDS("x20221123.rds")
names(x)

# shiny 控制台
dashboard()

# 全球总体汇总统计
global<-x$global
summary(global)
global

#致死率
global$deaths/global$cases

# 全球所有国家的最新一天数据
last<-x$latest
head(last["Global"],10)
plot(x$latest)

last[c("USA","India","China")]
last$detail[which(last$detail$country=="China"),]

# 全球所有国家的最新一天数据
hist<-x$historical
head(hist["China"],10)
tail(hist["China"],10)
head(hist['China','beijing'])

# #目前疫苗研发进展
# vac <- x$vaccine
# summary(vac)
# head(vac["all"])
# 
# # 科兴疫苗
# vac["all"][which(vac["all"]$candidate=='CoronaVac'),]
# vac[ID="id5"]
# 
# #目前治疗进展
# thera<- x$therapeutics
# summary(thera)
# head(thera["All"])
# thera[ID="id30"] 
# 
# # 不同国家的数据
# tmp <- hist["global"] %>%
#   group_by(country) %>%
#   arrange(country,date) %>%
#   mutate(diff = cases - lag(cases, default =  first(cases))) %>%
#   filter(country %in% c("Australia", "Japan", "Italy", "Germany",  "China", "USA")) 
# 
# ggplot(tmp,aes(date, diff, color=country)) + geom_line() +
#   labs(y="daily increase cases") + 
#   theme(axis.text = element_text(angle = 15, hjust = 1)) +
#   scale_x_date(date_labels = "%Y-%m-%d") + 
#   theme_minimal()
# 
# # 扩展维度可视化
# # 中国累计确诊
# china <- hist['China']
# china <- china[order(china$cases), ]
# 
# ggplot(china, 
#        aes(date, cases)) +
#   geom_col(fill = 'firebrick') + 
#   theme_minimal(base_size = 14) +
#   xlab(NULL) + ylab(NULL) + 
#   scale_x_date(date_labels = "%Y/%m/%d") +
#   labs(caption = paste("accessed date:", max(china$date)))

##################################

hist<-x$historical
beijing<-hist["China","beijing"]
head(beijing)
tail(beijing)

# 计算每日新增
beijing$daily<-c(0,diff(beijing$cases))
draw(beijing[,c("date", "cases" ,"deaths", "recovered","daily")])

# 本次疫情从4月1日到7月31日数据
bj202204<-beijing[which(beijing$date>=as.Date("2022-04-01") & beijing$date<=as.Date("2022-07-31")),]
bj202204$cum<-cumsum(bj202204$daily) # 累计
head(bj202204)
draw(bj202204[,c("date", "cum","daily")],FALSE)

bj202204now<-beijing[which(beijing$date>=as.Date("2022-04-01")),]
bj202204now$cum<-cumsum(bj202204now$daily) # 累计

# 本次疫情从10月1日到11月20日数据
bj202210<-beijing[which(beijing$date>=as.Date("2022-09-01")),]
bj202210$cum<-cumsum(bj202210$daily) # 累计
head(bj202210)
draw(bj202210[,c("date", "cum","daily")],FALSE)


##########
# 用SIR基本公式模型
##########

library(EpiModel)
param <- param.dcm(inf.prob = 0.5, act.rate = 0.5, rec.rate = 1/7)
init <- init.dcm(s.num = 30*1000*1000, i.num = 8, r.num = 0)
control <- control.dcm(type = "SIR", nsteps = 200, dt = 0.5)
mod <- dcm(param, init, control)
mod
par(mfrow = c(1, 1))
plot(mod)

moddf<-as.data.frame(mod)
moddf$date<-as.Date("2022-04-01")+0:(nrow(moddf)-1)

d1<-data.frame(date=moddf$date,i=moddf$i.num,variable="sir")
d2<-data.frame(date=bj202204now$date,i=bj202204now$cum,variable="real")
mdf<-rbind(d1,d2)

g<-ggplot(data = mdf, mapping = aes(x=date,y=i,colour=variable))
g<-g+geom_line() + geom_point()
g<-g+scale_y_log10()
g<-g+ggtitle("北京疫情统计")+xlab("日期")+ylab("Log(人数)")
g

#########################
# 模型拟合：202210
#########################
N <- 30*1000*1000         # 总人数
window<-7                 # 传染期

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


R0<-function(df){
  cases<-df$cum
  init <- c(s.num =N-cases[1],i.num=cases[1],r.num=0)
  
  # 定义损失函数
  LOSS <- function(parameters) {
    names(parameters) <- c("beta", "gamma")
    out <- ode(y = init, times = 1:window, func = SIR, parms = parameters)
    fit <- out[ , 3]
    -r2(cases,fit)
  }
  
  # 拟合计算
  Opt = optim(c(0.5, 0.5), LOSS, method = "L-BFGS-B", lower = c(0, 0), upper = c(Inf, 1))
  Opt_par = setNames(Opt$par, c("beta", "gamma"))
  data.frame(
    date = tail(df$date,1)+days(1),
    R0 = Opt_par["beta"] / Opt_par["gamma"],
    beta = Opt_par["beta"],
    gamma = Opt_par["gamma"],
    value = Opt$value)
}


##########################
# 202204

n<-nrow(bj202204)
rodf<-lapply(14:(n-window),function(s){
  df<-bj202204[s:(s+window-1),]
  print(df$date[1])
  R0(df)
})%>% ldply


d1<-merge(bj202204,rodf,by="date",all.x=TRUE)
d1$control<-0.01
d1$control[which(d1$date>=as.Date('2022-04-25') & d1$date<=as.Date('2022-05-25'))]<-0.5
d1$Re<-d1$R0*(1-d1$control) * (1-d1$cum/N)

par(mfrow = c(4, 1))
plot(cum~date,data=d1,type='b',col="blue")
plot(daily~date,data=d1,type='b',col="red")
plot(R0~date,data=d1,ylim=c(1,1.5),type='b',col="orange")  # R0
plot(Re~date,data=d1,ylim=c(0,1.5),type='b',col="green")  # Re


##########################
# 202210
n<-nrow(bj202210)
rodf<-lapply(14:(n-window),function(s){
  df<-bj202210[s:(s+window-1),]
  print(df$date[1])
  R0(df)
})%>% ldply

d1<-merge(bj202210,rodf,by="date",all.x=TRUE)
d1$control<-c(rep(0.01,n-4),rep(0.5,4))
d1$Re<-d1$R0*(1-d1$control) * (1-d1$cum/N)

par(mfrow = c(4, 1))
plot(cum~date,data=d1,type='b',col="blue")
plot(daily~date,data=d1,type='b',col="red")
plot(R0~date,data=d1,ylim=c(1,1.5),type='b',col="orange")  # R0
plot(Re~date,data=d1,ylim=c(0,1.5),type='b',col="green")  # Re

###

# par(mfrow = c(3, 1))
# plot(cum~date,data=bj202204,type='b',col="blue")
# plot(daily~date,data=bj202204,type='b',col="red")
# 
# d2<-data.frame(date=rodf$date,R0=rodf$R0)
# d1<-data.frame(date=as.Date('2022-04-01')+0:16,R0=NA)
# d2<-rbind(d1,d2)
# plot(R0~date,data=d2,ylim=c(1,1.5),type='b',col="orange")



