#传染病模型

library(deSolve)
library(magrittr)
library(ggplot2)
library(reshape2)
setwd("C:/work/R/covid19")

# 画图函数
draw<-function(sdf,title="SIR模型",ylab="感染人数",xlab="时间周期",ylog=FALSE){
  g<-ggplot(data = sdf, mapping = aes(x=time,y=value,colour=type))
  g<-g+geom_line() + geom_point()                                   #绘制线图和点图
  g<-g+scale_shape_manual(values = c(21,23))                        #自定义点形状
  if(ylog) g<-g+scale_y_log10()
  g<-g+ggtitle(title)+xlab(xlab)+ylab(ylab)
  g
}

# 按行合并
rbinddf<-function(ldf=list(),names=NULL){
  n <- length(ldf)
  res <- NULL
  if(n>0){
    for (i in seq(n)) {
      if(!is.null(names)){
        type <- names[i]
      }
      res <- rbind(res, data.frame(ldf[[i]],type=type))
    }
  }
  return(res)
}


##########################
# case1 自由增长模型
############################
calc1 <- function(rio = 0.3,         # 传染率系数,每个病人每天传染的人数
                  times = 0:20,      # 病情发展时间
                  I = 1) {           # 已感染者1人
  # 公式
  s_equations <- function(times, vars, params) {
    with(as.list(c(vars, params)), {
      dI <- rio * I
      return(list(dI))
    })
  }
  
  # 解微分方程
  ode(
    func = s_equations,
    y = c(I = I),
    times = times,
    parms = c(rio = rio)
  ) %>%  as.data.frame
}

s1 <- calc1(rio = 0.3)
s2 <- calc1(rio = 0.25)
s3 <- calc1(rio = 0.2)

sdf1<-rbinddf(list(s1,s2,s3),names=c("s1","s2","s3"))
head(sdf1)
names(sdf1)<-c("time","value","type")
draw(sdf1,title="自由增长模型",ylab="感染人数",xlab="时间周期")


##########################
# case2 SI模型
############################
calc2 <- function(rio = 0.3,                # 接触率,传染率系数
                 times = 0:100,             # 病情发展时间
                 i = 0.000001) {            # 已感染者,初始占比：百万分之一
  # 公式
  si_equations <- function(time, vars, params) {
    with(as.list(c(vars, params)), {
      di <- rio * i * (1 - i)
      return(list(di))
    })
  }
  
  # 解微分方程
  ode(
    func = si_equations,  
    y = c(i = i),         
    times = times,        
    parms = c(rio = rio)
  ) %>%  as.data.frame
}

si1 <- calc2(rio = 0.3)
si2 <- calc2(rio = 0.25)
si3 <- calc2(rio = 0.2)
si4 <- calc2(rio = 0.15)
si5 <- calc2(rio = 0.12)

sdf2<-rbinddf(list(si1,si2,si3,si4,si5),
             names=c("si1","si2","si3","si4","si5"))
names(sdf2)<-c("time","value","type")
draw(sdf2,title="SI模型",ylab="感染人数占比",xlab="时间周期")


##########################
# case3 SIS模型
############################
calc3 <- function(rio = 0.3,                 # 接触率,传染率系数
                  mu=0.1,                    # 治愈率
                  times = 0:150,             # 病情发展时间
                  i = 0.000001) {            # 已感染者初始占比：百万分之一
  sis_equations<-function(time,vars,params){
    with(as.list(c(vars,params)),{
      di<-rio*i*(1-i)-mu*i
      return(list(c(di)))
    })
  }
  
  ode(
    func=sis_equations,
    y=c(i=i),
    times=times,
    parms=c(rio=rio,mu=mu)
  ) %>% as.data.frame()
}


sis1 <- calc3(rio = 0.27,mu=0.1)
sis2 <- calc3(rio = 0.25,mu=0.1)
sis3 <- calc3(rio = 0.20,mu=0.1)
sis4 <- calc3(rio = 0.15,mu=0.1)
sis5 <- calc3(rio = 0.12,mu=0.1)

sdf<-rbinddf(list(sis1,sis2,sis3,sis4,sis5),
             names=c("sis1","sis2","sis3","sis4","sis5"))
names(sdf)<-c("time","value","type")
draw(sdf,title="SIS模型",ylab="感染人数占比",xlab="时间周期")



##########################
# case4 SIR模型
############################
calc4 <- function(rio = 0.3,                 # 接触率，传染率系数
                  mu=0.1,                    # 治愈率
                  times = 0:150,             # 病情发展时间
                  S=1000*1000*10,            # 易感染者1000万人
                  I=10,                      # 已感染者10人
                  R=5) {                     # 已移出者5人
  sir_equations<-function(time,vars,params){
    with(as.list(c(vars,params)),{
      N<-S+I+R
      dS<- -rio*I*S/N
      dI<-rio*I*S/N - mu*I
      dR<-mu*I
      return(list(c(dS,dI,dR)))
    })
  }
  ode(
    func=sir_equations,
    y=c(S=S,I=I,R=R),
    times=times,
    parms=c(rio=rio,mu=mu)
  ) %>% as.data.frame()
}

sir1<-calc4(rio=0.25,mu=0.1,times=0:150,S=1000*1000*10,I=10,R=100)
head(sir1)
sdf4<-melt(sir1,id.vars = c("time"))
head(sdf4)
names(sdf4)<-c("time","type","value")
draw(sdf4,title="SIR模型",ylab="感染人数",xlab="时间周期")























