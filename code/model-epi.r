library(EpiModel)
setwd("C:/work/R/covid19")

# SI基础模型
param <- param.dcm(inf.prob = 0.5, act.rate = 0.5)
init <- init.dcm(s.num = 500, i.num = 1)
control <- control.dcm(type = "SI", nsteps = 100)
mod <- dcm(param, init, control)
mod
plot(mod)
head(as.data.frame(mod),10)
comp_plot(mod, at = 20, digits = 2)
summary(mod, at = 50)

# SIR 人口分析
param <- param.dcm(inf.prob = 0.2, act.rate = 1, rec.rate = 1/20,
                   a.rate = 1/95, ds.rate = 1/100, di.rate = 1/80, dr.rate = 1/100)
init <- init.dcm(s.num = 1000, i.num = 1, r.num = 0)
control <- control.dcm(type = "SIR", nsteps = 500, dt = 0.5)
mod <- dcm(param, init, control)
mod
plot(mod)

par(mar = c(3.2, 3, 2, 1), mgp = c(2, 1, 0), mfrow = c(1, 2))
plot(mod, popfrac = FALSE, alpha = 0.5,
     lwd = 4, main = "Compartment Sizes")
plot(mod, y = "si.flow", lwd = 4, col = "firebrick",
     main = "Disease Incidence", legend = "n")

par(mfrow = c(1, 1))
comp_plot(mod, at = 49, digits = 1)

head(as.data.frame(mod))
tail(as.data.frame(mod)) # 总人数增加了


# SIS模型敏感性分析
inf<-c(0.1, 0.1, 0.2, 0.1, 0.2,0.2)
act<-c(0.2, 0.4, 0.4, 0.6, 0.6,1)
param <- param.dcm(inf.prob = inf, act.rate = act, rec.rate = 0.02)
init <- init.dcm(s.num = 500, i.num = 1)
control <- control.dcm(type = "SIS", nsteps = 350)
mod <- dcm(param, init, control)
mod
tail(as.data.frame(mod, run = 6)) #提取输出
head(as.data.frame(mod, run = 1)) #提取输出

par(mfrow = c(1,2), mar = c(3.2,3,2.5,1))
plot(mod, alpha = 1, main = "Disease Prevalence",legend = "full") #已感染者数量
plot(mod, y = "si.flow", col = "Greens", alpha = 0.8, 
     main = "Disease Incidence",legend = "full")


# 2组交叉模型
param <- param.dcm(inf.prob = 0.4,  inf.prob.g2 = 0.1, act.rate = 0.25, balance = "g1",
                   a.rate = 1/100, a.rate.g2 = NA, ds.rate = 1/100, ds.rate.g2 = 1/100,
                   di.rate = 1/50, di.rate.g2 = 1/50)
init <- init.dcm(s.num = 500, i.num = 1, s.num.g2 = 500, i.num.g2 = 0)
control <- control.dcm(type = "SI", nsteps = 500)
mod <- dcm(param, init, control)

par(mfrow = c(1, 1))
plot(mod)


# 模拟流行病传播模型

nw <-network.initialize(n = 1000, directed = FALSE)
nw <-set.vertex.attribute(nw, "race", rep(0:1, each = 500))
formation <-~edges + nodefactor("race") + nodematch("race") +concurrent
target.stats <- c(250, 375, 225, 100)
coef.diss <-dissolution_coefs(dissolution = ~offset(edges), duration = 25)
est1 <- netest(nw, formation,target.stats, coef.diss, edapprox = TRUE)
#模型诊断
dx <- netdx(est1, nsims = 5, nsteps =500,
            nwstats.formula = ~edges + nodefactor("race", base = 0) + nodematch("race") + concurrent)


param <- param.net(inf.prob = 0.1,act.rate = 5, rec.rate = 0.02)
status.vector <- c(rbinom(500, 1, 0.1),rep(0, 500)) # 二项分布
status.vector <- ifelse(status.vector ==1, "i", "s")
init <- init.net(status.vector =status.vector)
control <- control.net(type ="SIS", nsteps = 500, nsims = 10, epi.by = "race")

#模型构建
sim1 <- netsim(est1, param, init,control)

#数据结果
head(as.data.frame(sim1), 10)

#获得网络状态
nw <- get_network(sim1, sim = 1)

#获得每个时间点的时间状态
get_transmat(sim1, sim = 1)

#绘制任意时间点的网络状态
par(mfrow = c(1,2), mar = c(0,0,1,0))
plot(sim1, type = "network", at =1, col.status = TRUE,
     main = "Prevalence at t1")
plot(sim1, type = "network", at =500, col.status = TRUE,
     main = "Prevalence at t500")



# SEIR


SEIR <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Population size
    num <- s.num + e.num + i.num + r.num
    
    # Effective contact rate and FOI from a rearrangement of Beta * c * D
    ce <- R0 / i.dur
    lambda <- ce * i.num/num
    
    dS <- -lambda*s.num
    dE <- lambda*s.num - (1/e.dur)*e.num
    dI <- (1/e.dur)*e.num - (1 - cfr)*(1/i.dur)*i.num - cfr*(1/i.dur)*i.num
    dR <- (1 - cfr)*(1/i.dur)*i.num
    
    # Compartments and flows are part of the derivative vector
    # Other calculations to be output are outside the vector, but within the containing list
    list(c(dS, dE, dI, dR, 
           se.flow = lambda * s.num,
           ei.flow = (1/e.dur) * e.num,
           ir.flow = (1 - cfr)*(1/i.dur) * i.num,
           d.flow = cfr*(1/i.dur)*i.num),
         num = num,
         i.prev = i.num / num,
         ei.prev = (e.num + i.num)/num)
  })
}

param <- param.dcm(R0 = 1.9, e.dur = 10, i.dur = 14, cfr = c(0.5, 0.7, 0.9))
init <- init.dcm(s.num = 1e6, e.num = 10, i.num = 0, r.num = 0,
                 se.flow = 0, ei.flow = 0, ir.flow = 0, d.flow = 0)
control <- control.dcm(nsteps = 500, dt = 1, new.mod = SEIR)
mod <- dcm(param, init, control)
mod

par(mfrow = c(2, 2))
plot(mod, y = "i.num", run = 2, main = "患病率")
plot(mod, y = "se.flow", run = 2, main = "发病率")
plot(mod, y = "i.num", main = "感染人数")
plot(mod, y = "i.prev", main = "感染百分比", ylim = c(0, 0.5), legend = "full")


