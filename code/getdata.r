#######################
# 获取数据
######################

setwd("C:/work/R/covid19")
#remotes::install_github("YuLab-SMU/nCov2019")
library(nCov2019)
x <- query()
names(x)

# saveRDS(x,"x.rds")
# x<-readRDS("x.rds")

#######################
# 全球总体汇总统计
#######################
global<-x$global
summary(global)
global

#致死率
global$deaths/global$cases

#######################
# 全球所有国家的最新一天数据
#######################
last<-x$latest
head(last["Global"],10)

last[c("USA","India","China")]

last$detail[which(last$detail$country=="China"),]

#######################
# 全球所有国家的最新一天数据
#######################

hist<-x$historical
head(hist["China"],10)
tail(hist["China"],10)
head(hist['China','beijing'])

library(ggplot2)
library(reshape2)
beijing_df<-hist['China','beijing']
beijing<-melt(beijing_df[,c("date", "cases" ,"deaths", "recovered")],id.vars = c("date"))

g<-ggplot(data = beijing, mapping = aes(x=date,y=value,colour=variable))
g<-g+geom_line() + geom_point()                                   #绘制线图和点图
g<-g+scale_shape_manual(values = c(21,23))                        #自定义点形状
g<-g+scale_y_log10()
g<-g+ggtitle("北京疫情统计")+xlab("日期")+ylab("Log(人数)")
g

#########################
#目前疫苗研发进展
##########################

vac <- x$vaccine
summary(vac)
head(vac["all"])

# 科兴疫苗
vac["all"][which(vac["all"]$candidate=='CoronaVac'),]

vac[ID="id5"]

#########################
#目前治疗进展
##########################

thera<- x$therapeutics
summary(thera)
head(thera["All"])
thera[ID="id30"] 


#########################
#目前治疗进展
##########################

plot(last)
plot(last, type="tests",palette="Green")

library(ggplot2)
library(dplyr)
tmp <- hist["global"] %>%
  group_by(country) %>%
  arrange(country,date) %>%
  mutate(diff = cases - lag(cases, default =  first(cases))) %>%
  filter(country %in% c("Australia", "Japan", "Italy", "Germany",  "China")) 

ggplot(tmp,aes(date, log(diff+1), color=country)) + geom_line() +
  labs(y="Log2(daily increase cases)") + 
  theme(axis.text = element_text(angle = 15, hjust = 1)) +
  scale_x_date(date_labels = "%Y-%m-%d") + 
  theme_minimal()

# 历史数据画图
plot(hist, region="Global" ,date = "2020-08-01", type="cases")

# 动画效果
from = "2020-03-01"
to = "2020-04-01"
plot(hist, from = from, to=to)

#########################
# 扩展维度可视化
##########################

# 中国累计确诊
china <- hist['China']
china <- china[order(china$cases), ]

ggplot(china, 
       aes(date, cases)) +
  geom_col(fill = 'firebrick') + 
  theme_minimal(base_size = 14) +
  xlab(NULL) + ylab(NULL) + 
  scale_x_date(date_labels = "%Y/%m/%d") +
  labs(caption = paste("accessed date:", max(china$date)))


# 计算不同国家或地区新冠确珍人数每日增速

library(dplyr)
library(magrittr)
library(ggrepel)

country_list =  c("Italy","Brazil","Japan","China","USA","Mexico","Russia","India","Thailand")

hist[country_list]  %>%
  subset( date > as.Date("2020-10-01") ) %>%
  group_by(country) %>%
  arrange(country,date) %>%
  mutate(increase = cases - lag(cases, default =  first(cases))) -> df

ggplot(df, aes(x=date, y=increase, color=country  ))+
  geom_smooth() + 
  geom_label_repel(aes(label = paste(country,increase)), 
                   data = df[df$date == max(df$date), ], hjust = 1) + 
  labs(x=NULL,y=NULL)+ 
  theme_bw() + theme(legend.position = 'none') 


# 计算进度新冠确诊、治愈的趋势

library('tidyr')
country<-"India"
india<-hist[country]
india <- gather(india, curve, count, -date, -country)

ggplot(india, aes(date, count, color = curve)) + geom_point() + geom_line() + 
  labs(x=NULL,y=NULL,title=paste("Trend of cases, recovered and deaths in", country)) +
  scale_color_manual(values=c("#f39c12", "#dd4b39", "#00a65a")) +
  theme_bw() +   
  geom_label_repel(aes(label = paste(curve,count)), 
                   data = india[india$date == max(india$date), ], hjust = 1) + 
  theme(legend.position = "none", axis.text = element_text(angle = 15, hjust = 1)) +
  scale_x_date(date_labels = "%Y-%m-%d")


# heatmap

all <- hist["global"]
all <- all[all$cases > 0,]
length(unique(all$country))

all <- subset(all,date <= as.Date("2020-3-19"))
max_time <- max(all$date)
min_time <- max_time - 7
all <-  all[all$date >= min_time,]
all2 <- all[all$date == max(all$date,na.rm = TRUE),]

all$country <- factor(all$country, levels=unique(all2$country[order(all2$cases)]))
breaks = c(0,1000, 1000*10, 1000*100, 1000*1000)

ggplot(all, aes(date, country)) + 
  geom_tile(aes(fill = cases), color = 'black') + 
  scale_fill_viridis_c(trans = 'log',breaks = breaks, labels = breaks) + 
  xlab(NULL) + ylab(NULL) +
  scale_x_date(date_labels = "%Y-%m-%d") + theme_minimal()


# shiny 控制台

dashboard()