#### ��ǰ�� ####
#���ݵ���

spec <- read.csv("11-14-R.csv", header=TRUE) 
d <- spec[,2:length(spec)]
row.names(d) <- spec[,1]

#��ͼ1

library(tidyverse)
ggplot(data = d,mapping = aes( x = Altitude, y = s_leaf,group = Substrate)) +
  geom_point(mapping = aes(shape = Substrate))+
  geom_smooth(method = lm,se = FALSE)
ggsave(filename="al-ls.jpeg")
mapping = aes(group = matrix)





#��ͼ2��lowess��

library(gplots)
d<-na.omit(read.csv("C:\\Users\\1\\Desktop\\tx.csv"))
jpeg(filename="figtx1.jpg", width=12, height=12, units="in", res=600)
op<-par(mfrow=c(3,3), oma=c(1,1,1,1), mar=c(4,4,1,1),mgp=c(2.5,0.5,0))
for(i in 1:9){
  z1<-lowess(d$al, d[,8+i])
  plot(z1, xlab="Altitude")
}
par(op)
dev.off()



# Hier.part 



x <- d[,2:8]
y <- d[,9:22]
z <- list(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
c <- list(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
test <- list(1,2,3,4,5,6,7,8,9,10,11,12,13,14)

library(hier.part)

for(i in 1:14){ z[[i]] <-all.regs(y[,i], x, family = "gaussian", gof = "RMSPE", print.vars = FALSE)}
for(i in 1:14){ c[[i]] <-  partition(z[[i]], pcan = 7)}
for(i in 1:14){ test[[i]] <-rand.hp(y[,i], x, family = "gaussian", gof = "RMSPE", num.reps = 100)}

c[[8]]
test[[11]]

#��׼��
s <- matrix(nrow= ,ncol= )
for(i in 1:){s[,i] <- scale(, center=T,scale=T) 

#### ���ϵ��������ӻ� ####
### �������ݺͰ� ### 
rm(list = ls())
library(tidyverse)
# library(vegan)
library(psych)
library(Hmisc)
library(corrplot)               #  install.packages("corrplot")
library(PerformanceAnalytics)   # install.packages("PerformanceAnalytics")
library('lavaan')
library("DiagrammeR")
library('nloptr')
library("semPlot")
df_1 <- read.csv("D:\\Users\\Administrator\\Desktop\\DR.csv")
df_2 <- df_1[,3:42]
### ȫ������ ###
res1 <- rcorr(as.matrix(df_2),type="pearson")                      
res1_r <- res1$r
res1_p <- res1$P                                    # ����Pֵ��������
# corrplot����ͼ��
corrplot(res1_r,type = "upper",tl.col = "black",
         tl.srt = 45,p.mat = res1_p, sig.level = 0.05, 
         insig = "blank")               #����������ؾ�����ӻ�
### �������� ###
df_3 <- df_2 %>% select(-c("SR_mean","SR_spr","SR_sum","SR_aut",            
                            "SR_win","ai","d_abaxial_papilla","s_abaxial_papilla", 
                           "d_paraxial_papilla","s_paraxial_papilla","h_abaxial_papilla", 
                           "h_paraxial_papilla"))
res2 <- rcorr(as.matrix(df_3),type="pearson")                      
res2_r <- res2$r
res2_p <- res2$P                                    # ����Pֵ��������
# corrplot����ͼ��
corrplot(res2_r,type = "upper",tl.col = "black",
         tl.srt = 45,p.mat = res2_p, sig.level = 0.05, 
         insig = "blank")               # ����������ؾ�����ӻ�

#### CCA���� ####

### ��ֱ��� ###
sc <- df_1[,24:42]                     # ϸ����״����
rownames(sc) <- seq(1,38)
sc <- as.data.frame(scale(sc))
sl <- df_1[,c(19:23)]                  # ҶƬ��״����
rownames(sl) <- seq(1,38)
sl <- as.data.frame(scale(sl))
se <- df_1[,3:13]                       # �������ӱ���
rownames(se) <- seq(1,38)
se <- as.data.frame(scale(se))
se_1 <- cbind(se,df_1[,2])

### ҶƬ��״~�������� ###
sp1_1 <- rda(sl ~ . , se)
step_forward <- ordistep(rda(sl ~ 1,se),scope = formula(sp1_1),
                         direction = "forward",pstep = 1000)
summary(sp1_1)
plot(sp1_1)
# �������ͱ��� #
sp1_2 <- rda(sl ~  Altitude + 
               Temp + pet , se) 
summary(sp1_2)
plot(sp1_2)

### ϸ����״~�������� ###
sp2_1 <- rda(sc ~ ., se)
summary(sp2_1)
plot(sp2_1)
# �������ͱ���
sp2_2 <- rda(sc~Latitude + Altitude + Temp + 
               Precipitation + Longitude + pet + k,se)
summary(sp2_2)
plot(sp2_2)

## �� ##
sp2_3 <- rda(sc[,14:19]~.,se)
summary(sp2_3)
plot(sp2_3)

## ϸ����̬ ##
sp2_4 <- rda(sc[,1:13]~.,se)
summary(sp2_4)
plot(sp2_4)

### ҶƬ��״~ϸ����״ ###
sp3_1 <- rda(sl ~ . ,sc) 
summary(sp3_1)
plot(sp3_1)

## ϸ����̬ ##
sp3_4 <- rda(sl~.,sc[,1:13])
summary(sp3_4)
plot(sp3_4)


#### �������ӵ����ӷ��� ####
library(psych)
library(GPArotation)     # install.packages("GPArotation")
se_cor <- cor(se)
fa.parallel(se_cor,n.obs = 38,fa = "both",n.iter = 100)

## �������� ##
fa_2 <- fa(se_cor,nfactors = 2,rotate = "none",fm = "ml")
fa_2
summary(fa_2)
factor.plot(fa_2,labels = rownames(fa_2$loadings))
fa.diagram(fa_2,sample = FALSE)

## б����ת ##
fa_2_promax <- fa(se_cor,nfactors = 2,rotate = "promax",fm = "pa")
fa_2_promax
fa.diagram(fa_2_promax,sample = FALSE)

#### ҶƬ��״�����ӷ��� ####
sl_cor <- cor(sl)
fa.parallel(sl_cor,n.obs = 38,fa = "both",n.iter = 100)

## һ������ ##
fa_2 <- fa(sl_cor,nfactors = 1,rotate = "none",fm = "pa")
fa_2
summary(fa_2)
factor.plot(fa_2,labels = rownames(fa_2$loadings))
fa.diagram(fa_2,sample = FALSE)

#### ϸ����״�����ӷ��� ####

### ϸ������ ###
sc_cor <- cor(sc[,10:19])
fa.parallel(sc_cor,n.obs = 38,fa = "both",n.iter = 100)

## �������� ##
fa_2 <- fa(sc_cor,nfactors = 2,rotate = "none",fm = "ml")
fa_2
summary(fa_2)
factor.plot(fa_2,labels = rownames(fa_2$loadings))
fa.diagram(fa_2,sample = FALSE)

## б����ת ##
fa_2_promax <- fa(sc_cor,nfactors = 2,rotate = "promax",fm = "pa")
fa_2_promax
fa.diagram(fa_2_promax,sample = FALSE)

### ������ ###
sc_cor_1 <- cor(sc[,1:9])
fa.parallel(sc_cor_1,n.obs = 38,fa = "both",n.iter = 100)

## һ������ ##
fa_2_1 <- fa(sc_cor_1,nfactors = 1,rotate = "none",fm = "pa")
fa_2_1
summary(fa_2_1)
factor.plot(fa_2_1,labels = rownames(fa_2$loadings))
fa.diagram(fa_2_1,sample = FALSE)

### ϸ��������״ ###
sc_cor <- cor(sc)
fa.parallel(sc_cor,n.obs = 38,fa = "both",n.iter = 100)

## �������� ##
fa_2 <- fa(sc_cor,nfactors = 3,rotate = "none",fm = "pa")
fa_2
summary(fa_2)
factor.plot(fa_2,labels = rownames(fa_2$loadings))
fa.diagram(fa_2,sample = FALSE)

## б����ת ##
fa_2_promax <- fa(sc_cor,nfactors = 2,rotate = "promax",fm = "pa")
fa_2_promax
fa.diagram(fa_2_promax,sample = FALSE)

#### �������ӷ���ģ�� ####
### �������Ӳ���ģ�� ###
model_1 <- ' 
PA_1 =~ Altitude + Latitude
PA_2 =~ Precipitation + k + pet
PA_3 =~ Temp
PA_3 ~ PA_1 + PA_2
PA_2 ~ PA_1
Temp ~~ pet
Precipitation ~~ k
pet ~~ k
'
model_1_1 <- ' 
PA_1 =~ Precipitation + k + pet
PA_2 =~ Temp
PA_2 ~~ PA_1
PA_1 ~ Latitude + Altitude
PA_2 ~ Latitude + Altitude
Temp ~~ pet
Precipitation ~~ k
pet ~~ k
'
model_2 <- ' 
PA_2 =~ Precipitation + k + pet
PA_3 =~ Temp
PA_3 ~~  PA_2
Temp ~~ pet
Precipitation ~~ k
'
model_3 <- ' 
PA_2 =~ Precipitation + k 
PA_3 =~ Temp + pet
PA_3 ~~  PA_2
Temp ~~ pet
Precipitation ~~ k
'
model_4 <- ' 
PA_1 =~ Altitude + Latitude
PA_2 =~ Precipitation + k 
PA_3 =~ Temp + pet
PA_3 ~ PA_1 
PA_2 ~ PA_1
PA_3 ~~ PA_2
Temp ~~ pet
Precipitation ~~ k
'
model_5 <- ' 
PA_1 =~ Longitude + Altitude + Latitude
PA_2 =~ Precipitation + k + pet + Temp
PA_2 ~ PA_1
Temp ~~ pet
Precipitation ~~ k
'
model_6 <- ' 
PA_1 =~ Longitude 
PA_2 =~ Altitude + Latitude
PA_3 =~ Precipitation + k + pet + Temp
PA_3 ~ PA_1
PA_3 ~ PA_2
Temp ~~ pet
Precipitation ~~ k
'
model_7 <- ' 
pet ~ Temp + Latitude + Altitude
Temp ~ Latitude + Altitude
'

model_8 <- ' 
pet ~  Latitude + Altitude
Temp ~ Latitude + Altitude
Temp ~~ pet
'

fit_1 <- lavaan::sem(model_1_1, data= se)
fit_2 <- lavaan::sem(model_2, data= se)
fit_3 <- lavaan::sem(model_3, data= se)
fit_4 <- lavaan::sem(model_4, data= se)
fit_5 <- lavaan::sem(model_5, data= se)
fit_6 <- lavaan::sem(model_6, data= se)
fit_7 <- lavaan::sem(model_7, data= se)
fit_8 <- lavaan::sem(model_8, data= se)

summary(fit_1, fit.measures=TRUE)
semPaths(fit_1, intercept = FALSE, whatLabel = "est",
         residuals = FALSE, exoCov = FALSE)

### sem�� ###
se_cov <- cov(se[,-c(4:5,8,10)])
model.se.1 <- specifyModel(text ="
PA_2 -> PA_3,2_3,NA
PA_1 -> PA_3,1_3,NA
PA_1 -> PA_2,1_2,NA
PA_1 -> Latitude,1_La,NA
PA_1 -> Longitude,1_Lo,NA
PA_1 -> Altitude,1_A,NA
PA_2 -> Precipitation,2_P,NA
PA_2 -> k,2_k,NA
PA_2 -> pet,2_p,NA
PA_3 -> Temp,3_T,NA
Precipitation <-> k,P_k,NA
pet <-> Temp,p_T,NA
")
sem.se.1 <- sem::sem(model.se.1,se_cov,266)
summary(sem.se.1)
semPaths(sem.se.1, intercept = FALSE, whatLabel = "est",
         residuals = FALSE, exoCov = FALSE)
#### ��������ͨ������ ####

model_1 <- ' 
Temp ~ Latitude + Altitude
pet ~ Latitude + Altitude
pet ~ Temp
Latitude ~~ Altitude
'
fit_1 <- lavaan::sem(model_1, data= se)
summary(fit_1, fit.measures=TRUE)
semPaths(fit_1, intercept = FALSE, whatLabel = "est",
         residuals = FALSE, exoCov = FALSE)

model_2 <- ' 
Precipitation ~ Longitude 
k ~  Longitude + Precipitation
'
fit_2 <- lavaan::sem(model_2, data= se)
summary(fit_2, fit.measures=TRUE)
semPaths(fit_2, intercept = FALSE, whatLabel = "est",
         residuals = FALSE, exoCov = FALSE)

model_3 <- ' 
Precipitation ~ Longitude 
k ~  Longitude + Precipitation
Temp ~ Latitude + Altitude
pet ~ Latitude + Altitude + Temp
pet ~ k
'
fit_3 <- lavaan::sem(model_3, data= se)
summary(fit_3, fit.measures=TRUE)
semPaths(fit_3, intercept = FALSE, whatLabel = "est",
         residuals = FALSE, exoCov = FALSE)

#### ���ݱȶ� ####
chart.Correlation(se, histogram=TRUE, pch=19)
chart.Correlation(PoliticalDemocracy, histogram=TRUE, pch=19)
res1 <- rcorr(as.matrix(se),type="pearson")                      
res1_r <- res1$r
res1_p <- res1$P                                    # ����Pֵ��������
# corrplot����ͼ��
corrplot(res1_r,type = "upper",tl.col = "black",
         tl.srt = 45,p.mat = res1_p, sig.level = 0.05, 
         insig = "blank") 
res1 <- rcorr(as.matrix(PoliticalDemocracy),type="pearson")                      
res1_r <- res1$r
res1_p <- res1$P                                    # ����Pֵ��������
# corrplot����ͼ��
corrplot(res1_r,type = "upper",tl.col = "black",
         tl.srt = 45,p.mat = res1_p, sig.level = 0.05, 
         insig = "blank") 
res1 <- rcorr(as.matrix(demoOneFactor),type="pearson")                      
res1_r <- res1$r
res1_p <- res1$P                                    # ����Pֵ��������
# corrplot����ͼ��
corrplot(res1_r,type = "upper",tl.col = "black",
         tl.srt = 45,p.mat = res1_p, sig.level = 0.05, 
         insig = "blank")
res1 <- rcorr(as.matrix(se[,c(1,3,6,11)]),type="pearson")                      
res1_r <- res1$r
res1_p <- res1$P                                    # ����Pֵ��������
# corrplot����ͼ��
corrplot(res1_r,type = "upper",tl.col = "black",
         tl.srt = 45,p.mat = res1_p, sig.level = 0.05, 
         insig = "blank")
chart.Correlation(se[,c(1,3,6,11)], histogram=TRUE, pch=19)