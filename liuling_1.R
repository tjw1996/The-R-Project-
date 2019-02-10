#### 载入数据和包 ####
library(tidyverse)
library(ade4)
library(vegan)
library(PerformanceAnalytics)
library(gclus)
library(cluster)
library(RColorBrewer)
df <- read.csv("D:\\Users\\Administrator\\Desktop\\DR.csv")
df <- arrange(df,df$Altitude)

### 数据转化 ###
scale_1 <- function(x) (x-min(x))/(max(x)-min(x))  # 0-1标准化函数

env <- df[,3:13]                                    # 环境矩阵
env.sc <- scale(env)                                # 标准化
env.sc_1 <- apply(env,MARGIN = 2,scale_1)           # 0-1标准化
env.sc.norm <- decostand(env.sc,"nor")              # 基于标准化的弦转化

cha <- df[,19:42]                                     # 性状矩阵
cha.hel <- decostand(cha,"hellinger")                 # hellinger转化
cha.chi <- decostand(cha,"chi.square")                # 卡方转化
cha.sc <- scale(cha)                                  # 标准化
cha.sc.norm <- decostand(cha.sc,"nor")                # 基于标准化的弦转化
cha.sc_1 <- apply(cha,MARGIN = 2,scale_1)             # 0-1标准化
cha.sc_1.norm <- decostand(cha.sc_1,"nor")            # 基于0-1标准化的弦转化
cha.sc_1.hel <- decostand(cha.sc_1,"hellinger")       # 基于0-1标准化的hellinger转化
cha.sc_1.chi <- decostand(cha.sc_1,"chi.square")      # 基于0-1标准化的卡方转化

#### 关联矩阵 ####
### 关联矩阵（Q模式） ###
## 环境矩阵 ##
source("coldiss.R")               # 绘图函数
env.de <- dist(env.sc)            # 欧式距离矩阵
coldiss(env.de,diag = T)
env.dn <- dist(env.sc.norm)       # 弦距离矩阵
coldiss(env.dn,diag = T)

## 性状矩阵 ##
cha.dh <- dist(cha.hel)              # 基于原始数据的hellinger距离矩阵
coldiss(cha.dh,diag = T)
cha.sc.de <- dist(cha.sc)            # 基于标准化的欧式距离矩阵
coldiss(cha.sc.de,diag = T)
cha.sc_1.de <- dist(cha.sc_1)         # 基于0-1标准化的欧式距离矩阵
coldiss(cha.sc_1.de,diag = T) 
cha.sc_1.dn <- dist(cha.sc_1.norm)    # 基于0-1标准化的弦距离矩阵
coldiss(cha.sc_1.dn,diag = T) 
cha.sc_1.dh <- dist(cha.sc_1.hel)     # 基于0-1标准化的hellinger距离矩阵
coldiss(cha.sc_1.dh,diag = T)
cha.sc_1.db <- vegdist(cha.sc_1)      # 基于0-1标准化的Bray-Curtis相异矩阵
coldiss(cha.sc_1.db,diag = T)

### 关联矩阵（R模式） ###
## 环境矩阵 ##
env.sc_1.t <- t(env.sc_1)               # 基于0-1标准化的卡方距离 
env.sc_1.t.chi <- decostand(env.sc_1.t,"chi.square")
env.sc_1.t.dc <- dist(env.sc_1.t.chi)   # 基于0-1标准化的卡方距离矩阵 
coldiss(env.sc_1.t.dc,diag = T)

## 性状矩阵 ##
cha.sc_1.t <- t(cha.sc_1)            # 基于0-1标准化的卡方距离 
cha.sc_1.t.chi <- decostand(cha.sc_1.t,"chi.square") 
cha.sc_1.t.dc <- dist(cha.sc_1.t.chi)   # 基于0-1标准化的卡方距离矩阵 
coldiss(cha.sc_1.t.dc,diag = T)

#### 基于性状变量对样方进行聚类 ####
## 基于0-1标准化的弦距离 ##
cha.sc_1.dn.single <- hclust(cha.sc_1.dn,method = "single")  # 单连接矩阵
plot(cha.sc_1.dn.single)
cha.sc_1.dn.complete <- hclust(cha.sc_1.dn,method = "complete") # 完全连接
plot(cha.sc_1.dn.complete)
cha.sc_1.dn.UPGMA <- hclust(cha.sc_1.dn,method = "average") # UPGMA连接
plot(cha.sc_1.dn.UPGMA)
cha.sc_1.dn.ward <- hclust(cha.sc_1.dn,method = "ward") # ward最小方差
plot(cha.sc_1.dn.ward)

## 同表型相关 ##
cha.sc_1.dn.single.coph <- cophenetic(cha.sc_1.dn.single)
cor(cha.sc_1.dn,cha.sc_1.dn.single.coph)
cha.sc_1.dn.complete.coph <- cophenetic(cha.sc_1.dn.complete)
cor(cha.sc_1.dn,cha.sc_1.dn.complete.coph)
cha.sc_1.dn.UPGMA.coph <- cophenetic(cha.sc_1.dn.UPGMA)
cor(cha.sc_1.dn,cha.sc_1.dn.UPGMA.coph)
cha.sc_1.dn.ward.coph <- cophenetic(cha.sc_1.dn.ward)
cor(cha.sc_1.dn,cha.sc_1.dn.ward.coph)
cor(cha.sc_1.dn,cha.sc_1.dn.single.coph,method = "spearman")
cor(cha.sc_1.dn,cha.sc_1.dn.complete.coph,method = "spearman")
cor(cha.sc_1.dn,cha.sc_1.dn.UPGMA.coph,method = "spearman")
cor(cha.sc_1.dn,cha.sc_1.dn.ward.coph,method = "spearman")
# Gower（1983）距离
gow.dist.single <- sum((cha.sc_1.dn-cha.sc_1.dn.single.coph)^2)
gow.dist.comp <- sum((cha.sc_1.dn-cha.sc_1.dn.complete.coph)^2)
gow.dist.UPGMA <- sum((cha.sc_1.dn-cha.sc_1.dn.UPGMA.coph)^2)
gow.dist.ward <- sum((cha.sc_1.dn-cha.sc_1.dn.ward.coph)^2)
gow.dist.single
gow.dist.comp
gow.dist.UPGMA
gow.dist.ward

## 寻找聚类簇 ##
par(mfrow=c(2,2))
plot(cha.sc_1.dn.single$height, nrow(cha):2, type="S", main="融合水平值-弦距离-单连接", ylab="k （聚类簇数量）", xlab="h (节点高度）", col="grey")
text(cha.sc_1.dn.single$height, nrow(cha):2, nrow(cha):2, col="red", cex=0.8)
plot(cha.sc_1.dn.complete$height, nrow(cha):2, type="S", main="融合水平值-弦距离-完全连接", ylab="k （聚类簇数量）", xlab="h (节点高度）", col="grey")
text(cha.sc_1.dn.complete$height, nrow(cha):2, nrow(cha):2, col="red", cex=0.8)
plot(cha.sc_1.dn.UPGMA$height, nrow(cha):2, type="S", main="融合水平值-弦距离-UPGMA", ylab="k （聚类簇数量）", xlab="h (节点高度）", col="grey")
text(cha.sc_1.dn.UPGMA$height, nrow(cha):2, nrow(cha):2, col="red", cex=0.8)
plot(cha.sc_1.dn.ward$height, nrow(cha):2, type="S", main="融合水平值-弦距离-Ward", ylab="k （聚类簇数量）", xlab="h (节点高度）", col="grey")
text(cha.sc_1.dn.ward$height, nrow(cha):2, nrow(cha):2, col="red", cex=0.8)

## 轮廓宽度图 ##
par(mfrow = c(1,2))
# 首先产生一个长度等于样方数量的空向量asw
asw <- numeric(nrow(cha))
# 其次循环获得轮廓宽度值并依次填入asw向量
for (k in 2:(nrow(cha)-1)) {
  sil <- silhouette(cutree(cha.sc_1.dn.complete, k=k), cha.sc_1.dn)
  asw[k] <- summary(sil)$avg.width
}

# 选择最佳（最大）轮廓宽度值
k.best <- which.max(asw)
# 利用cluster程序包内函数plot.silhouette（）绘制轮廓宽度值k
plot(1:nrow(cha), asw, type="h", 
     main="轮廓宽度-最优聚类簇数（Ward聚类）", 
     xlab="k (组数）", ylab="平均轮廓宽度")
axis(1, k.best, paste("最优",k.best,sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
    "with an average silhouette width of", max(asw), "\n")
# 最终分组的轮廓图
# 选择聚类簇的数量
k <- 4
# 轮廓图
cutg <- cutree(cha.sc_1.dn.complete, k=k)
sil <- silhouette(cutg, cha.sc_1.dn)
silo <- sortSilhouette(sil)
rownames(silo) <- row.names(cha)[attr(silo,"iOrd")]
plot(silo, main="轮廓宽度图-弦距离（complete聚类）", 
     cex.names=0.8, col=cutg+1, nmax.lab=100,border="white"
)
# 最终聚类树 
par(mfrow = c(1,1))
# 根据聚类结果重排距离矩阵
cha.scco <- reorder.hclust(cha.sc_1.dn.complete, cha.sc_1.dn)
# 绘制聚类树
source("hcoplot.R")	       
hcoplot(cha.sc_1.dn.complete, cha.sc_1.dn, k=k)
# 用聚类结果重排距离矩阵的热图
dend <- as.dendrogram(cha.scco)
heatmap(as.matrix(cha.sc_1.dn), Rowv=dend, symm=TRUE, margin=c(3,3))
# 基于聚类树的双排列群落表格的热图
heatmap(t(cha.sc_1), Rowv=NA, Colv=dend,
        col=c("white", brewer.pal(5,"Greens")), scale="none", margin=c(4,4), 
        ylab="性状", xlab="样方")

#### 基于性状变量对样方进行排序 ####
## PCA ##
cha.sc_1.pca <- rda(cha.sc_1.hel)
summary(cha.sc_1.pca)
summary(cha.sc_1.pca, scaling=1)
# 使用biplot（）函数绘制排序图 
par(mfrow=c(1,2))
biplot(cha.sc_1.pca, scaling=1, main="PCA-1型标尺")
biplot(cha.sc_1.pca, main="PCA-2型标尺")  # 默认 scaling = 2
# 使用cleanplot.pca（）函数绘图
source("cleanplot.pca.R")
cleanplot.pca(cha.sc_1.pca)              

## CA ## 
cha.sc_1.ca <- cca(cha.sc_1)
summary(cha.sc_1.ca)
plot(cha.sc_1.ca)

#### 基于环境变量对样方进行聚类 ####
env.dn.single <- hclust(env.dn,method = "single")  # 单连接矩阵
plot(env.dn.single)
env.dn.complete <- hclust(env.dn,method = "complete") # 完全连接
plot(cha.sc_1.dn.complete)
env.dn.UPGMA <- hclust(env.dn,method = "average") # UPGMA连接
plot(env.dn.UPGMA)
env.dn.ward <- hclust(env.dn,method = "ward") # ward最小方差
plot(env.dn.ward)

## 同表型相关 ##
env.dn.single.coph <- cophenetic(env.dn.single)
cor(env.dn,env.dn.single.coph)
env.dn.complete.coph <- cophenetic(env.dn.complete)
cor(env.dn,env.dn.complete.coph)
env.dn.UPGMA.coph <- cophenetic(env.dn.UPGMA)
cor(env.dn,env.dn.UPGMA.coph)
env.dn.ward.coph <- cophenetic(env.dn.ward)
cor(env.dn,env.dn.ward.coph)
cor(env.dn,env.dn.single.coph,method = "spearman")
cor(env.dn,env.dn.complete.coph,method = "spearman")
cor(env.dn,env.dn.UPGMA.coph,method = "spearman")
cor(env.dn,env.dn.ward.coph,method = "spearman")

## 轮廓宽度图 ##
par(mfrow = c(1,2))
# 首先产生一个长度等于样方数量的空向量asw
asw <- numeric(nrow(env))
# 其次循环获得轮廓宽度值并依次填入asw向量
for (k in 2:(nrow(env)-1)) {
  sil <- silhouette(cutree(env.dn.ward, k=k), env.dn)
  asw[k] <- summary(sil)$avg.width
}

# 选择最佳（最大）轮廓宽度值
k.best <- which.max(asw)
# 利用cluster程序包内函数plot.silhouette（）绘制轮廓宽度值k
plot(1:nrow(env), asw, type="h", 
     main="轮廓宽度-最优聚类簇数（Ward聚类）", 
     xlab="k (组数）", ylab="平均轮廓宽度")
axis(1, k.best, paste("最优",k.best,sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
    "with an average silhouette width of", max(asw), "\n")

# 最终分组的轮廓图
# 选择聚类簇的数量
k <- 3
# 轮廓图
cutg <- cutree(env.dn.ward, k=k)
sil <- silhouette(cutg,env.dn)
silo <- sortSilhouette(sil)
rownames(silo) <- row.names(env)[attr(silo,"iOrd")]
plot(silo, main="轮廓宽度图-弦距离（ward聚类）", 
     cex.names=0.8, col=cutg+1, nmax.lab=100,border="white"
)
# 最终聚类树 
par(mfrow = c(1,1))
# 根据聚类结果重排距离矩阵
env.dnwo <- reorder.hclust(env.dn.ward, env.dn)
# 绘制聚类树
source("hcoplot.R")	       
hcoplot(env.dn.ward, env.dn, k=k)
# 用聚类结果重排距离矩阵的热图
dend <- as.dendrogram(env.dnwo)
heatmap(as.matrix(env.dn), Rowv=dend, symm=TRUE, margin=c(3,3))
# 基于聚类树的双排列群落表格的热图
heatmap(t(env.sc), Rowv=NA, Colv=dend,
        col=c("white", brewer.pal(5,"Greens")), scale="none", margin=c(4,4), 
        ylab="环境", xlab="样方")

#### 基于环境变量对样方进行排序 ####
env2 <- cbind(env.sc,as.factor(df$Substrate))
colnames(env2)[12] <- "Substrate"
rownames(env2) <- seq(1,38)
env2 <- as.data.frame(env2)

## PCA ##
env.pca <- rda(env2)
summary(env.pca)
summary(env.pca, scaling=1)
# 使用biplot（）函数绘制排序图 
par(mfrow=c(1,2))
biplot(env.pca, scaling=1, main="PCA-1型标尺")
biplot(env.pca, main="PCA-2型标尺")  # 默认 scaling = 2

#### 约束排序 ####
### RDA ###
cha.rda <- rda(cha.sc_1.hel~.,env2) 
summary(cha.rda) 
# 从rda的结果中提取未校正R2 
(R2 <- RsquareAdj(cha.rda)$r.squared)
# 从rda的结果中提取校正R2
(R2adj <- RsquareAdj(cha.rda)$adj.r.squared)

# 1型标尺：距离三序图
plot(cha.rda, scaling=1, main="RDA三序图：cha.sc_1.hel～env2 - 1型标尺- 加权和样方坐标")
#此排序图同时显示所有的元素：样方、物种、定量解释变量（用箭头表示）
#和因子变量的形心。为了与定量解释变量区分，物种用不带箭头的线表示。
cha1.sc <- scores(cha.rda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, cha1.sc[, 1], cha1.sc[, 2], length=0, lty=1, col="red")
# 2型标尺（默认）：相关三序图
plot(cha.rda, main="RDA三序图：cha.sc_1.hel～env2 - 2型标尺- 加权和样方坐标")
cha2.sc <- scores(cha.rda, choices=1:2, display="sp")
arrows(0, 0, cha2.sc[, 1], cha2.sc[, 2], length=0, lty=1, col="red")
# 样方坐标是环境因子线性组合 
# 1型标尺
plot(cha.rda, scaling=1, display=c("sp", "lc", "cn"), 
     main="RDA三序图：cha.sc_1.hel～env2 - 1型标尺- 拟合的样方坐标")
arrows(0, 0, cha1.sc[, 1], cha1.sc[, 2], length=0, lty=1, col="red")
# 2型标尺
plot(cha.rda, display=c("sp", "lc", "cn"), 
     main="RDA三序图：cha.sc_1.hel～env2 - 2型标尺- 拟合的样方坐标")
arrows(0, 0, cha2.sc[,1], cha2.sc[,2], length=0, lty=1, col="red")
# RDA所有轴置换检验
anova.cca(cha.rda, step=1000)	
# 每个典范轴逐一检验
anova.cca(cha.rda, by="axis", step=1000)
# 使用vegan包内ordistep（）函数前向选择解释变量
cha.rda.all <- rda(cha.sc_1.hel ~ ., data=env2)
vif.cca(cha.rda.all)
step.forward <- ordistep(rda(cha.sc_1.hel ~ 1, data=env2), scope = formula(cha.rda.all ), 
                         direction="forward", pstep = 1000)
#简约的RDA分析
cha.rda.pars <- rda(cha.sc_1.hel ~ Temp + aspect + k + Altitude, data=env2)
cha.rda.pars
anova.cca(cha.rda.pars, step=1000)
anova.cca(cha.rda.pars, step=1000, by="axis")
vif.cca(cha.rda.pars)
RsquareAdj(cha.rda.all)$adj.r.squared
RsquareAdj(cha.rda.pars)$adj.r.squared

## 逆向RDA ##
cha2 <- as.data.frame(cha.sc_1)
env.rda <- rda(env.sc.norm~.,cha2) 
summary(env.rda) 
RsquareAdj(env.rda)$r.squared
RsquareAdj(env.rda)$adj.r.squared
# RDA所有轴置换检验
anova.cca(env.rda, step=1000)	
# 每个典范轴逐一检验
anova.cca(env.rda, by="axis", step=1000)
# 使用vegan包内ordistep（）函数前向选择解释变量
vif.cca(env.rda)
step.forward <- ordistep(rda(env.sc.norm ~ 1, data=cha2), scope = formula(env.rda ), 
                         direction="forward", pstep = 1000)
#简约的RDA分析
env.rda.pars <- rda( env.sc.norm ~ cell + d_cell + t_cell + s_cell + 
                      w_leaf + s_leaf + l_leaf + 厚z + 厚b + 长z + 长q + 宽q + 
                       宽z + d_paraxial_papilla + d_abaxial_papilla, data=cha2)
env.rda.pars
anova.cca(env.rda.pars, step=1000)
anova.cca(env.rda.pars, step=1000, by="axis")
vif.cca(env.rda.pars)
RsquareAdj(env.rda)$adj.r.squared
RsquareAdj(env.rda.pars)$adj.r.squared

## 修正数据RDA ##
cha.1 <- select(as.data.frame(cha.sc_1.hel),cell,d_cell,t_cell,s_cell,w_leaf,s_leaf,
                               l_leaf,厚z,厚b,长z,长q,宽q,宽z,
                               d_paraxial_papilla,d_abaxial_papilla)
### RDA ###
cha.rda <- rda(cha.1~.,env2) 
summary(cha.rda) 
# 从rda的结果中提取未校正R2 
(R2 <- RsquareAdj(cha.rda)$r.squared)
# 从rda的结果中提取校正R2
(R2adj <- RsquareAdj(cha.rda)$adj.r.squared)
# 1型标尺：距离三序图
plot(cha.rda, scaling=1, main="RDA三序图：cha.1～env2 - 1型标尺- 加权和样方坐标")
#此排序图同时显示所有的元素：样方、物种、定量解释变量（用箭头表示）
#和因子变量的形心。为了与定量解释变量区分，物种用不带箭头的线表示。
cha1.sc <- scores(cha.rda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, cha1.sc[, 1], cha1.sc[, 2], length=0, lty=1, col="red")
# 2型标尺（默认）：相关三序图
plot(cha.rda, main="RDA三序图：cha.1～env2 - 2型标尺- 加权和样方坐标")
cha2.sc <- scores(cha.rda, choices=1:2, display="sp")
arrows(0, 0, cha2.sc[, 1], cha2.sc[, 2], length=0, lty=1, col="red")
# RDA所有轴置换检验
anova.cca(cha.rda, step=1000)	
# 每个典范轴逐一检验
anova.cca(cha.rda, by="axis", step=1000)
# 使用vegan包内ordistep（）函数前向选择解释变量
step.forward <- ordistep(rda(cha.1 ~ 1, data=env2), scope = formula(cha.rda), 
                         direction="forward", pstep = 1000)
step.both <- ordistep(rda(cha.1 ~ 1, data=env2), scope = formula(cha.rda), 
                         direction="both", pstep = 1000)
step.backward <- ordistep(rda(cha.1 ~ ., data=env2), scope = formula(cha.rda), 
                      direction="backward", pstep = 1000)
#简约的RDA分析
cha.rda.pars <- rda(cha.1 ~ Temp + aspect + k + Altitude, data=env2)
cha.rda.pars
anova.cca(cha.rda.pars, step=1000)
anova.cca(cha.rda.pars, step=1000, by="axis")
vif.cca(cha.rda.pars)
RsquareAdj(cha.rda.pars)$r.squared
RsquareAdj(cha.rda.pars)$adj.r.squared

## 偏相关 ##
cha.prda <- rda(cha.1 ~ Temp + Precipitation + NDVI + k + ai + pet +
                     Condition(Latitude + Longitude + Altitude + aspect +
                                 slope + Substrate),data=env2)
summary(cha.prda) 
# 从rda的结果中提取未校正R2 
(R2 <- RsquareAdj(cha.prda)$r.squared)
# 从rda的结果中提取校正R2
(R2adj <- RsquareAdj(cha.prda)$adj.r.squared)
# 1型标尺：距离三序图
plot(cha.prda, scaling=1, main="偏RDA三序图：cha.1～env2 - 1型标尺- 加权和样方坐标")
#此排序图同时显示所有的元素：样方、物种、定量解释变量（用箭头表示）
#和因子变量的形心。为了与定量解释变量区分，物种用不带箭头的线表示。
cha1.sc <- scores(cha.prda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, cha1.sc[, 1], cha1.sc[, 2], length=0, lty=1, col="red")
# 2型标尺（默认）：相关三序图
plot(cha.prda, main="偏RDA三序图：cha.1～env2 - 2型标尺- 加权和样方坐标")
cha2.sc <- scores(cha.prda, choices=1:2, display="sp")
arrows(0, 0, cha2.sc[, 1], cha2.sc[, 2], length=0, lty=1, col="red")
