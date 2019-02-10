#### �������ݺͰ� ####
library(tidyverse)
library(ade4)
library(vegan)
library(PerformanceAnalytics)
library(gclus)
library(cluster)
library(RColorBrewer)
df <- read.csv("D:\\Users\\Administrator\\Desktop\\DR.csv")
df <- arrange(df,df$Altitude)

### ����ת�� ###
scale_1 <- function(x) (x-min(x))/(max(x)-min(x))  # 0-1��׼������

env <- df[,3:13]                                    # ��������
env.sc <- scale(env)                                # ��׼��
env.sc_1 <- apply(env,MARGIN = 2,scale_1)           # 0-1��׼��
env.sc.norm <- decostand(env.sc,"nor")              # ���ڱ�׼������ת��

cha <- df[,19:42]                                     # ��״����
cha.hel <- decostand(cha,"hellinger")                 # hellingerת��
cha.chi <- decostand(cha,"chi.square")                # ����ת��
cha.sc <- scale(cha)                                  # ��׼��
cha.sc.norm <- decostand(cha.sc,"nor")                # ���ڱ�׼������ת��
cha.sc_1 <- apply(cha,MARGIN = 2,scale_1)             # 0-1��׼��
cha.sc_1.norm <- decostand(cha.sc_1,"nor")            # ����0-1��׼������ת��
cha.sc_1.hel <- decostand(cha.sc_1,"hellinger")       # ����0-1��׼����hellingerת��
cha.sc_1.chi <- decostand(cha.sc_1,"chi.square")      # ����0-1��׼���Ŀ���ת��

#### �������� ####
### ��������Qģʽ�� ###
## �������� ##
source("coldiss.R")               # ��ͼ����
env.de <- dist(env.sc)            # ŷʽ�������
coldiss(env.de,diag = T)
env.dn <- dist(env.sc.norm)       # �Ҿ������
coldiss(env.dn,diag = T)

## ��״���� ##
cha.dh <- dist(cha.hel)              # ����ԭʼ���ݵ�hellinger�������
coldiss(cha.dh,diag = T)
cha.sc.de <- dist(cha.sc)            # ���ڱ�׼����ŷʽ�������
coldiss(cha.sc.de,diag = T)
cha.sc_1.de <- dist(cha.sc_1)         # ����0-1��׼����ŷʽ�������
coldiss(cha.sc_1.de,diag = T) 
cha.sc_1.dn <- dist(cha.sc_1.norm)    # ����0-1��׼�����Ҿ������
coldiss(cha.sc_1.dn,diag = T) 
cha.sc_1.dh <- dist(cha.sc_1.hel)     # ����0-1��׼����hellinger�������
coldiss(cha.sc_1.dh,diag = T)
cha.sc_1.db <- vegdist(cha.sc_1)      # ����0-1��׼����Bray-Curtis�������
coldiss(cha.sc_1.db,diag = T)

### ��������Rģʽ�� ###
## �������� ##
env.sc_1.t <- t(env.sc_1)               # ����0-1��׼���Ŀ������� 
env.sc_1.t.chi <- decostand(env.sc_1.t,"chi.square")
env.sc_1.t.dc <- dist(env.sc_1.t.chi)   # ����0-1��׼���Ŀ���������� 
coldiss(env.sc_1.t.dc,diag = T)

## ��״���� ##
cha.sc_1.t <- t(cha.sc_1)            # ����0-1��׼���Ŀ������� 
cha.sc_1.t.chi <- decostand(cha.sc_1.t,"chi.square") 
cha.sc_1.t.dc <- dist(cha.sc_1.t.chi)   # ����0-1��׼���Ŀ���������� 
coldiss(cha.sc_1.t.dc,diag = T)

#### ������״�������������о��� ####
## ����0-1��׼�����Ҿ��� ##
cha.sc_1.dn.single <- hclust(cha.sc_1.dn,method = "single")  # �����Ӿ���
plot(cha.sc_1.dn.single)
cha.sc_1.dn.complete <- hclust(cha.sc_1.dn,method = "complete") # ��ȫ����
plot(cha.sc_1.dn.complete)
cha.sc_1.dn.UPGMA <- hclust(cha.sc_1.dn,method = "average") # UPGMA����
plot(cha.sc_1.dn.UPGMA)
cha.sc_1.dn.ward <- hclust(cha.sc_1.dn,method = "ward") # ward��С����
plot(cha.sc_1.dn.ward)

## ͬ������� ##
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
# Gower��1983������
gow.dist.single <- sum((cha.sc_1.dn-cha.sc_1.dn.single.coph)^2)
gow.dist.comp <- sum((cha.sc_1.dn-cha.sc_1.dn.complete.coph)^2)
gow.dist.UPGMA <- sum((cha.sc_1.dn-cha.sc_1.dn.UPGMA.coph)^2)
gow.dist.ward <- sum((cha.sc_1.dn-cha.sc_1.dn.ward.coph)^2)
gow.dist.single
gow.dist.comp
gow.dist.UPGMA
gow.dist.ward

## Ѱ�Ҿ���� ##
par(mfrow=c(2,2))
plot(cha.sc_1.dn.single$height, nrow(cha):2, type="S", main="�ں�ˮƽֵ-�Ҿ���-������", ylab="k �������������", xlab="h (�ڵ�߶ȣ�", col="grey")
text(cha.sc_1.dn.single$height, nrow(cha):2, nrow(cha):2, col="red", cex=0.8)
plot(cha.sc_1.dn.complete$height, nrow(cha):2, type="S", main="�ں�ˮƽֵ-�Ҿ���-��ȫ����", ylab="k �������������", xlab="h (�ڵ�߶ȣ�", col="grey")
text(cha.sc_1.dn.complete$height, nrow(cha):2, nrow(cha):2, col="red", cex=0.8)
plot(cha.sc_1.dn.UPGMA$height, nrow(cha):2, type="S", main="�ں�ˮƽֵ-�Ҿ���-UPGMA", ylab="k �������������", xlab="h (�ڵ�߶ȣ�", col="grey")
text(cha.sc_1.dn.UPGMA$height, nrow(cha):2, nrow(cha):2, col="red", cex=0.8)
plot(cha.sc_1.dn.ward$height, nrow(cha):2, type="S", main="�ں�ˮƽֵ-�Ҿ���-Ward", ylab="k �������������", xlab="h (�ڵ�߶ȣ�", col="grey")
text(cha.sc_1.dn.ward$height, nrow(cha):2, nrow(cha):2, col="red", cex=0.8)

## ��������ͼ ##
par(mfrow = c(1,2))
# ���Ȳ���һ�����ȵ������������Ŀ�����asw
asw <- numeric(nrow(cha))
# ���ѭ�������������ֵ����������asw����
for (k in 2:(nrow(cha)-1)) {
  sil <- silhouette(cutree(cha.sc_1.dn.complete, k=k), cha.sc_1.dn)
  asw[k] <- summary(sil)$avg.width
}

# ѡ����ѣ������������ֵ
k.best <- which.max(asw)
# ����cluster������ں���plot.silhouette����������������ֵk
plot(1:nrow(cha), asw, type="h", 
     main="��������-���ž��������Ward���ࣩ", 
     xlab="k (������", ylab="ƽ����������")
axis(1, k.best, paste("����",k.best,sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
    "with an average silhouette width of", max(asw), "\n")
# ���շ��������ͼ
# ѡ�����ص�����
k <- 4
# ����ͼ
cutg <- cutree(cha.sc_1.dn.complete, k=k)
sil <- silhouette(cutg, cha.sc_1.dn)
silo <- sortSilhouette(sil)
rownames(silo) <- row.names(cha)[attr(silo,"iOrd")]
plot(silo, main="��������ͼ-�Ҿ��루complete���ࣩ", 
     cex.names=0.8, col=cutg+1, nmax.lab=100,border="white"
)
# ���վ����� 
par(mfrow = c(1,1))
# ���ݾ��������ž������
cha.scco <- reorder.hclust(cha.sc_1.dn.complete, cha.sc_1.dn)
# ���ƾ�����
source("hcoplot.R")	       
hcoplot(cha.sc_1.dn.complete, cha.sc_1.dn, k=k)
# �þ��������ž���������ͼ
dend <- as.dendrogram(cha.scco)
heatmap(as.matrix(cha.sc_1.dn), Rowv=dend, symm=TRUE, margin=c(3,3))
# ���ھ�������˫����Ⱥ��������ͼ
heatmap(t(cha.sc_1), Rowv=NA, Colv=dend,
        col=c("white", brewer.pal(5,"Greens")), scale="none", margin=c(4,4), 
        ylab="��״", xlab="����")

#### ������״������������������ ####
## PCA ##
cha.sc_1.pca <- rda(cha.sc_1.hel)
summary(cha.sc_1.pca)
summary(cha.sc_1.pca, scaling=1)
# ʹ��biplot����������������ͼ 
par(mfrow=c(1,2))
biplot(cha.sc_1.pca, scaling=1, main="PCA-1�ͱ��")
biplot(cha.sc_1.pca, main="PCA-2�ͱ��")  # Ĭ�� scaling = 2
# ʹ��cleanplot.pca����������ͼ
source("cleanplot.pca.R")
cleanplot.pca(cha.sc_1.pca)              

## CA ## 
cha.sc_1.ca <- cca(cha.sc_1)
summary(cha.sc_1.ca)
plot(cha.sc_1.ca)

#### ���ڻ����������������о��� ####
env.dn.single <- hclust(env.dn,method = "single")  # �����Ӿ���
plot(env.dn.single)
env.dn.complete <- hclust(env.dn,method = "complete") # ��ȫ����
plot(cha.sc_1.dn.complete)
env.dn.UPGMA <- hclust(env.dn,method = "average") # UPGMA����
plot(env.dn.UPGMA)
env.dn.ward <- hclust(env.dn,method = "ward") # ward��С����
plot(env.dn.ward)

## ͬ������� ##
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

## ��������ͼ ##
par(mfrow = c(1,2))
# ���Ȳ���һ�����ȵ������������Ŀ�����asw
asw <- numeric(nrow(env))
# ���ѭ�������������ֵ����������asw����
for (k in 2:(nrow(env)-1)) {
  sil <- silhouette(cutree(env.dn.ward, k=k), env.dn)
  asw[k] <- summary(sil)$avg.width
}

# ѡ����ѣ������������ֵ
k.best <- which.max(asw)
# ����cluster������ں���plot.silhouette����������������ֵk
plot(1:nrow(env), asw, type="h", 
     main="��������-���ž��������Ward���ࣩ", 
     xlab="k (������", ylab="ƽ����������")
axis(1, k.best, paste("����",k.best,sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
    "with an average silhouette width of", max(asw), "\n")

# ���շ��������ͼ
# ѡ�����ص�����
k <- 3
# ����ͼ
cutg <- cutree(env.dn.ward, k=k)
sil <- silhouette(cutg,env.dn)
silo <- sortSilhouette(sil)
rownames(silo) <- row.names(env)[attr(silo,"iOrd")]
plot(silo, main="��������ͼ-�Ҿ��루ward���ࣩ", 
     cex.names=0.8, col=cutg+1, nmax.lab=100,border="white"
)
# ���վ����� 
par(mfrow = c(1,1))
# ���ݾ��������ž������
env.dnwo <- reorder.hclust(env.dn.ward, env.dn)
# ���ƾ�����
source("hcoplot.R")	       
hcoplot(env.dn.ward, env.dn, k=k)
# �þ��������ž���������ͼ
dend <- as.dendrogram(env.dnwo)
heatmap(as.matrix(env.dn), Rowv=dend, symm=TRUE, margin=c(3,3))
# ���ھ�������˫����Ⱥ��������ͼ
heatmap(t(env.sc), Rowv=NA, Colv=dend,
        col=c("white", brewer.pal(5,"Greens")), scale="none", margin=c(4,4), 
        ylab="����", xlab="����")

#### ���ڻ��������������������� ####
env2 <- cbind(env.sc,as.factor(df$Substrate))
colnames(env2)[12] <- "Substrate"
rownames(env2) <- seq(1,38)
env2 <- as.data.frame(env2)

## PCA ##
env.pca <- rda(env2)
summary(env.pca)
summary(env.pca, scaling=1)
# ʹ��biplot����������������ͼ 
par(mfrow=c(1,2))
biplot(env.pca, scaling=1, main="PCA-1�ͱ��")
biplot(env.pca, main="PCA-2�ͱ��")  # Ĭ�� scaling = 2

#### Լ������ ####
### RDA ###
cha.rda <- rda(cha.sc_1.hel~.,env2) 
summary(cha.rda) 
# ��rda�Ľ������ȡδУ��R2 
(R2 <- RsquareAdj(cha.rda)$r.squared)
# ��rda�Ľ������ȡУ��R2
(R2adj <- RsquareAdj(cha.rda)$adj.r.squared)

# 1�ͱ�ߣ���������ͼ
plot(cha.rda, scaling=1, main="RDA����ͼ��cha.sc_1.hel��env2 - 1�ͱ��- ��Ȩ����������")
#������ͼͬʱ��ʾ���е�Ԫ�أ����������֡��������ͱ������ü�ͷ��ʾ��
#�����ӱ��������ġ�Ϊ���붨�����ͱ������֣������ò�����ͷ���߱�ʾ��
cha1.sc <- scores(cha.rda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, cha1.sc[, 1], cha1.sc[, 2], length=0, lty=1, col="red")
# 2�ͱ�ߣ�Ĭ�ϣ����������ͼ
plot(cha.rda, main="RDA����ͼ��cha.sc_1.hel��env2 - 2�ͱ��- ��Ȩ����������")
cha2.sc <- scores(cha.rda, choices=1:2, display="sp")
arrows(0, 0, cha2.sc[, 1], cha2.sc[, 2], length=0, lty=1, col="red")
# ���������ǻ�������������� 
# 1�ͱ��
plot(cha.rda, scaling=1, display=c("sp", "lc", "cn"), 
     main="RDA����ͼ��cha.sc_1.hel��env2 - 1�ͱ��- ��ϵ���������")
arrows(0, 0, cha1.sc[, 1], cha1.sc[, 2], length=0, lty=1, col="red")
# 2�ͱ��
plot(cha.rda, display=c("sp", "lc", "cn"), 
     main="RDA����ͼ��cha.sc_1.hel��env2 - 2�ͱ��- ��ϵ���������")
arrows(0, 0, cha2.sc[,1], cha2.sc[,2], length=0, lty=1, col="red")
# RDA�������û�����
anova.cca(cha.rda, step=1000)	
# ÿ���䷶����һ����
anova.cca(cha.rda, by="axis", step=1000)
# ʹ��vegan����ordistep��������ǰ��ѡ����ͱ���
cha.rda.all <- rda(cha.sc_1.hel ~ ., data=env2)
vif.cca(cha.rda.all)
step.forward <- ordistep(rda(cha.sc_1.hel ~ 1, data=env2), scope = formula(cha.rda.all ), 
                         direction="forward", pstep = 1000)
#��Լ��RDA����
cha.rda.pars <- rda(cha.sc_1.hel ~ Temp + aspect + k + Altitude, data=env2)
cha.rda.pars
anova.cca(cha.rda.pars, step=1000)
anova.cca(cha.rda.pars, step=1000, by="axis")
vif.cca(cha.rda.pars)
RsquareAdj(cha.rda.all)$adj.r.squared
RsquareAdj(cha.rda.pars)$adj.r.squared

## ����RDA ##
cha2 <- as.data.frame(cha.sc_1)
env.rda <- rda(env.sc.norm~.,cha2) 
summary(env.rda) 
RsquareAdj(env.rda)$r.squared
RsquareAdj(env.rda)$adj.r.squared
# RDA�������û�����
anova.cca(env.rda, step=1000)	
# ÿ���䷶����һ����
anova.cca(env.rda, by="axis", step=1000)
# ʹ��vegan����ordistep��������ǰ��ѡ����ͱ���
vif.cca(env.rda)
step.forward <- ordistep(rda(env.sc.norm ~ 1, data=cha2), scope = formula(env.rda ), 
                         direction="forward", pstep = 1000)
#��Լ��RDA����
env.rda.pars <- rda( env.sc.norm ~ cell + d_cell + t_cell + s_cell + 
                      w_leaf + s_leaf + l_leaf + ��z + ��b + ��z + ��q + ��q + 
                       ��z + d_paraxial_papilla + d_abaxial_papilla, data=cha2)
env.rda.pars
anova.cca(env.rda.pars, step=1000)
anova.cca(env.rda.pars, step=1000, by="axis")
vif.cca(env.rda.pars)
RsquareAdj(env.rda)$adj.r.squared
RsquareAdj(env.rda.pars)$adj.r.squared

## ��������RDA ##
cha.1 <- select(as.data.frame(cha.sc_1.hel),cell,d_cell,t_cell,s_cell,w_leaf,s_leaf,
                               l_leaf,��z,��b,��z,��q,��q,��z,
                               d_paraxial_papilla,d_abaxial_papilla)
### RDA ###
cha.rda <- rda(cha.1~.,env2) 
summary(cha.rda) 
# ��rda�Ľ������ȡδУ��R2 
(R2 <- RsquareAdj(cha.rda)$r.squared)
# ��rda�Ľ������ȡУ��R2
(R2adj <- RsquareAdj(cha.rda)$adj.r.squared)
# 1�ͱ�ߣ���������ͼ
plot(cha.rda, scaling=1, main="RDA����ͼ��cha.1��env2 - 1�ͱ��- ��Ȩ����������")
#������ͼͬʱ��ʾ���е�Ԫ�أ����������֡��������ͱ������ü�ͷ��ʾ��
#�����ӱ��������ġ�Ϊ���붨�����ͱ������֣������ò�����ͷ���߱�ʾ��
cha1.sc <- scores(cha.rda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, cha1.sc[, 1], cha1.sc[, 2], length=0, lty=1, col="red")
# 2�ͱ�ߣ�Ĭ�ϣ����������ͼ
plot(cha.rda, main="RDA����ͼ��cha.1��env2 - 2�ͱ��- ��Ȩ����������")
cha2.sc <- scores(cha.rda, choices=1:2, display="sp")
arrows(0, 0, cha2.sc[, 1], cha2.sc[, 2], length=0, lty=1, col="red")
# RDA�������û�����
anova.cca(cha.rda, step=1000)	
# ÿ���䷶����һ����
anova.cca(cha.rda, by="axis", step=1000)
# ʹ��vegan����ordistep��������ǰ��ѡ����ͱ���
step.forward <- ordistep(rda(cha.1 ~ 1, data=env2), scope = formula(cha.rda), 
                         direction="forward", pstep = 1000)
step.both <- ordistep(rda(cha.1 ~ 1, data=env2), scope = formula(cha.rda), 
                         direction="both", pstep = 1000)
step.backward <- ordistep(rda(cha.1 ~ ., data=env2), scope = formula(cha.rda), 
                      direction="backward", pstep = 1000)
#��Լ��RDA����
cha.rda.pars <- rda(cha.1 ~ Temp + aspect + k + Altitude, data=env2)
cha.rda.pars
anova.cca(cha.rda.pars, step=1000)
anova.cca(cha.rda.pars, step=1000, by="axis")
vif.cca(cha.rda.pars)
RsquareAdj(cha.rda.pars)$r.squared
RsquareAdj(cha.rda.pars)$adj.r.squared

## ƫ��� ##
cha.prda <- rda(cha.1 ~ Temp + Precipitation + NDVI + k + ai + pet +
                     Condition(Latitude + Longitude + Altitude + aspect +
                                 slope + Substrate),data=env2)
summary(cha.prda) 
# ��rda�Ľ������ȡδУ��R2 
(R2 <- RsquareAdj(cha.prda)$r.squared)
# ��rda�Ľ������ȡУ��R2
(R2adj <- RsquareAdj(cha.prda)$adj.r.squared)
# 1�ͱ�ߣ���������ͼ
plot(cha.prda, scaling=1, main="ƫRDA����ͼ��cha.1��env2 - 1�ͱ��- ��Ȩ����������")
#������ͼͬʱ��ʾ���е�Ԫ�أ����������֡��������ͱ������ü�ͷ��ʾ��
#�����ӱ��������ġ�Ϊ���붨�����ͱ������֣������ò�����ͷ���߱�ʾ��
cha1.sc <- scores(cha.prda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, cha1.sc[, 1], cha1.sc[, 2], length=0, lty=1, col="red")
# 2�ͱ�ߣ�Ĭ�ϣ����������ͼ
plot(cha.prda, main="ƫRDA����ͼ��cha.1��env2 - 2�ͱ��- ��Ȩ����������")
cha2.sc <- scores(cha.prda, choices=1:2, display="sp")
arrows(0, 0, cha2.sc[, 1], cha2.sc[, 2], length=0, lty=1, col="red")