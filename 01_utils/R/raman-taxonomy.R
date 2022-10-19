library(reticulate)
library(umap)
library(gplots)
library(plotrix)
library(RColorBrewer)
library(lsa)

setwd('/Users/ryanyork/Documents/Research/github/raman-taxonomy/')
source('01_utils/R/arcadia.pal.R')
np <- import("numpy")

#############################
#####Load and clean data#####
#############################
#Load Ho et al. 2019 data
x = np$load('ho_et_al_2019/data/X_finetune.npy')
y = np$load('ho_et_al_2019/data/y_finetune.npy')

#Load config (contains strain name info, etc.)
config = source_python('ho_et_al_2019/config.py')

#Get names by ordering strains (have to +1 order since its 0 indexed as in python; R wants to start with 1)
n = unlist(STRAINS)

#Load taxonomic groups
taxa = read.csv('ho_et_al_2019/data/taxonomic_groups.csv')

#Calculate strain means
m = lapply(split(as.data.frame(x), y), function(z) colMeans(z))
m = do.call(rbind, m)
rownames(m) = n

#Calculate strain standard error
se = lapply(split(as.data.frame(x), y), function(z){
  apply(z, 2, function(w) plotrix::std.error(w))
})
se = do.call(rbind, se)
rownames(se) = n

##################
#####Compare######
##################
#Correlate means
corr = cor(t(m))

#Heatmap
heatmap.2(corr,
          trace = 'none',
          col = rev(brewer.pal(11, 'RdBu')),
          symbreaks = TRUE)

#PCA
pca = prcomp(x)

#Variance explained
plot(cumsum((pca$sdev)^2/sum(((pca$sdev)^2)))[1:100]*100,
     ylab = '% Variance explained',
     xlab = 'n PCs', 
     pch = 20,
     cex.axis = 1.5, 
     cex.lab = 1.5, 
     bty = 'n', 
     ylim = c(0,100))

#Plot
par(mfrow = c(1,3))
plot(pca$x[,1:2])
plot(pca$x[,2:3])
plot(pca$x[,3:4])

#UMAP on all data
u_all = umap(x, verbose = TRUE)

#UMAP on means
u_means = umap(m, verbose = TRUE)

#####################################################
#####GLMs on PCs comparing taxonomy as predictor#####
#####################################################
#Set up predictor matrix
pred = as.data.frame(cbind(apply(taxa, 2, function(x) rep(x, each = 100))))
pred$out = pca$x[,1]

#GLM on PCA using indvidual taxonomic units
mods = list()
for(i in 1:(ncol(pred)-1)){
  mods[[colnames(pred)[i]]] = lm(pred$out~pred[,i])
}

#Add comparison with all units
mods[['All']] = lm(out~., data = pred)

#Calculate BIC
bics = lapply(mods, function(x) BIC(x))

#Plot
bics = sort(unlist(bics))

o = order(bics)
cols = c(arcadia.pal(n = 6, name = 'Accent'), arcadia.pal(n = 6, name = 'Lighter_accents'))[1:length(bics)]
plot(bics[o],
     pch = 20,
     cex = 2,
     col = cols,
     ylab = 'BIC',
     xlab = '',
     cex.lab = 1.5,
     cex.axis = 1.5,
     xaxt = 'n',
     bty = 'n')
axis(1, 1:length(bics), names(bics)[o], cex.axis = 1.5, las = 2)

######################################################
#####Compare cosine similarity by taxonomic group#####
######################################################
#All by all cosine similarity
d = lsa::cosine(t(x))

#Distribution
hist(unlist(as.data.frame(d)))

#Compare taxonomic units
res = list()
for(i in 1:(ncol(pred)-1)){
  
  tmp = split(as.data.frame(d), pred[,i])
  for(j in 1:length(tmp)){
    tmp[[j]] = tmp[[j]][,pred[,i]%in%names(tmp)[j]]
    tmp[[j]] = tmp[[j]][!tmp[[j]] == 1]
    tmp[[j]] = mean(tmp[[j]])
  }
  res[[colnames(pred)[i]]] = unlist(tmp)
}

#Plot
m1 = unlist(lapply(res, function(y) mean(y)))
se1 = unlist(lapply(res, function(y) plotrix::std.error(y)))

o = order(m1, decreasing = TRUE)
cols = c(arcadia.pal(n = 6, name = 'Accent'), arcadia.pal(n = 6, name = 'Lighter_accents'))[1:length(m1)]
plot(m1[o],
     pch = 20,
     cex = 2,
     col = cols,
     ylim = c(0.92, 0.97),
     ylab = 'Cosine similarity',
     xlab = '',
     cex.lab = 1.5,
     cex.axis = 1.5,
     xaxt = 'n',
     bty = 'n')
axis(1, 1:length(m1), names(m1)[o], cex.axis = 1.5, las = 2)
for(i in 1:length(se1)){
  segments(i, m1[o][i]-se1[o][i], i, m1[o][i]+se1[o][i], col = cols[i], lwd = 1.5)
}


