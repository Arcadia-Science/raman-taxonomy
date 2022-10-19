library(reticulate)
library(umap)
library(gplots)
library(plotrix)
library(RColorBrewer)
library(lsa)
library(ape)
library(phytools)
library(phylosignal)
library(scales)

setwd('/Users/ryanyork/Documents/Research/github/raman-taxonomy/')
source('01_utils/R/arcadia.pal.R')
source('01_utils/R/darken_color.R')
np <- import("numpy")

#############################
#####Load and clean data#####
#############################
#Load Ho et al. 2019 data
x = np$load('00_data/ho_et_al_2019/X_finetune.npy')
y = np$load('00_data/ho_et_al_2019/y_finetune.npy')

#Load config (contains strain name info, etc.)
config = source_python('ho_et_al_2019/config.py')

#Get names by ordering strains (have to +1 order since its 0 indexed as in python; R wants to start with 1)
n = unlist(STRAINS)

#Load taxonomic groups
taxa = read.csv('00_data/ho_et_al_2019/taxonomic_groups.csv')

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

######################
#####Plot spectra#####
######################
##Mean plot with standard error
#Get colors (to color by family)
cols = c(arcadia.pal(n = 6, name = 'Accent'), arcadia.pal(n = 6, name = 'Lighter_accents'))[1:length(unique(taxa$Family))]
names(cols) = unique(taxa$Family)
cols = cols[match(taxa$Family, names(cols))]

#Set up plot
par(mfrow = c(15,2), mar = c(1,1,1,1))

#Loop through and plot mean spectra and standard error
for(i in 1:nrow(m)){
  
  #Initiate plot
  plot(m[i,], bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', type = 'n')
  
  #Plot standard error as polygons
  polygon(c(seq(1, length(m[i,]), 1), 
            rev(seq(1, length(m[i,]), 1))),
          c(m[i,]-se[i,],
            rev(m[i,]+se[i,])),
          col = cols[i], 
          border = FALSE)
  
  lines(m[i,], col = darken_color(cols[i]), lwd = 0.75)
}

##Mean plot with all spectra 
#Set up predictor matrix
pred = as.data.frame(cbind(apply(taxa, 2, function(x) rep(x, each = 100))))

#Get colors (to color by family)
cols = c(arcadia.pal(n = 6, name = 'Accent'), arcadia.pal(n = 6, name = 'Lighter_accents'))[1:length(unique(taxa$Genus))]
names(cols) = unique(taxa$Genus)
cols = cols[match(taxa$Genus, names(cols))]

#Set up plot
par(mfrow = c(15,2), mar = c(1,1,1,1))

#Loop through and plot
for(i in 1:length(unique(pred$Strain))){
  
  #Extract spectra
  tmp = x[grep(unique(pred$Strain)[i], pred$Strain),]
  
  #Initiate plot
  plot(tmp[1,], col = alpha(cols[i], 0.05), type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
  
  #Add lines
  for(j in 2:nrow(tmp)){
    lines(tmp[j,], col = alpha(cols[i], 0.05))
  }
  
  #Add mean
  lines(m[i,], lwd = 2, col = darken_color(cols[i], factor = 2))
  
  #Add label
  title(main = unique(pred$Strain)[i], font.main = 1)
  
}

##########################################################
#####Compare overall structure of data via PCA, UMAP######
##########################################################
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

#####################################
#####Correlate taxonomy with pcs#####
#####################################
#PCA
pca = prcomp(x)

#Set up predictor matrix
pred = as.data.frame(cbind(apply(taxa, 2, function(x) rep(x, each = 100))))

#Get all comparisons of pcs and predictors
toTest = expand.grid(1:50, 1:ncol(pred))

#Linear models
res = list()
for(i in 1:nrow(toTest)){
  res[[paste(toTest[i,1], toTest[i,2], sep = '_')]] = summary(lm(pca$x[,toTest[i,1]]~pred[,toTest[i,2]]))[[9]]
}
res = unlist(res)

#Split on pc
res = split(res, toTest[,1])

#Recombine
res = do.call(rbind, res)

#Plot
image(res)

##################################
#####Phylogenetic comparisons#####
##################################
#Load phylogeny
phylo = read.newick('00_data/ho_et_al_2019/ho_2019_species_list_for_timetree.nwk')

#Generate matrix for calculate means by species
m2 = m
rownames(m2) = paste(taxa$Genus, taxa$Species, sep = '_')

#Split by species
m2 = split(as.data.frame(m2), rownames(m2))

#Calculate mean by species
m2 = lapply(m2, function(x) colMeans(x))

#Recombine
m2 = do.call(rbind, m2)

#Generate phylogeny from spectra
spectra = as.phylo(nj(dist(m2)))

#Cophyloplot
obj = cophylo(phylo, spectra)
plot(obj,
     link.type="curved",
     link.lwd=2,
     link.col=make.transparent("grey",0.7))

#By PCs
pca = prcomp(m2)

#Generate phylogeny from spectra
spectra = as.phylo(nj(dist(pca$x[,1:10])))

#Cophyloplot
obj = cophylo(phylo, spectra)
plot(obj,
     link.type="curved",
     link.lwd=2,
     link.col=make.transparent("grey",0.7))

#Phylogenetic signal of PCs
phylosig_pcs = apply(pca$x[match(phylo$tip.label, rownames(pca$x)),], 2, function(x) phylosig(phylo, x))

#Phylogenetic signal of spectra
phylosig_spectra = apply(m2[match(phylo$tip.label, rownames(m2)),], 2, function(x) phylosig(phylo, x))
