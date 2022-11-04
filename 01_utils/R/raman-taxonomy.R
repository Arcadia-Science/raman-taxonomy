##TO DO:
#All by all comparisons between species using binned phylogenetic signal -> compare distance between species as a function of tree time (x vs y)
#Calculate phylogenetic signal using windows of each spectra e.g. 50 wavenumbers per window (can do both mean per window and/or pca on each window)
#Try removing S. agalactiae
#Plot phylosig via spectra alongside mean species traces

#Set working directory
setwd('/Users/ryanyork/Documents/Research/github/raman-taxonomy/')

#Source utilities
source('01_utils/R/raman-taxonomy-utils.R')

###################################################
#####Load reference bio molecule spectral data#####
###################################################
##Downloaded from https://www.raman.ugent.be/reference-database-raman-spectra-biological-molecules-0
#List .spc files to load
files = list.files('00_data/degelder_et_al_2007/spcs/')

#Create empty list to load into
spcs = list()

#Loop through and load
for(i in 1:length(files)){
  spcs[[gsub('.spc', '', files[i])]] = find_spectral_peaks(paste('00_data/degelder_et_al_2007/spcs/', files[i], sep = ''),
                                                           plot = TRUE)
}

#############################
#####Load and clean data#####
#############################
#Load Ho et al. 2019 data
x = np$load('00_data/ho_et_al_2019/X_finetune.npy')
y = np$load('00_data/ho_et_al_2019/y_finetune.npy')

#Reverse
x = x[,ncol(x):1]

#Load config (contains strain name info, etc.)
config = source_python('00_data/ho_et_al_2019/config.py')

#Get names by ordering strains (have to +1 order since its 0 indexed as in python; R wants to start with 1)
n = unlist(STRAINS)

#Load taxonomic groups
taxa = read.csv('00_data/ho_et_al_2019/taxonomic_groups.csv')

#Load wavenumbers
wav = rev(np$load('00_data/ho_et_al_2019/wavenumbers.npy'))

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
pred = as.data.frame(cbind(apply(taxa, 2, function(x) rep(x, each = 100))))
cols = WGCNA::standardColors(length(unique(pred$Strain)))
#cols = c(arcadia.pal(n = 6, name = 'Accent'), arcadia.pal(n = 6, name = 'Lighter_accents'))[1:length(unique(taxa$Genus))]
#names(cols) = unique(taxa$Genus)
names(cols) = unique(pred$Strain)
#cols = cols[match(pred$Genus, names(cols))]
cols = cols[match(pred$Strain, names(cols))]

par(mfrow = c(1,3))
plot(pca$x[,1:2], pch = 20, col = cols, cex = 0.5)
plot(pca$x[,2:3], pch = 20, col = cols, cex = 0.5)
plot(pca$x[,3:4], pch = 20, col = cols, cex = 0.5)

#UMAP on all data
u_all = umap(x, verbose = TRUE)

#Plot
plot(u_all$layout, 
     pch = 20, 
     col = cols, 
     cex = 0.5,
     xlab = 'Dim 1',
     ylab = 'Dim 2',
     cex.axis = 1.5, 
     cex.lab = 1.5)

#UMAP on means
u_means = umap(m, verbose = TRUE)

#ICA
ica_res = ica::ica(m, nc = 10)

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
names(cols) = c('Strain', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Division', 'Domain')
cols = cols[match(names(nic[o]), names(cols))]
cols = c(cols[1], 'black', cols[2:length(cols)])
plot(bics[o],
     pch = 20,
     cex = 3,
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
for(i in 1:(ncol(pred))){
  
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
names(cols) = c('Strain', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Division', 'Domain')
cols = cols[match(names(m1[o]), names(cols))]
plot(m1[o],
     pch = 20,
     cex = 3,
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
  segments(i, m1[o][i]-se1[o][i], i, m1[o][i]+se1[o][i], col = cols[i], lwd = 2)
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
phylo = drop.tip(phylo, 'Streptococcus_agalactiae')

#Generate matrix for calculate means by species
m2 = m
rownames(m2) = paste(taxa$Genus, taxa$Species, sep = '_')

#Split by species
m2 = split(as.data.frame(m2), rownames(m2))

#Calculate mean by species
m2 = lapply(m2, function(x) colMeans(x))

#Recombine
m2 = do.call(rbind, m2)
m2 = m2[!rownames(m2) == 'Streptococcus_agalactiae',]

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

#Distance between trees
TreeDistance(phylo, spectra)

#Phylogenetic signal of PCs
phylosig_pcs = apply(pca$x[match(phylo$tip.label, rownames(pca$x)),], 2, function(x) phylosig(phylo, x, method = 'lambda')[1]$lambda)

#Phylogenetic signal of spectra
phylosig_spectra = apply(m2[match(phylo$tip.label, rownames(m2)),], 2, function(x) phylosig(phylo, x, method = 'lambda')[1]$lambda)

#By windows
len = 100
win = split_with_overlap(t(m2), len, len-1)[1:(ncol(m2)-len)]
win_phy = list()
for(i in 1:length(win)){
  p = prcomp(t(win[[i]]))
  p = p$x[match(phylo$tip.label, rownames(p$x)),2]
  
  mn = colMeans(win[[i]])
  mn = phylosig(phylo, mn, method = 'lambda')[1]$lambda
  p = phylosig(phylo, p, method = 'lambda')[1]$lambda
  l = list(mn, p)
  names(l) = c('mean', 'pca')
  win_phy[[as.character(i)]] = l
}

pc = unlist(lapply(win_phy, function(x) x$pca))
mn = unlist(lapply(win_phy, function(x) x$mean))

#Plot
#Set up plot
par(mfrow = c(2,1))

n = c(78, which(round(wav)%in% seq(500, 1750, 250)))
plot(phylosig_spectra,
     type = 'l',
     xaxt = 'n',
     ylab = 'Phylogenetic signal',
     cex.axis = 1.5,
     cex.lab = 1.5,
     bty = 'n',
     xlab = 'Wavenumber',
     ylim = c(0,1))
axis(1, n, seq(500, 1750, 250), cex.axis = 1.5)

#Order
ord = rev(phylo$tip.label)

#Get colors (to color by family)
g = unlist(lapply(strsplit(ord, '_'), function(x) x[1]))
cols = c(arcadia.pal(n = 6, name = 'Accent'), arcadia.pal(n = 6, name = 'Lighter_accents'))[1:length(unique(g))]
names(cols) = unique(g)
cols = cols[match(g, names(cols))]

#Extract spectra
tmp = x[grep(ord[1], paste(pred$Genus, pred$Species, sep = '_')),]

#Initiate plot
plot(tmp[1,][1:length(phylosig_spectra)], 
     col = alpha(cols[1], 0.05), 
     type = 'l', 
     bty = 'n', 
     xaxt = 'n', 
     yaxt = 'n', 
     cex.axis = 1.5,
     xlab = 'Wavenumber', 
     ylab = '')

#Loop through and plot
for(i in 2:length(ord)){
  
  #Extract spectra
  tmp = x[grep(ord[i], paste(pred$Genus, pred$Species, sep = '_')),]
  
  #Add mean
  lines(m2[grep(ord[i], rownames(m2)),][1:length(phylosig_spectra)], lwd = 2, col = darken_color(cols[i], factor = 2))
}
axis(1, n, seq(500, 1750, 250), cex.axis = 1.5)

#################################################################
#####Comparing phylogenetic distance and phenotypic distance#####
#################################################################
#Calculate cophenic distances
cophen = cophenetic.phylo(phylo)

#Convert to three column matrix
cophen = as.data.frame(as.table(cophen))

#Calculate cosine similarity
d = cosine(t(m2))
d = as.data.frame(as.table(d))

#Calculate by time
cophen_time = list()
for(i in 0:round(max(cophen$Freq))){
  
  #Filter on cophenetic distances
  z = d[which(cophen$Freq>=i),]
  
  #Remove self comparisons
  z = z[!z[,1] == z[,2],]
  
  #Add to list
  cophen_time[[as.character(i)]] = z
}

#Calculate mean
p_mean = unlist(lapply(cophen_time, function(x) mean(x$Freq)))
p_se = unlist(lapply(cophen_time, function(x) std.error(x$Freq)))

psd<-p_mean+p_se
nsd<-p_mean-p_se

#Plot
plot(1:length(p_mean),
     p_mean,
     type = 'l',
     ylab = 'Spectral similarity',
     cex.axis = 1.5,
     cex.lab = 1.5,
     xaxt = 'n',
     bty = 'n',
     lwd = 1.5,
     xlab = 'Million years')
polygon(x=c(p_mean, rev(p_mean)), y=c(psd, rev(nsd)), col="lightblue", density = 40, angle=90)
axis(1, )

#####################################
#####Plot phylogeny with spectra#####
#####################################
#Plot phylogeny
plot(phylo)

##Plot spectra
#Order
ord = rev(phylo$tip.label)

#Get colors (to color by family)
g = unlist(lapply(strsplit(ord, '_'), function(x) x[1]))
cols = c(arcadia.pal(n = 6, name = 'Accent'), arcadia.pal(n = 6, name = 'Lighter_accents'))[1:length(unique(g))]
names(cols) = unique(g)
cols = cols[match(g, names(cols))]

#Set up plot
par(mfrow = c(length(ord),1), mar = c(0.25, 0.25, 0.25, 0.25))

#Loop through and plot
for(i in 1:length(ord)){
  
  #Extract spectra
  tmp = x[grep(ord[i], paste(pred$Genus, pred$Species, sep = '_')),]
  
  #Initiate plot
  plot(tmp[1,], col = alpha(cols[i], 0.05), type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
  
  #Add lines
  for(j in 2:nrow(tmp)){
    lines(tmp[j,], col = alpha(cols[i], 0.05))
  }
  
  #Add mean
  lines(m2[grep(ord[i], rownames(m2)),], lwd = 2, col = darken_color(cols[i], factor = 2))
  
  #Add label
  #title(main = ord[i], font.main = 1)
  
}

########################################
#####Compare spectra to genome size#####
########################################
#Load genome sizes
size = read.csv('00_data/ho_et_al_2019/ho_2019_genome_statistics.csv')

#Match genome size matrix to order of spectral data
size = size[match(rownames(pca$x), gsub(' ', '_', size$Genome)),]

#Calculate pca
pca = prcomp(m2)

#Linear model
summary(lm(pca$x[,2]~size$media_gc_content+size$median_genome_size+size$median_protein_count))

#By wavenumber
par(mfrow = c(1,3))
n = c(78, which(round(wav)%in% seq(500, 1750, 250)))
plot(apply(m2, 2, function(x) cor(x, size$median_genome_size)),
     type = 'l',
     xaxt = 'n',
     ylab = 'Correlation coeff.',
     cex.axis = 1.5,
     cex.lab = 1.5,
     xlab = 'Wavenumber',
     ylim = c(-1,1))
title(main = 'Genome size', font.main = 1, cex.main = 1.5)
axis(1, n, seq(500, 1750, 250), cex.axis = 1.5)

plot(apply(m2, 2, function(x) cor(x, size$median_protein_count)),
     type = 'l',
     xaxt = 'n',
     ylab = 'Correlation coeff.',
     cex.axis = 1.5,
     cex.lab = 1.5,
     xlab = 'Wavenumber',
     ylim = c(-1,1))
title(main = 'Protein count', font.main = 1, cex.main = 1.5)
axis(1, n, seq(500, 1750, 250), cex.axis = 1.5)

plot(apply(m2, 2, function(x) cor(x, size$media_gc_content)),
     type = 'l',
     xaxt = 'n',
     ylab = 'Correlation coeff.',
     cex.axis = 1.5,
     cex.lab = 1.5,
     xlab = 'Wavenumber',
     ylim = c(-1,1))
title(main = 'GC %', font.main = 1, cex.main = 1.5)
axis(1, n, seq(500, 1750, 250), cex.axis = 1.5)

#Compare to phylogenetic signal
par(mfrow = c(2,1))

plot(apply(m2, 2, function(x) cor(x, size$media_gc_content)),
     type = 'l',
     xaxt = 'n',
     ylab = 'Correlation coeff.',
     cex.axis = 1.5,
     cex.lab = 1.5,
     xlab = 'Wavenumber',
     ylim = c(-1,1))
title(main = 'Genome size', font.main = 1, cex.main = 1.5)
axis(1, n, seq(500, 1750, 250), cex.axis = 1.5)

plot(wav, phylosig_spectra, type = 'l')

cor(apply(m2, 2, function(x) cor(x, size$median_protein_count)), 
    phylosig_spectra)

################################################################
#####Comparing bio molecule spectra to phylogenetic results#####
################################################################
#Associate phylogenetic signal 'peaks' with bio molecules
p = pracma::findpeaks(phylosig_spectra, minpeakheight = 0.2, minpeakdistance = 10)
p = p[order(p[,2]),]

#Plot
plot(phylosig_spectra, type = 'l')
points(p[,2], p[,1], pch = 20, col = 'red')

#Find closest peaks to each 
mol_peaks = list()
for(i in 1:nrow(p)){
  r = list()
  for(j in 1:length(spcs)){
    z = spcs[[j]][abs(p[i,2]-spcs[[j]]$peak)<=10,]
    if(nrow(z)>0){
      r[[names(spcs)[j]]] = z
    }
  }
  if(length(r)>0){
    r = do.call(rbind, r)
    r = r[order(r$intensity, decreasing = TRUE),]
    mol_peaks[[as.character(wav[p[i,2]])]] = r
  }
}

##Get windows corresponding to phylogenetic signal 'peaks'
#Split into peak 'windows'
x = which(abs(diff(phylosig_spectra))>0)
idx <- c(0, cumsum(abs(diff(x)) > 1))
z = split(x, idx)

#Filter on size
z = z[lapply(z, function(x) length(x))>10]

#Convert to wave number
z = lapply(z, function(x) round(wav[x]))



#Split wavenumbers on



