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

#Scale, if desired
spcs_s = spcs
for(i in 1:length(spcs_s)){
  spcs_s[[i]]$intensity = scale(spcs_s[[i]]$intensity)  
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

#Save 
saveRDS(se, '01_utils/ho_et_al_spectra_se.RDS')
saveRDS(m, '01_utils/ho_et_al_spectra_means.RDS')

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
#cols = WGCNA::standardColors(length(unique(pred$Strain)))
cols = c(arcadia.pal(n = 6, name = 'Accent'), arcadia.pal(n = 6, name = 'Lighter_accents'))[1:length(unique(taxa$Genus))]
names(cols) = unique(taxa$Genus)
#names(cols) = unique(pred$Strain)
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
cols = c(arcadia.pal(n = 6, name = 'Accent'), arcadia.pal(n = 6, name = 'Lighter_accents'))[1:length(bics)]
names(cols) = c('Species', 'Strain', 'Family', 'Domain', 'Division', 'Order', 'Class', 'Genus', 'All')
cols = cols[match(names(bics), names(cols))]
cols[2] = 'black'
plot(bics[o],
     pch = 20,
     cex = 4,
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
v1 = unlist(lapply(res, function(y) var(y)))

o = order(v1)
cols = c(arcadia.pal(n = 6, name = 'Accent'), arcadia.pal(n = 6, name = 'Lighter_accents'))[1:length(v1)]
names(cols) = c('Species', 'Strain', 'Family', 'Domain', 'Division', 'Order', 'Class', 'Genus')
cols = cols[match(names(v1[o]), names(cols))]
plot(v1[o],
     pch = 20,
     cex = 3,
     col = cols,
     #ylim = c(0.92, 0.97),
     ylab = 'Cosine similarity (variance)',
     xlab = '',
     cex.lab = 1.5,
     cex.axis = 1.5,
     xaxt = 'n',
     bty = 'n')
axis(1, 1:length(v1), names(v1)[o], cex.axis = 1.5, las = 2)
#for(i in 1:length(se1)){
#  segments(i, m1[o][i]-se1[o][i], i, m1[o][i]+se1[o][i], col = cols[i], lwd = 2)
#}

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
phylosig_spectra_k = apply(m2[match(phylo$tip.label, rownames(m2)),], 2, function(x) phylosig(phylo, x, method = 'K')[1])

#By windows
len = 25
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

##Plot window-based phylogenetic signal
bands = as.data.frame(readxl::read_xlsx('00_data/shipp_2017_raman_bands.xlsx', col_names = FALSE))
colnames(bands) = c('wavenumber', 'end', 'molecule', 'type')
plot(wav[1:length(mn)],
     mn,
     type = 'n',
     bty = 'n',
     ylim = c(0,0.8),
     cex.axis = 1.5,
     cex.lab = 1.5,
     ylab = 'Phylogenetic signal',
     xlab = 'Wavenumber (cm -1)')
#for(i in 1:nrow(bands)){
#  abline(v = bands[i,]$wavenumber, col = 'gray40')
#}
#text(bands$wavenumber, rep(0.5, nrow(bands)), bands$molecule, srt = 90, adj = 0, col = 'gray40')
lines(wav[1:length(mn)],
      mn,
      lwd = 1.5)

##Plot with common bands
bands = as.data.frame(readxl::read_xlsx('00_data/cui_et_al_2022_raman_bands.xlsx', col_names = FALSE))
colnames(bands) = c('wavenumber', 'molecule', 'bond')
plot(wav,
     phylosig_spectra,
     type = 'n',
     bty = 'n',
     ylim = c(0,1.5),
     cex.axis = 1.5,
     cex.lab = 1.5,
     ylab = 'Phylogenetic signal',
     xlab = 'Wavenumber (cm -1)')
for(i in 1:nrow(bands)){
  abline(v = bands[i,]$wavenumber, col = 'gray40')
}
text(bands$wavenumber, rep(1, nrow(bands)), bands$molecule, srt = 90, adj = 0, col = 'gray40')
lines(wav,
      phylosig_spectra,
      lwd = 1.5)

##Plot with variance
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

#Variance
plot(apply(m2, 2, function(x) sd(x)),
     type = 'l',
     xaxt = 'n',
     ylab = 'Standard deviation',
     cex.axis = 1.5,
     cex.lab = 1.5,
     bty = 'n',
     xlab = 'Wavenumber')
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

#Get max divergence times for each taxonomic grouping
times = matrix(ncol = 4)
for(i in 1:nrow(cophen)){
  for(j in 1:ncol(cophen)){
    x = taxa2[grep(rownames(cophen)[i], 
                   paste(taxa2$Genus, taxa2$Species, sep = '_')),][1,]
    y = taxa2[grep(rownames(cophen)[j], 
                   paste(taxa2$Genus, taxa2$Species, sep = '_')),][1,]
    z = names(x)[max(which(x%in%y==FALSE))]
    times = rbind(times, cbind(z, rownames(cophen)[i], colnames(cophen)[j], as.numeric(cophen[i,j])))
  }
}
times = as.data.frame(times)
a = max(as.numeric(times[,4]), na.rm = TRUE)/max(phylo$edge.length)
times[,4] = as.numeric(times[,4])/a

times2 = split(times, times[,1])
lapply(times2, function(x) mean(as.numeric(x[,4]), na.rm = TRUE))

#Convert to three column matrix
cophen = as.data.frame(as.table(cophen))

#Add taxonomic relationship
taxa2 = cbind(taxa$Species, 
              taxa$Genus, 
              taxa[,4:ncol(taxa)])
colnames(taxa2)[1] = 'Species'
colnames(taxa2)[2] = 'Genus'
cophen$relationship = rep(NA, nrow(cophen))
for(i in 1:nrow(cophen)){
  x = taxa2[grep(as.character(cophen[i,1]), 
                 paste(taxa2$Genus, taxa2$Species, sep = '_')),][1,]
  y = taxa2[grep(as.character(cophen[i,2]), 
                 paste(taxa2$Genus, taxa2$Species, sep = '_')),][1,]
  cophen$relationship[i] = names(x)[max(which(x%in%y==FALSE))]
}

#Calculate times for each taxa
z = split(as.data.frame(cophen), cophen$relationship)
z = unlist(lapply(z, function(x) max(x$Freq)))

#Standardize to mya
a = max(cophen$Freq)/max(phylo$edge.length)
cophen$Freq = cophen$Freq/a

#Calculate cosine similarity
d = cosine(t(m2))
d = as.data.frame(as.table(d))

summary(lm(cophen$Freq~d$Freq))

#Calculate cosine similarity for regions with phylogenetic signal
d2 = cosine(t(m2[,which(mn>0.1)]))
d2 = as.data.frame(as.table(d2))

summary(lm(cophen$Freq~d2$Freq))
plot(cophen$Freq, d2$Freq)

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
     xlab = 'Million years ago')
axis(1, at = c(1, seq(500, 3500, 500))+36, labels = rev(seq(0, 3.5, 0.5)), cex.axis = 1.5)

##By window
d_win = list()
step = 50
for(i in seq(1, 1000-step, step)){
  
  #Calculate cosine distance
  d = as.data.frame(as.table(cosine(t(m2[,i:(i+step)]))))
  
  phen = list()
  for(j in 0:round(max(cophen$Freq))){
    
    #Filter on cophenetic distances
    z = d[which(cophen$Freq>=j),]
    
    #Remove self comparisons
    z = z[!z[,1] == z[,2],]
      
    #Add to list
    phen[[as.character(j)]] = z
  }

  #Add to list
  d_win[[as.character(i)]] = unlist(lapply(phen, function(x) mean(x$Freq)))
}

z = do.call(rbind, lapply(d_win, function(x) unlist(x)))
me = colMeans(z)
se = apply(z, 2, function(x) std.error(x))

plot(1:length(me),
     me,
     ylab = 'Spectral similarity',
     cex.axis = 1.5,
     cex.lab = 1.5,
     xaxt = 'n',
     bty = 'n',
     ylim = c(0.99,1),
     xlab = 'Million years',
     type = 'n')
lines(me+se, lty = 2, col = 'grey70', lwd = 1.5)
lines(me-se, lty = 2, col = 'grey70', lwd = 1.5)
lines(me, lwd = 2)
axis(1, at = c(1, seq(500, 3500, 500))+36, labels = rev(seq(0, 3.5, 0.5)), cex.axis = 1.5)
for(i in 1:length(times)){
  abline(v = length(me) - min(times[[i]][,4]))
  abline(v = length(me) - max(times[[i]][,4]))
}

####################################################################
#####Spectral similarity as a function of position (by windows)#####
####################################################################
len = 25
win = split_with_overlap(t(m2), len, len-1)[1:(ncol(m2)-len)]
win_cosine = list()
for(i in 1:length(win)){
  p = as.matrix(dist(t(win[[i]])))

  win_cosine[[as.character(i)]] = unlist(as.data.frame(p))
}

win_cosine = do.call(cbind, win_cosine)

#Sort by divergence time
p = as.matrix(dist(t(win[[1]])))
p = as.data.frame(as.table(p))

z = as.matrix(dist(t(cophen)))
z = as.data.frame(as.table(z))
a = max(z$Freq)/max(phylo$edge.length)
z$Freq = z$Freq/a

p = p[match(paste(z[,1], z[,2]), paste(p[,1], p[,2])),]
p$Freq = z$Freq

win_cosine = win_cosine[order(p$Freq),]
image(t(win_cosine), col = rev(viridis::viridis_pal(option = 'A')(100)))

png('~/Desktop/test.png', width = 1200, height = 400)
image(t(win_cosine[20:nrow(win_cosine),]), xaxt = 'n', yaxt = 'n', bty = 'n')
dev.off()

png('~/Desktop/test2.png', width = 1200, height = 1200)
heatmap(win_cosine, Colv = NA)
dev.off()

#Tree variation as a function of time
cophen_times = list()
for(i in 1:length(win)){
  res = list()
  d = as.matrix(dist(t(win[[i]])))
  d = as.data.frame(as.table(d))
  for(j in 0:round(max(p$Freq))){
    
    #Filter on cophenetic distances
    z = d[which(p$Freq>=j),]
    
    #Remove self comparisons
    z = z[!z[,1] == z[,2],]
    
    #Add to list
    res[[as.character(j)]] = z
  }
  
  #Calculate mean
  p_mean = unlist(lapply(res, function(x) mean(x$Freq)))
  #plot(p_mean, type = 'l')
  
  cophen_times[[as.character(i)]] = p_mean
}

plot(unlist(lapply(cophen_times, function(x) which.max(x))))

#Heatmap
time_mat = do.call(cbind, cophen_times)

png('~/Desktop/spectral_divergence_time_matrix_020123.png', width = 1200, height = 400)
image(t(time_mat), xaxt = 'n', yaxt = 'n', bty = 'n')
dev.off()

#Correlation with baseline
apply(time_mat, 2, function(x) cor(x, p_mean, use = 'complete.obs'))

###################################################################
#####Spectral similarity in phylogenetically conserved regions#####
###################################################################
#Find peaks in window-based phylogenetic signal
a = which(mn>0.01)
idx <- c(0, cumsum(abs(diff(a)) > 2))

indices = split(a, idx)

split(m2,cumsum(a>0.01))

#Tree variation as a function of time
len = 25
win = split_with_overlap(t(m2), len, len-1)[1:(ncol(m2)-len)]
win_times = list()
for(i in 1:length(indices)){
  res = list()
  
  d = as.matrix(dist(m2[,min(indices[[i]]):max(indices[[i]])]))
  d = as.data.frame(as.table(d))
  for(j in 0:round(max(p$Freq))){
    
    #Filter on cophenetic distances
    z = d[which(p$Freq>=j),]
    
    #Remove self comparisons
    z = z[!z[,1] == z[,2],]
    
    #Add to list
    res[[as.character(j)]] = z
  }
  
  #Calculate mean
  p_mean = unlist(lapply(res, function(x) mean(x$Freq)))
  #plot(p_mean, type = 'l')
  
  win_times[[as.character(i)]] = p_mean
}

#Plot
cols = c(arcadia.pal(n = 6, 'Accent'), arcadia.pal(n = 6, 'Lighter_accents'))[1:length(win_times)]
plot(1:length(win_times[[1]]),
     win_times[[1]]/max(win_times[[1]], na.rm = TRUE),
     type = 'l',
     col = cols[1],
     ylab = 'Spectral distance',
     cex.axis = 1.5,
     cex.lab = 1.5,
     xaxt = 'n',
     bty = 'n',
     lwd = 1.5,
     ylim = c(0.3,1),
     xlab = 'Million years ago')
for(i in 2:length(win_times)){
  lines(win_times[[i]]/max(win_times[[i]], na.rm = TRUE),
        col = cols[i],
        lwd = 1.5)
}
lines(1:length(win_times[[1]]), p_mean/max(p_mean, na.rm = TRUE), lwd = 2)
axis(1, at = c(1, seq(500, 3500, 500))+36, labels = rev(seq(0, 3.5, 0.5)), cex.axis = 1.5)

par(mfrow = c(2,1))
plot(apply(win_cosine, 2, function(x) p$Freq[order(p$Freq)][which.max(x)]), type = 'l')
image(t(win_cosine[20:nrow(win_cosine),]), xaxt = 'n', yaxt = 'n', bty = 'n')

#Distance to "baseline"?
d = as.matrix(dist(m2))
d = unlist(as.data.frame(d))

out = apply(win_cosine, 2, function(x) cosine(x, d))

plot(mn[1:length(out)]/max(mn), type = 'l', col = 'red')
lines(out/max(out), type = 'l')

#Time x spectral similarity relationship as a function of position
d_win = list()
step = 25
for(i in seq(1, 1000-step, 1)){
  
  #Calculate cosine distance
  d = as.data.frame(as.table(cosine(t(m2[,i:(i+step)]))))
  
  phen = list()
  for(j in 0:round(max(cophen$Freq))){
    
    #Filter on cophenetic distances
    z = d[which(cophen$Freq>=j),]
    
    #Remove self comparisons
    z = z[!z[,1] == z[,2],]
    
    #Add to list
    phen[[as.character(j)]] = z
  }
  
  #Add to list
  d_win[[as.character(i)]] = unlist(lapply(phen, function(x) mean(x$Freq)))
}

##TO DO: Make a spectral similarity x time plot for phylogenetically conserved windows
##TO DO: This could also be represented as a heatmap for each bin across the spectrum too, comparing the spectra of all species at that point

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

#Remove 'Streptococcus_agalactiae'
size = size[!size$Genome == 'Streptococcus agalactiae',]

#Calculate pca
pca = prcomp(m2)

#Linear model
summary(lm(pca$x[,1]~size$media_gc_content+size$median_genome_size+size$median_protein_count))

#By wavenumber
par(mfrow = c(1,3))
n = c(78, which(round(wav)%in% seq(500, 1750, 250)))
plot(apply(m2, 2, function(x) cor(x, size$median_genome_size)),
     type = 'l',
     xaxt = 'n',
     ylab = 'Correlation coefficient',
     cex.axis = 1.5,
     cex.lab = 1.5,
     xlab = 'Wavenumber',
     ylim = c(-1,1))
title(main = 'Genome size', font.main = 1, cex.main = 1.5)
axis(1, n, seq(500, 1750, 250), cex.axis = 1.5)

plot(apply(m2, 2, function(x) cor(x, size$median_protein_count)),
     type = 'l',
     xaxt = 'n',
     ylab = 'Correlation coefficient',
     cex.axis = 1.5,
     cex.lab = 1.5,
     xlab = 'Wavenumber',
     ylim = c(-1,1))
title(main = 'Protein count', font.main = 1, cex.main = 1.5)
axis(1, n, seq(500, 1750, 250), cex.axis = 1.5)

plot(apply(m2, 2, function(x) cor(x, size$media_gc_content)),
     type = 'l',
     xaxt = 'n',
     ylab = 'Correlation coefficient',
     cex.axis = 1.5,
     cex.lab = 1.5,
     xlab = 'Wavenumber',
     ylim = c(-1,1))
title(main = 'GC %', font.main = 1, cex.main = 1.5)
axis(1, n, seq(500, 1750, 250), cex.axis = 1.5)

#Compare to phylogenetic signal
v = apply(m2, 2, function(x) sd(x))
mg = apply(m2, 2, function(x) minerva::mine(x, size$median_genome_size)$MIC)
mp = apply(m2, 2, function(x) minerva::mine(x, size$median_protein_count)$MIC)
gc = apply(m2, 2, function(x) minerva::mine(x, size$media_gc_content)$MIC)

minerva::mine(phylosig_spectra, mp)

mod = lm(phylosig_spectra~v+mg+mp+gc)

#Linear models by window
len = 25
mods = list()
for(i in 1:(ncol(m2)-len)){
  
  d = m2[,i:(i+len)]
  pca = prcomp(d)
  
  g = summary(lm(pca$x[,1]~size$median_genome_size))[9]
  p = summary(lm(pca$x[,1]~size$median_protein_count))[9]
  gc = summary(lm(pca$x[,1]~size$media_gc_content))[9]
  
  l = list(g, p, gc)
  names(l) = c('genome', 'protein', 'gc')
  
  mods[[as.character(i)]] = l
}

plot(unlist(lapply(mods, function(x) x$gc)), type = 'l', ylim = c(-0.1,0.6))
lines(unlist(lapply(mods, function(x) x$genome)), col = 'red')
lines(unlist(lapply(mods, function(x) x$protein)), col = 'cyan4')



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
    z = spcs_s[[j]][abs(p[i,2]-spcs_s[[j]]$peak)<=2,]
    if(nrow(z)>0){
      r[[names(spcs)[j]]] = z
    }
  }
  if(length(r)>0){
    r = do.call(rbind, r)
    r = r[order(r$intensity, decreasing = TRUE),]
    r$peak_wav = wav[r$peak]
    mol_peaks[[as.character(wav[p[i,2]])]] = r
  }
}

##Get windows corresponding to phylogenetic signal 'peaks'
#Split into peak 'windows'
x = which(abs(diff(phylosig_spectra))>0)
idx <- c(0, cumsum(abs(diff(x)) > 1))
z = split(x, idx)

#Filter on size
z = z[lapply(z, function(x) length(x))>4]

#Convert to wave number
z = lapply(z, function(x) round(wav[x]))

#Find intersections with bio molecule peaks
s = do.call(rbind, spcs_s)
res = list()
for(i in 1:length(z)){
  fishers = list()
  for(j in 1:length(spcs)){
    
    b = s[s$peak%in%z[[i]],]
    b = b[order(b$intensity, decreasing = TRUE),]
    
    a = s[-grep(names(spcs)[j], rownames(s)),]
    
    x = sum(spcs[[j]]$peak%in%z[[i]])
    y = sum(!spcs[[j]]$peak%in%z[[i]])
    
    m = sum(a$peak%in%z[[i]])
    n = sum(!a$peak%in%z[[i]])
    
    f = fisher.test(as.matrix(cbind(c(x,y), c(m,n))))
    fishers[[names(spcs)[j]]] = f
  }
  l = list(fishers, b)
  names(l) = c('fishers', 'peaks')
  res[[names(z)[i]]] = l
}

#Extract results
out = list()
#for(i in 1:length(res)){
#  if(nrow(res[[i]]$peaks)>=5){
#    out[[i]] = res[[i]]$peaks[1:5,]
#  }else{
#    out[[i]] = res[[i]]$peaks
#  }
#}
out = do.call(rbind, lapply(res, function(x) as.data.frame(x$peaks[1:2,])))
out$name =  rownames(out)
out$name = unlist(lapply(strsplit(out$name, '\\.'), function(x) x[2]))
out = out[out$intensity>0.5,]

#Plot
plot(wav,
     phylosig_spectra,
     type = 'n',
     bty = 'n',
     ylim = c(0,1.5),
     cex.axis = 1.5,
     cex.lab = 1.5,
     ylab = 'Phylogenetic signal',
     xlab = 'Wavenumber (cm -1)')
for(i in 1:nrow(out)){
  segments(out[i,]$peak, 0, out[i,]$peak, 1, col = 'gray40')
}
text(out$peak, rep(1, nrow(out)), out$name, srt = 90, adj = 0, col = 'gray40')
lines(wav,
      phylosig_spectra,
      lwd = 1.5)
