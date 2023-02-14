library(here)
here::here('../') # Should be set to the root of the repository
setwd(here::here())

suppressPackageStartupMessages(source('01_utils/R/raman-taxonomy-utils.R'))

# Load Ho et al. 2019 spectral data
x <- np$load("00_data/ho_et_al_2019/X_finetune.npy")

# Reverse to ensure wavenumbers are ascending
x <- x[, ncol(x):1]

# Load config (contains strain name info, etc.)
config <- source_python("00_data/ho_et_al_2019/config.py")

# Get names by ordering strains (have to +1 order since its 0 indexed as in python; R wants to start with 1)
n <- unlist(STRAINS)

# Load wavenumbers
wav <- rev(np$load("00_data/ho_et_al_2019/wavenumbers.npy"))

# Load taxonomic groups (collected from NCBI)
taxa <- read.csv("00_data/ho_et_al_2019/taxonomic_groups.csv")

# Calculate the mean spectra for all strains by splitting on strain and extracting column means
m <- lapply(split(as.data.frame(x), y), function(z)
  colMeans(z))

# Recombine
m <- do.call(rbind, m)

# Add rownames corresponding to strain ID
rownames(m) <- n

# Calculate spectral standard error for all strains
se <- lapply(split(as.data.frame(x), y), function(z) {
  apply(z, 2, function(w)
    plotrix::std.error(w))
})

# Recombine
se <- do.call(rbind, se)

# Add rownames corresponding to strain ID
rownames(se) <- n

saveRDS(se, '01_utils/ho_et_al_spectra_se.RDS')
saveRDS(m, '01_utils/ho_et_al_spectra_means.RDS')

pred <-
  as.data.frame(cbind(apply(taxa, 2, function(x)
    rep(x, each = 100))))

cols <-
  c(arcadia.pal(n = 6, name = "Accent"),
    arcadia.pal(n = 6, name = "Lighter_accents"))[1:length(unique(taxa$Genus))]
names(cols) <- unique(taxa$Genus)
cols <- cols[match(taxa$Genus, names(cols))]

# Set up plot
par(mfrow = c(15, 2), mar = c(1, 1, 1, 1))

# Loop through and plot
for (i in 1:length(unique(pred$Strain))) {
  # Extract spectra
  tmp <- x[grep(unique(pred$Strain)[i], pred$Strain),]
  
  # Initiate plot
  plot(
    tmp[1,],
    col = alpha(cols[i], 0.05),
    type = "l",
    bty = "n",
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = ""
  )
  
  # Add lines
  for (j in 2:nrow(tmp)) {
    lines(tmp[j,], col = alpha(cols[i], 0.05))
  }
  
  # Add mean
  lines(m[i,], lwd = 2, col = darken_color(cols[i], factor = 2))
  
  # Add label
  title(main = unique(pred$Strain)[i], font.main = 1)
}

pca <- prcomp(x)

plot(
  cumsum((pca$sdev) ^ 2 / sum(((
    pca$sdev
  ) ^ 2)))[1:100] * 100,
  ylab = "% Variance explained",
  xlab = "n PCs",
  pch = 20,
  cex = 3,
  cex.axis = 1.5,
  cex.lab = 1.5,
  bty = "n",
  xlim = c(0, 20),
  ylim = c(0, 100)
)

# Generate empty list to save results
mods <- list()

# Loop through and construct GLMs
for (i in 1:(ncol(pred))) {
  mods[[colnames(pred)[i]]] <- lm(pca$x[, 1] ~ pred[, i])
}

bics <- lapply(mods, function(x)
  BIC(x))

# Order the BIC values
bics <- unlist(bics)
bics <- bics[order(bics)]

# Set up colors
cols <-
  c(arcadia.pal(n = 6, name = "Accent"),
    arcadia.pal(n = 6, name = "Lighter_accents"))[1:length(bics)]
names(cols) <-
  c("Species",
    "Strain",
    "Family",
    "Domain",
    "Division",
    "Order",
    "Class",
    "Genus")
cols <- cols[match(names(bics), names(cols))]

# Plot
plot(
  bics,
  pch = 20,
  cex = 4,
  col = cols,
  ylab = "BIC",
  xlab = "",
  cex.lab = 1.5,
  cex.axis = 1.5,
  xaxt = "n",
  bty = "n"
)

# Add axis
axis(1,
     1:length(bics),
     names(bics),
     cex.axis = 1.5,
     las = 2)

d <- lsa::cosine(t(x))

hist(
  unlist(as.data.frame(d)),
  cex.lab = 1.5,
  cex.axis = 1.5,
  xlab = "Cosine similarity",
  main = ""
)

# Generate empty list to save results in
res <- list()

# Loop through the taxonimc group matrix ('pred') and extract cosine similarites for all associated strains
for (i in 1:(ncol(pred))) {
  # Split cosine similarity matrix on taxonomic unit
  tmp <- split(as.data.frame(d), pred[, i])
  
  # Loop through and calculate mean per strain
  for (j in 1:length(tmp)) {
    tmp[[j]] <- tmp[[j]][, pred[, i] %in% names(tmp)[j]]
    tmp[[j]] <- tmp[[j]][!tmp[[j]] == 1]
    tmp[[j]] <- mean(tmp[[j]])
  }
  
  # Add to results list
  res[[colnames(pred)[i]]] <- unlist(tmp)
}

# Mean per group
m1 <- unlist(lapply(res, function(y)
  mean(y)))

# Standard error
se1 <- unlist(lapply(res, function(y)
  plotrix::std.error(y)))

# Variance
v1 <- unlist(lapply(res, function(y)
  var(y)))

# Order variance measures
o <- order(v1)

# Set up colors
cols <-
  c(arcadia.pal(n = 6, name = "Accent"),
    arcadia.pal(n = 6, name = "Lighter_accents"))[1:length(v1)]
names(cols) <-
  c("Species",
    "Strain",
    "Family",
    "Domain",
    "Division",
    "Order",
    "Class",
    "Genus")
cols <- cols[match(names(v1[o]), names(cols))]

# Plot
plot(
  v1[o],
  pch = 20,
  cex = 3,
  col = cols,
  # ylim = c(0.92, 0.97),
  ylab = "Cosine similarity (variance)",
  xlab = "",
  cex.lab = 1.5,
  cex.axis = 1.5,
  xaxt = "n",
  bty = "n"
)

# Add axis
axis(1,
     1:length(v1),
     names(v1)[o],
     cex.axis = 1.5,
     las = 2)

phylo <-
  read.newick("00_data/ho_et_al_2019/ho_2019_species_list_for_timetree.nwk")

phylo <- drop.tip(phylo, "Streptococcus_agalactiae")

# Split the mean spetrum matrix on species name
m2 <-
  split(as.data.frame(m), paste(taxa$Genus, taxa$Species, sep = "_"))

# Calculate mean spectra by species
m2 <- lapply(m2, function(x)
  colMeans(x))

# Recombine
m2 <- do.call(rbind, m2)

# Remove Streptococcus agalactiae
m2 <- m2[!rownames(m2) == "Streptococcus_agalactiae",]

pca <- prcomp(m2)

# Set length
len <- 25

# Extract spectral windows of desired size using the function 'split_with_overlap'
win <- split_with_overlap(t(m2), len, len - 1)[1:(ncol(m2) - len)]

# Generate empty list for saving results
win_phy <- list()

# Loop and calculate phylogenetic signal
for (i in 1:length(win)) {
  # Reduce the dimensions of spectra in the window from size n to 1 via PCA
  p <- prcomp(t(win[[i]]))
  
  # Match PCA results with phylogeny ordering
  p <- p$x[match(phylo$tip.label, rownames(p$x)), 2]
  
  # Reduce the dimensions of spectra in the window from size n to 1 via calculating the mean
  mn <- colMeans(win[[i]])
  
  # Compute phylogenetic signal using spectral means
  mn <- phylosig(phylo, mn, method = "lambda")[1]$lambda
  
  # Compute phylogenetic signal using PC1
  p <- phylosig(phylo, p, method = "lambda")[1]$lambda
  
  # Combine into a list
  l <- list(mn, p)
  
  # Add names
  names(l) <- c("mean", "pca")
  
  # Add to results list
  win_phy[[as.character(i)]] <- l
}

# Extract results for PC ('pc') and spectral mean ('mn') metrics
pc <- unlist(lapply(win_phy, function(x)
  x$pca))
mn <- unlist(lapply(win_phy, function(x)
  x$mean))

options(repr.plot.width = 12, repr.plot.height = 4)
plot(
  wav[1:length(mn)],
  mn,
  type = "l",
  lwd = 1.5,
  bty = "n",
  ylim = c(0, 0.8),
  cex.axis = 1.5,
  cex.lab = 1.5,
  ylab = "Phylogenetic signal",
  xlab = "Wavenumber (cm -1)"
)

# Set length
len <- 25

# Extract spectral windows of desired size using the function 'split_with_overlap'
win <- split_with_overlap(t(m2), len, len - 1)[1:(ncol(m2) - len)]

# Generate empty list for saving results
win_dist <- list()

# Loop and calculate phylogenetic signal
for (i in 1:length(win)) {
  # Calculate distance between strains
  p <- as.matrix(dist(t(win[[i]])))
  
  # Add to results list
  win_dist[[as.character(i)]] <- unlist(as.data.frame(p))
}

# Combine results in matrix
win_dist <- do.call(cbind, win_dist)

# Generate empty list to save results
cophen_times <- list()

# Loop and calculate cophenic distance in windows
for (i in 1:length(win)) {
  # Empty list to save results
  res <- list()
  
  # Calculate distance in window
  d <- as.matrix(dist(t(win[[i]])))
  
  # Convert to three column matrix
  d <- as.data.frame(as.table(d))
  
  # Looop through and calculate as a function of time
  for (j in 0:round(max(p$Freq))) {
    # Filter on cophenetic distances
    z <- d[which(p$Freq >= j),]
    
    # Remove self comparisons
    z <- z[!z[, 1] == z[, 2],]
    
    # Add to list
    res[[as.character(j)]] <- z
  }
  
  # Calculate mean
  p_mean <- unlist(lapply(res, function(x)
    mean(x$Freq)))
  
  # Add to list
  cophen_times[[as.character(i)]] <- p_mean
}

# Combine into matrix
time_mat <- do.call(cbind, cophen_times)

options(repr.plot.width=12, repr.plot.height=4)

# Heatmap
image(t(time_mat),
      xaxt = "n",
      yaxt = "n",
      bty = "n")

# Add line corresponding to max distance
y <-
  t(unlist(lapply(cophen_times, function(x)
    which.max(x) / length(x))))
x <- (1:length(y) / length(y))
lines(x, y, lwd = 3)

# Set length
len <- 25

# Generate empty list to save results
win_tree <- list()

# Loop and calculate tree distances
for (i in 1:length(win)) {
  # Generate a tree for given window using neighbor joining
  z <- as.phylo(nj(dist(t(win[[i]]))))
  
  # Compare to the phylogenetic tree using the 'TreeDistance' function (Robinson-Foulds metric)
  win_tree[[as.character(i)]] <- TreeDist::TreeDistance(phylo, z)
}

options(repr.plot.width=12, repr.plot.height=4)

plot(
  wav[1:length(win_tree)],
  unlist(win_tree),
  type = "l",
  lwd = 1.5,
  ylab = "Distance to genetic tree",
  cex.axis = 1.5,
  cex.lab = 1.5,
  bty = "n",
  xlab = "Wavenumber"
)

options(repr.plot.width = 12, repr.plot.height = 12)
par(mfrow = c(3, 1))

# Phylogenetic signal
plot(
  wav[1:length(mn)],
  mn,
  type = "l",
  lwd = 1.5,
  bty = "n",
  ylim = c(0, 0.8),
  cex.axis = 1.5,
  cex.lab = 1.5,
  ylab = "Phylogenetic signal",
  xlab = "Wavenumber (cm -1)"
)

# Spectral distance
image(t(time_mat),
      xaxt = "n",
      yaxt = "n",
      bty = "n")
y <-
  t(unlist(lapply(cophen_times, function(x)
    which.max(x) / length(x))))
x <- (1:length(y) / length(y))
lines(x, y, lwd = 3)

# Tree distance
plot(
  wav[1:length(win_tree)],
  unlist(win_tree),
  type = "l",
  lwd = 1.5,
  ylab = "Distance to genetic tree",
  cex.axis = 1.5,
  cex.lab = 1.5,
  bty = "n",
  xlab = "Wavenumber"
)

# Load genome statistics
size <-
  read.csv("00_data/ho_et_al_2019/ho_2019_expanded_genome_statistics.csv")

# Clean and convert to data frame
size <-
  as.data.frame(apply(size, 2, function(x)
    gsub("\\,", "", x)))

# Match genome size matrix to order of spectral data
size <- size[match(rownames(pca$x), gsub(" ", "_", size$Genome)),]

# Remove 'Streptococcus_agalactiae'
size <- size[!size$Genome == "Streptococcus agalactiae",]                    

# Set length
len <- 25

# Generate empty list to save results
mods <- list()

# Loop through and calculate
for (i in 1:(ncol(m2) - len)) {
  # Extract mean strain spectra for given window
  d <- m2[, i:(i + len)]
  
  # PCA
  pca <- prcomp(d)
  
  # Fit GLMs for each genomic feature within the given window
  o <-
    apply(size[, 2:ncol(size)], 2, function(x)
      summary(lm(pca$x[, 1] ~ as.numeric(x)))[9])
  
  # Add to list
  mods[[as.character(i)]] <- o
}

# Plot phylogenetic signal
plot(
  wav[1:length(mn)],
  mn,
  type = "l",
  bty = "n",
  col = "#3B9886",
  ylim = c(-0.1, 0.8),
  cex.axis = 1.5,
  cex.lab = 1.5,
  lwd = 1.5,
  ylab = "Phylogenetic signal",
  xlab = "Wavenumber (cm -1)"
)

# Add GC line
lines(wav[1:length(mn)],
      lwd = 1.5,
      unlist(lapply(mods, function(x)
        x$GC.)),
      col = "#F898AE")

mod <- lm(
  mn ~ unlist(lapply(mods, function(x)
    x$Size..Mb.)) +
    unlist(lapply(mods, function(x)
      x$GC.)) +
    unlist(lapply(mods, function(x)
      x$Protein)) +
    unlist(lapply(mods, function(x)
      x$rRNA)) +
    unlist(lapply(mods, function(x)
      x$tRNA)) +
    unlist(lapply(mods, function(x)
      x$Other.RNA)) +
    unlist(lapply(mods, function(x)
      x$Gene)) +
    unlist(lapply(mods, function(x)
      x$Pseudogene))
)

summary(mod)

# Generate empty vector to save results
cors <- c()

# Loop and correlate
for (i in 1:length(mods[[1]])) {
  cors <- c(cors, cor(mn, unlist(lapply(mods, function(x)
    x[i]))))
}

# Add genome feature names
names(cors) <- names(mods[[1]])

options(repr.plot.width=8, repr.plot.height=8)

# Plot
bp <- barplot(
  rev(cors),
  horiz = TRUE,
  xlim = c(-0.8, 0.8),
  las = 1,
  xlab = "Correlation coefficient",
  ylab = "",
  cex.lab = 1.5,
  cex.axis = 1.5,
  names.arg = "",
  col = c(rep("gray80", 6), "#F898AE", "gray80"),
  border = c(rep("gray80", 6), "#F898AE", "gray80")
)

# Add labels
lab <-
  c(
    "Genome size",
    "GC content",
    "Protein #",
    "rRNA #",
    "tRNA #",
    "Other RNA #",
    "Gene #",
    "Pseudogene #"
  )
text(0, bp, rev(lab), cex = 1.5, pos = 4)
