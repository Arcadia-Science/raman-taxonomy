library(reticulate)
library(gplots)
library(plotrix)
library(lsa)
library(ape)
library(phytools)
library(phylosignal)
library(scales)
library(ArcadiaColorBrewer)
library(TreeDist)
library(hyperSpec)
library(pracma)

#Import numpy
np <- import("numpy")

#Get windows function
split_with_overlap <- function(vec, seg_length, overlap) {
  starts <- seq(1, nrow(vec), by = seg_length - overlap)
  ends <- starts + seg_length - 1
  ends[ends > nrow(vec)] <- nrow(vec)
  
  lapply(1:length(starts), function(i) {
    vec[starts[i]:ends[i], ]
  })
}

#Function to load spectra (in .spc files) and find peaks
find_spectral_peaks = function(file,
                               plot = FALSE,
                               min_peak_height = 0.1){
  
  #Load
  tmp = read.spc(filename = file)
  
  #Extract spectrum
  p = tmp$spc[1,]
  
  #Normalize
  p = p/max(p)
  
  #Find peaks
  z = pracma::findpeaks(p, minpeakheight = min_peak_height)
  
  #Plot, if desired
  if(plot == TRUE){
    plot(p, 
         type = 'l',
         xlab = 'Wavenumber',
         ylab = 'Intensity (a.u.)',
         cex.axis = 1.5,
         cex.lab = 1.5,
         xaxt = 'n',
         bty = 'n')
    axis(1, seq(100, 700, 200), round(as.numeric(names(p))[seq(100, 700, 200)]), cex.axis = 1.5)
    points(z[,2], z[,1], pch = 20, col = 'red')
  }
  
  #Create peak dataframe
  peaks = data.frame(peak = round(as.numeric(names(tmp$spc[1,]))[z[,2]]),
                     intensity = z[,1])
  
  #Return
  return(peaks)
}