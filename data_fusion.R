#----------------------------------------------------------------------------IMPORT PACKAGES-----------------------------------------------------------------------------------------


library(data.table)
library(dplyr)
library(stats)
library(Rcpp)
library(plotly)


#-------------------------------------------------------------------------------GET DATA-------------------------------------------------------------------------------------------

#import MS and RAMAN data by reading .csv files and store in data tables
#RAMAN data


raman_path1 <- "project/RamanSpectra/220511_Proteine_BatchSEC5_10015.csv" # path to RAMAN .csv files
raman_path2 <- "project/RamanSpectra/220511_Proteine_BatchSEC5_11735.csv"
raman_path3 <- "project/RamanSpectra/220511_Proteine_BatchSEC5_12460.csv"


#ms features and ms spectra data


ms_feat_path <- "project/MassSpectra/ms_features.csv" # path to MS .csv files
ms_spec_path <- "project/MassSpectra/ms_spectra.csv" 
spec_lc_path <- "project/MassSpectra/spectra_lc_ms.csv"  



#function to read samples stored as .csv files. fread() reads csv and stores it in a dataTable.

get_csv <- function(path) {
  
  file <- fread(path)
  return(file)
  
}

#--SELECT import file from list above----------------------------------------------------------------------------------------------------------------------------------------
spec_lc <- get_csv(spec_lc_path) # read spec_lc into data.table
raman_data_1 <- get_csv(raman_path1) # read raman 1

spec_lc_list <- split(spec_lc, spec_lc$analysis) # split on basis of analyses

#----------------------------------------------------------------------------END OF GET DATA-----------------------------------------------------------------------------------




#------------------------------------------------------------------------------PLOT DATA 1------------------------------------------------------------------------------------

# visualize (raw)data before processing


# prepare MS data
df <- spec_lc_list[[1]]
df <- df[, `:=`(intensity = mean(intensity)), by = c("rt")][]
df <- unique(df[, c("rt", "intensity"), with = FALSE])
plot(df$rt, df$intensity)

# prepare RAMAN data
df1 <- unique(raman_data_1[, c("shift", "intensity"), with = FALSE])
plot(df1$shift, df1$intensity)


# INTERACTIVE plots
interPlot_data <- function(data){
  
  if ("shift" %in% colnames(data) | "Features_ID" %in% colnames(data)) {
    
    plot <- plot_ly(data,
                    x = if ("shift" %in% colnames(data)) {
                          ~shift } 
                        else {
                          ~Features_ID},
                    y = ~intensity, 
                    type = "scatter", 
                    mode = "lines+markers", 
                    text = ~paste("shift",
                                  if ("shift" %in% colnames(data)) {
                                    ~shift } 
                                  else {
                                    ~Features_ID},
                                  '<br>intensity:',
                                  intensity),
                    color = ~intensity, size = ~intensity)
    
  } else {
    
    plot <- plot_ly(data, x = ~rt, y = ~intensity,
                    type = "scatter",
                    mode = "lines+markers", 
                    text = ~paste("rt", rt, '<br>intensity:', intensity),
                    color = ~intensity, size = ~intensity)
    
  }
  plot
  
}

interPlot_data(df)
interPlot_data(df1)


#--------------------------------------------------------------------------END OF PLOT DATA 1------------------------------------------------------------------------------------------
# kernel density plots give probability density of different intensities

# d <- density(raman_data_1$intensity)
# plot(d, lwd = 2, main = "raman data")
# 
# d1 <- density(spec_lc$intensity)
# plot(d1, lwd = 2, main = "spec_lc")

#----------------------------------------------------------------------------PROCESS DATA------------------------------------------------------------------------------------------------


# preprocessing step 1 : breaks to bin respective data and unify x-axis


# ******************************NOTE : ONLY MS data requires binning*****************************


#bin function to create unified x-axis and facilitate meaningful concatenation of samples.
#set bin size according to data and user input

set_bin <- function(data, rt_bin_size, mz_bin_size) {

  max_rt <- max(data$rt)
  min_rt <- min(data$rt)

  max_mz <- max(data$mz)
  min_mz <- min(data$mz)

  rts_all <- seq(round(min_rt, digits = 0), round(max_rt , digits = 0), rt_bin_size)#last bin contains all entries larger than the minimum value of the corresponding bin
  mzs_all <- seq(round(min_mz, digits = 0), round(max_mz , digits = 0), mz_bin_size)
  
  bins_number <- length(rts_all) + length(mzs_all)
  bins <- matrix(rep(NA, bins_number), nrow = max(length(rts_all), length(mzs_all)), ncol = 2)
  bins[1:length(rts_all), 1] <- rts_all
  bins[1:length(mzs_all), 2] <- mzs_all
  
  return(bins)
}


#--GET INPUT to tweak bin size---------------------------------------------------------------------------------------------------------------------------------------------------


rt_bin_size <- 5
mz_bin_size <- 5


#--END OF GET INPUT--------------------------------------------------------------------------------------------------------------------------------------------------------------


#create bins and necessary values to be fed into binning function.

make_bin <- function(data){
  
  bins <- set_bin(data, rt_bin_size, mz_bin_size) 
  
  
  rts_all <- bins[ , 1]
  rts_all <- na.omit(rts_all)
  
  mzs_all <- bins[ , 2]
  mzs_all <- na.omit(mzs_all)
  
  
  bins_number <- length(rts_all) * length(mzs_all)
  bins_id <- rep(NA_character_, bins_number)
  mat <- matrix(rep(1, bins_number * 2), nrow = 2, ncol = bins_number)
  
  counter <- 0
  for (i in seq_len(length(rts_all))) {
    for (j in seq_len(length(mzs_all))) {
      bins_id[counter + j] <- paste0(rts_all[i], "-", mzs_all[j])
      mat[1, counter + j] <- rts_all[i]
      mat[2, counter + j] <- mzs_all[j]
    }
    counter <- counter + j
  }
  
  mat2_rowNames <- character()
  
  if (uniqueN(data$analysis) > 1) {
    
    spec_lc_list <- split(data, data$analysis) # split on basis of analyses.
    mat2_rowNames <- names(spec_lc_list)
    
    mat2 <- matrix(rep(1, bins_number * length(spec_lc_list)), nrow = length(spec_lc_list), ncol = bins_number, dimnames = list(
      "rows" = mat2_rowNames,# unique analyses
      "cols" = bins_id
    ))
    
  } else { # data is an individual MS sample
    
    mat2 <- matrix(rep(1, bins_number), nrow = 1, ncol = bins_number, dimnames = list(
      "rows" = "sample 1",
      "cols" = bins_id
    ))
    
  }
  
  half_rt_bin <- rt_bin_size/2
  half_mz_bin <- mz_bin_size/2
  
  Rcpp::sourceCpp("C:/Users/PC0118/Desktop/DataFusion/project/fillMat.cpp")
  
  if (length(mat2_rowNames) > 0) {
    
    for (i in mat2_rowNames) {
      
      comparator <- (data$analysis == i) 
      ints <- data$intensity[comparator]
      rts <- data$rt[comparator]
      mzs <- data$mz[comparator]
      ints_out <- fillMat(ints, rts, mzs, rts_all, mzs_all, half_rt_bin, half_mz_bin)

      mat2[i, ] <- ints_out
      
    }
    
  } else {
    
    ints <- data$intensity
    rts <- data$rt
    mzs <- data$mz
    ints_out <- fillMat(ints, rts, mzs, rts_all, mzs_all, half_rt_bin, half_mz_bin)

    mat2[1, ] <- ints_out
    
  }
  
  mat2[is.na(mat2)] <- 0
  return(mat2)
  
}

#--INPUT sample to be binned--------------------------------------------------------------------------------------------------------------------------------------------------------

mat2 <- make_bin(spec_lc)
# mat3 <- make_bin(ms_spec)


bins_id <- colnames(mat2)

ms_sample_1 <- mat2[1, ]
ms_sample_2 <- mat2[2, ]

ms_index <- c(1:length(bins_id))

# ms_sample_1 <- data.table(ms_index, ms_sample_1)
# ms_sample_2 <- data.table(ms_index, ms_sample_2)
# colnames(ms_sample_1) <- c("rt", "intensity")
# colnames(ms_sample_2) <- c("rt", "intensity")

plot(ms_sample_1)
plot(ms_sample_2)

# plot_data(raman_data_1)
interPlot_data(raman_data_1)

#interactive 3-D plot of Mat2

fig <- plot_ly(z = mat2) %>% add_surface()
fig


#--------------------------------------------------------------------------END OF PROCESS DATA----------------------------------------------------------------------------------
# d1 <- density(ms_sample_1)
# plot(d1, lwd = 2, main = "spec_lc")

#heatmap(mat2)


#-----------------------------------------------------------------------------SCALE DATA----------------------------------------------------------------------------------------


# preprocessing step 2 : mean-centering, scaling and block-weighting to respective RAMAN and MS samples.
# this ensures uniform distribution of data points and weighting such that a larger dataset does not
# cause a bias in represented data.

#standardize and center ms/raman/fused dataset: inbuilt scale() applies mean centering and 
#standardization to unit variance

scaled_ms_1 <- scale(ms_sample_1, center = TRUE, scale = TRUE) # scale() performs mean-centering and scaling to unit variation by default
scaled_ms_2 <- scale(ms_sample_2)                            
scaled_raman_1 <- data.table(raman_data_1$shift, scale(raman_data_1$intensity)) 

colnames(scaled_raman_1) <- c("shift", "intensity")

#block-weighting : since one sample may have more datapoints than the other,
#it may completely crowd the output space and lead to a bias while visualizing the readings. Block-weighting 
#multiplies the sample by the inverse square of its size(no. of datapoints).

put_weight <- function(data1, data2) {
  
  part1 <- nrow(data1)
  part2 <- nrow(data2)
  full <- part1 + part2
  prob1 <- (part1/full) * 100
  prob2 <- (part2/full) * 100
  contrib <- 50 #comparing 2 samples, contribution is 50% for each
  weight <- min(contrib/prob1, contrib/prob2)
  newRows <- max(part1, part2) * weight 
  return(newRows)
  
}

newRows <- round(put_weight(scaled_ms_1, scaled_raman_1))

#we give weights to larger(MS) dataset and pick the first x columns given by the calculated weight, perform fusion and do the same for the next x datapoints in the sample

scaled_ms_1 <- scaled_ms_1[1:newRows, ] #similar number of rows as RAMAN, better representation of each sample
scaled_ms_2 <- scaled_ms_2[1:newRows, ]
#unified x-axis, MS data will be concatenated to RAMAN data with max(shift) + 1 = min(features)
#here, features are the bins that represent individual rt-mz pairs for datapoints.

cutoff <- max(scaled_raman_1[ , 1])
x_axis_sample_1 <- ms_index[1:length(scaled_ms_1)] + cutoff
x_axis_sample_2 <- x_axis_sample_1

Sample_ms_1 <- data.table(x_axis_sample_1, scaled_ms_1)
colnames(Sample_ms_1) <- c("Features_ID", "intensity")

Sample_ms_2 <- data.table(x_axis_sample_2, scaled_ms_2)
colnames(Sample_ms_2) <- c("Features_ID", "intensity")


plot(scaled_raman_1)
plot(Sample_ms_1)
plot(Sample_ms_2)
interPlot_data(Sample_ms_1)
#--------------------------------------------------------------------------END OF SCALE DATA------------------------------------------------------------------------------------

d1 <- density(Sample_ms_1$intensity)
plot(d1, lwd = 2, main = "spec_lc")


#----------------------------------------------------------------------------SELECT DATA----------------------------------------------------------------------------------------

#---------check 5 bins, form grid to overlay on raw data
# choose raw data as grid backdrop

spec_lc <- spec_lc
select_data <- function(bins_id, newRows){
  
  selected_bins <- bins_id[1:newRows]
  bin_rt <- substr(selected_bins, 1, 1)
  bin_mz <- substr(selected_bins, 3, 4)
  data <- spec_lc[spec_lc$analysis == names(spec_lc_list)[[1]]]
  data <- data[(data$rt <= max(bin_rt) & data$rt >= min(bin_rt)) &
               (data$mz <= max(bin_mz) & data$mz >= min(bin_mz))]
  return(data)
  
}

selected_data <- select_data(bins_id, newRows)
View(selected_data)

#--------------------------------------------------------------------------END OF SELECT DATA----------------------------------------------------------------------------------

selected_ms_1 <- data.table(c(1:nrow(selected_data)), scale(selected_data$intensity))
colnames(selected_ms_1) <- c("Features_ID", "intensity")
cutoff <- max(scaled_raman_1[ , 1])
x_axis_sample_1 <- ms_index[1:nrow(selected_ms_1)] + cutoff
selected_ms_1$Features_ID <- x_axis_sample_1

#-----------------------------------------------------------------------------FUSE DATA-----------------------------------------------------------------------------------------


# preprocessing step 3 : concatenating RAMAN and MS datasets after respective block-scaling. 
# MS data would be attached to the end of RAMAN data, where bin 901-200 begins after max(shift).
# bins are given an index starting at 1. Depends on bin size, tweaked through function that takes
# user-defined bin size as argument.

#concatenate each sample to same RAMAN sample to visualize and perform PCA.

fuse_data <- function(data1, data2){
  
  colnames(data1) <- c("Features_ID", "intensity")
  colnames(data2) <- c("Features_ID", "intensity")
  fused_data <- rbind(data1, data2)
  return(fused_data)

}


fused_data_1 <- fuse_data(scaled_raman_1, Sample_ms_1)
fused_data_2 <- fuse_data(scaled_raman_1, selected_ms_1)


fwrite(fused_data_1, "fusion.csv")


#--------------------------------------------------------------------------END OF FUSE DATA---------------------------------------------------------------------------------------





#----------------------------------------------------------------------------PLOT DATA 2-----------------------------------------------------------------------------------------


#plots to visualize merged datasets after necessary preprocessing(binning, scaling, etc.)
df <- fused_data_1
df1 <- fused_data_2

plot(scaled_raman_1)
plot(Sample_ms_1)
plot(Sample_ms_2)
plot(x = df$Features_ID, y = df$intensity, type = "l", main = "Features_ID vs. Intensity of fused sample")

interPlot_data(scaled_raman_1)
interPlot_data(df)
interPlot_data(df1)


#----------------------------------------------------------------------------END OF PLOT DATA 2---------------------------------------------------------------------------------- 





#---------------------------------------------------------------------------------PCA-------------------------------------------------------------------------------------------


# preprocessing step 4 : PCA is conducted on processed datasets.


#function to perform PCA on datasets and extract relevant info

pca <- function(data){

  data_pca <- prcomp(data, center = TRUE, scale. = TRUE)
  return(data_pca)

}

ms_pca <- pca(ms_data)
ms_pca2 <- pca(ms_data2)

raman_pca <- pca(raman_data_1)

fusion_pca <- pca(fused_data1)

summary(ms_pca2)

pca_components <- ms_pca$rotation
pca_scores <- ms_pca$x



#------------------------------------------------------------------------------END OF PCA---------------------------------------------------------------------------------------