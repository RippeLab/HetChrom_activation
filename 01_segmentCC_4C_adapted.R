# Script title: segmentCC_RW
# Overview: thresholding-based segmentation strategies to generate nucleus (Nuc), chromocenter (CC) and nucleoplasm masks. Intensities and pixel area are quantified for the masks to extract relevant image features.
# Use: initially created for analysis of iMEF VPR recruitment to CC + IF (H3K27ac + HP1)

# load libraries and functions
library(EBImage)
library(abind)
source("/Volumes/sd17B002-1/2026_HetChrom/Fig2/BC_HP1H3K9me3stay/_scripts_4Channels/functions/makeCCMask.R")
source("/Volumes/sd17B002-1/2026_HetChrom/Fig2/BC_HP1H3K9me3stay/_scripts_4Channels/functions/makeNucMask.R")

####################################################################################################################
############################################### ADJUSTABLE PARAMETERS ##############################################
####################################################################################################################

# location of the input (maximum-intensity projection .tif files)
folder <- "/Volumes/sd17B002/2026_HetChrom/Fig2/wt/1618_p65/"

# adjustable adaptive thresholding offset for nucleus segmentation
# default is 0.0008
threshold <- 0.0008

# select the channel in which the nuclei will be segmented
# "DNA", "dCas9", "H3K9me3", "HP1"
nucleus_seg_channel <- "DNA"

# remove incomplete nuclei at the image borders ?
# "y" = yes, "n" = no
remove_bordercells <- "y"

# thresholding parameter for CC segmentation in the dCas9 channel using the makeCCMask function
# intensity threshold above which a pixel is considered to belong to a chromocenter (based on dCas9 enrichment) is given by threshold = median(nuclear intensity) + sat_cutoff * (maximum(nuclear intensity) - median(nuclear intensity))
# default is 0.1
sat_cutoff <- 0.15

# CC segmentation will only be carried out on cells that express dCas9 ("positive cells")
# cell is considered "positive" if its mean nuclear dCas9 intensity is "mean_cutoff"-fold above the median (background) of the total image intensity
# this parameter heavily depends on the type of construct used and the signal-to-background ratio
# default is 1.5
mean_cutoff <- 1.5

# size-filter for nucleus segmentation (pixels); needs to be adjusted to cell size in sample
# default is 8000
nucSize_cutoff <- 8000

# maximum number of nuclei expected per image; required for variable initialization
# default is 20
mn <- 20

####################################################################################################################
#################################################### FIXED #########################################################
####################################################################################################################

# create timestamp variable  for this run
date_time <- gsub(" ", "_", gsub(":","-", substr(Sys.time(),1,19)))

# create results folders
dir.create(paste0(folder,date_time,"_segmentCC_results"))
results_folder <- paste0(folder,date_time,"_segmentCC_results")
dir.create(paste0(results_folder,"/",date_time,"_nuc_shape_features"))

# document adjustable segmentation parameters used in this run
parameters_df <- data.frame(threshold, nucleus_seg_channel, remove_bordercells, sat_cutoff, mean_cutoff, nucSize_cutoff, mn)
colnames(parameters_df) <- c("threshold","nucleus_seg_channel","remove_bordercells","sat_cutoff","mean_cutoff","nucSize_cutoff","mn")
write.table(parameters_df, file = paste0(results_folder,"/",date_time,"_segmentation_parameters.csv"), quote = F, sep = "\t", row.names = F)

# generate and write a position list to document the files that were processed in this run
positions <- list.files(path = folder, pattern = "tif")
write.table(cbind(seq(1,length(positions), by = 1), positions), file = paste0(results_folder,"/",date_time,"_positions.csv"), quote = F, sep = "\t", row.names = F)

# initialization of variables and empty vectors for storing quantification results
mean_DNA_intensity <- rep(0, mn*length(positions))
sd_DNA_intensity <- rep(0, mn*length(positions))
mean_DNA_cc_intensity <- rep(0, mn*length(positions))
sd_DNA_cc_intensity <- rep(0, mn*length(positions))
mean_DNA_np_intensity <- rep(0, mn*length(positions))
sd_DNA_np_intensity <- rep(0, mn*length(positions))
mean_dCas9_intensity <- rep(0, mn*length(positions))
sd_dCas9_intensity <- rep(0, mn*length(positions))
mean_dCas9_cc_intensity <- rep(0, mn*length(positions))
sd_dCas9_cc_intensity <- rep(0, mn*length(positions))
mean_dCas9_np_intensity <- rep(0, mn*length(positions))
sd_dCas9_np_intensity <- rep(0, mn*length(positions))
mean_H3K9me3_intensity <- rep(0, mn*length(positions))
sd_H3K9me3_intensity <- rep(0, mn*length(positions))
mean_H3K9me3_cc_intensity <- rep(0, mn*length(positions))
sd_H3K9me3_cc_intensity <- rep(0, mn*length(positions))
mean_H3K9me3_np_intensity <- rep(0, mn*length(positions))
sd_H3K9me3_np_intensity <- rep(0, mn*length(positions))
mean_HP1_intensity <- rep(0, mn*length(positions))
sd_HP1_intensity <- rep(0, mn*length(positions))
mean_HP1_cc_intensity <- rep(0, mn*length(positions))
sd_HP1_cc_intensity <- rep(0, mn*length(positions))
mean_HP1_np_intensity <- rep(0, mn*length(positions))
sd_HP1_np_intensity <- rep(0, mn*length(positions))
area <- rep(0, mn*length(positions))
area_cc <- rep(0, mn*length(positions))
area_np <- rep(0, mn*length(positions))
pos <- rep(0, mn*length(positions))
num_cc <- rep(0, mn*length(positions))
cnt <- 1

# analyze nuclei
for(i in 1:length(positions)) {
  data <- readImage(paste0(folder,"/",positions[i]))
  
  # create image objects from data
  DNA_projection <- data[,,seq(3*numberOfFrames(data)/4+1, 4*numberOfFrames(data)/4, by = 1)]    # channel pos 4: DNA
  dCas9_projection <- data[,,seq(2*numberOfFrames(data)/4+1, 3*numberOfFrames(data)/4, by = 1)]   # channel pos 3: dCas9-GFP(-VPR)
  H3K9me3_projection <- data[,,seq(numberOfFrames(data)/4+1, 2*numberOfFrames(data)/4, by = 1)]   # channel pos 2: H3K9me3
  HP1_projection <- data[,,seq(1, numberOfFrames(data)/4, by = 1)]                                # channel pos 1: HP1
  
  # generate nucleus masks using the above-defined threshold offset and channel selection
  if(nucleus_seg_channel=="DNA"){nucmask <- makeNucMask(gblur(DNA_projection, sigma = 3), threshold, nucSize_cutoff, remove_bordercells)}
  if(nucleus_seg_channel=="dCas9"){nucmask <- makeNucMask(gblur(dCas9_projection, sigma = 3), threshold, nucSize_cutoff, remove_bordercells)}
  if(nucleus_seg_channel=="H3K9me3"){nucmask <- makeNucMask(gblur(H3K9me3_projection, sigma = 3), threshold, nucSize_cutoff, remove_bordercells)}
  if(nucleus_seg_channel=="HP1"){nucmask <- makeNucMask(gblur(HP1_projection, sigma = 3), threshold, nucSize_cutoff, remove_bordercells)}
  
  # if no nucleus could be found for this position:
  # create a dummy image with one squared "nucleus" of value 1 in the top left corner (40000 px area) and background intensity 0
  # this dummy_image will be processed as all other images, thereby maintaining results for all positions and marking out-of-focus positions for later filtering
  if(max(nucmask)==0){
    dummy_image <- Image(0,dim(DNA_projection)[1:2])
    dummy_image[50:249,50:249] <- 1
    nucmask <- makeNucMask(dummy_image, threshold, nucSize_cutoff, remove_bordercells)
  }
  
  # initialize variables for CC segmentation
  allcc <- matrix(0, nrow = dim(nucmask)[1], ncol = dim(nucmask)[2])
  allnp <- matrix(0, nrow = dim(nucmask)[1], ncol = dim(nucmask)[2])
  
  # for the above-identified nuclear masks, quantify intensities, SD and area in each mask throughout all channels
    rois <- as.numeric(names(table(nucmask)[2:length(table(nucmask))]))
    
    for(j in 1:length(rois)) {
      mean_DNA_intensity[cnt] <- mean(DNA_projection[nucmask==rois[j]])
      sd_DNA_intensity[cnt] <- sd(DNA_projection[nucmask==rois[j]])
      mean_dCas9_intensity[cnt] <- mean(dCas9_projection[nucmask==rois[j]])
      sd_dCas9_intensity[cnt] <- sd(dCas9_projection[nucmask==rois[j]])
      mean_H3K9me3_intensity[cnt] <- mean(H3K9me3_projection[nucmask==rois[j]])
      sd_H3K9me3_intensity[cnt] <- sd(H3K9me3_projection[nucmask==rois[j]])
      mean_HP1_intensity[cnt] <- mean(HP1_projection[nucmask==rois[j]])
      sd_HP1_intensity[cnt] <- sd(HP1_projection[nucmask==rois[j]])
      area[cnt] <- sum(nucmask==rois[j])
      pos[cnt] <- i
      
      # compute nucleus shape features and positions, export them to a file (they are required later for manual curation)
      nuc_shape_features <- cbind(i,j,computeFeatures.moment(nucmask==rois[j]),computeFeatures.shape(nucmask==rois[j]))
      write.table(nuc_shape_features, file = paste0(paste0(results_folder,"/",date_time,"_nuc_shape_features"),"/nuc_shape_features_pos_",i,"_cell",j), quote = F, sep = "\t", row.names = F)
      
      # identify if current nucleus is dCas9-"positive" (if yes, proceed with CC segmentation for this nucleus)
      if(mean(gblur(dCas9_projection, sigma = 1)[nucmask==rois[j]]) > mean_cutoff * median(gblur(dCas9_projection,sigma = 1))) {
        
        # generate CC masks
        mask <- nucmask
        mask[mask!=rois[j]] <- 0
        nuc <- gblur(dCas9_projection, sigma = 1)
        nuc[mask!=rois[j]] <- NA
        ccmask <- makeCCMask(nuc, median(nuc, na.rm = T) + sat_cutoff * (max(nuc, na.rm = T) - median(nuc, na.rm = T)), mask)
        ccmask[nuc < (median(nuc, na.rm = T) + sat_cutoff * (max(nuc, na.rm = T) - median(nuc, na.rm = T)))] <- 0
        num_cc[cnt] <- length(table(ccmask)) - 1
        ccmask[ccmask!=0] <- 1
        allcc[ccmask>0] <- 1
        
        # make mask for nucleoplasm in the current nucleus (= inverse of CC masks)
        npmask <- nuc<(median(nuc, na.rm = T) + sat_cutoff * (max(nuc, na.rm = T) - median(nuc, na.rm = T)))
        npmask[is.na(npmask)] <- 0
        allnp[npmask>0] <- 1
        
        # CC region and nucleoplasm quantification
        mean_DNA_cc_intensity[cnt] <- mean(DNA_projection[ccmask>0])
        sd_DNA_cc_intensity[cnt] <- sd(DNA_projection[ccmask>0])
        mean_DNA_np_intensity[cnt] <- mean(DNA_projection[npmask>0])
        sd_DNA_np_intensity[cnt] <- sd(DNA_projection[npmask>0])
        mean_dCas9_cc_intensity[cnt] <- mean(dCas9_projection[ccmask>0])
        sd_dCas9_cc_intensity[cnt] <- sd(dCas9_projection[ccmask>0])
        mean_dCas9_np_intensity[cnt] <- mean(dCas9_projection[npmask>0])
        sd_dCas9_np_intensity[cnt] <- sd(dCas9_projection[npmask>0])
        mean_H3K9me3_cc_intensity[cnt] <- mean(H3K9me3_projection[ccmask>0])
        sd_H3K9me3_cc_intensity[cnt] <- sd(H3K9me3_projection[ccmask>0])
        mean_H3K9me3_np_intensity[cnt] <- mean(H3K9me3_projection[npmask>0])
        sd_H3K9me3_np_intensity[cnt] <- sd(H3K9me3_projection[npmask>0])
        mean_HP1_cc_intensity[cnt] <- mean(HP1_projection[ccmask>0])
        sd_HP1_cc_intensity[cnt] <- sd(HP1_projection[ccmask>0])
        mean_HP1_np_intensity[cnt] <- mean(HP1_projection[npmask>0])
        sd_HP1_np_intensity[cnt] <- sd(HP1_projection[npmask>0])
        area_cc[cnt] <- sum(ccmask>0)
        area_np[cnt] <- sum(npmask>0)
      }
      # increase total nucleus counter
      cnt <- cnt + 1
    }
  
  # write out images with nucleus and CC masks overlayed onto DNA/dCas9/H3K9me3/HP1 for visualization & later curation
  # nucleus mask: yellow outlines
  # cc mask: blue outlines
  # np mask: opaque red area
  DNA_projection_with_mask <- toRGB(DNA_projection/max(DNA_projection))
  DNA_projection_with_mask <- paintObjects(nucmask, DNA_projection_with_mask, opac=c(1,0), col=c("yellow","yellow"))
  DNA_projection_with_mask <- paintObjects(allcc, DNA_projection_with_mask, opac=c(1,0), col=c("blue","blue"))
  DNA_projection_with_mask <- paintObjects(allnp, DNA_projection_with_mask, opac=c(0,0.2), col=c("red","red"))
  writeImage(combine(DNA_projection_with_mask, toRGB(DNA_projection/max(DNA_projection))), bits.per.sample=8, files=paste0(results_folder,"/",date_time,"_",gsub(".tif","_DNA_with_masks.tif",positions[i])))
  
  dCas9_projection_with_mask <- toRGB(dCas9_projection/max(dCas9_projection))
  dCas9_projection_with_mask <- paintObjects(nucmask, dCas9_projection_with_mask, opac=c(1,0), col=c("yellow","yellow"))
  dCas9_projection_with_mask <- paintObjects(allcc, dCas9_projection_with_mask, opac=c(1,0), col=c("blue","blue"))
  dCas9_projection_with_mask <- paintObjects(allnp, dCas9_projection_with_mask, opac=c(0,0.2), col=c("red","red"))
  writeImage(combine(dCas9_projection_with_mask, toRGB(dCas9_projection/max(dCas9_projection))), bits.per.sample=8, files=paste0(results_folder,"/",date_time,"_",gsub(".tif","_dCas9_with_masks.tif",positions[i])))
  
  H3K9me3_projection_with_mask <- toRGB(H3K9me3_projection/max(H3K9me3_projection))
  H3K9me3_projection_with_mask <- paintObjects(nucmask, H3K9me3_projection_with_mask, opac=c(1,0), col=c("yellow","yellow"))
  H3K9me3_projection_with_mask <- paintObjects(allcc, H3K9me3_projection_with_mask, opac=c(1,0), col=c("blue","blue"))
  H3K9me3_projection_with_mask <- paintObjects(allnp, H3K9me3_projection_with_mask, opac=c(0,0.2), col=c("red","red"))
  writeImage(combine(H3K9me3_projection_with_mask, toRGB(H3K9me3_projection/max(H3K9me3_projection))), bits.per.sample=8, files=paste0(results_folder,"/",date_time,"_",gsub(".tif","_H3K9me3_with_masks.tif",positions[i])))
  
  HP1_projection_with_mask <- toRGB(HP1_projection/max(HP1_projection))
  HP1_projection_with_mask <- paintObjects(nucmask, HP1_projection_with_mask, opac=c(1,0), col=c("yellow","yellow"))
  HP1_projection_with_mask <- paintObjects(allcc, HP1_projection_with_mask, opac=c(1,0), col=c("blue","blue"))
  HP1_projection_with_mask <- paintObjects(allnp, HP1_projection_with_mask, opac=c(0,0.2), col=c("red","red"))
  writeImage(combine(HP1_projection_with_mask, toRGB(HP1_projection/max(HP1_projection))), bits.per.sample=8, files=paste0(results_folder,"/",date_time,"_",gsub(".tif","_HP1_with_masks.tif",positions[i])))
  
  # introduce position-dependent cell numbering
  cell_number_per_pos <- vector()
  pos <- pos[1:length(grep(FALSE, pos[]==0))]
  for (t in 1:max(pos)) {
    cell_number_per_pos <- c(cell_number_per_pos, c(1:table(pos)[[t]]))
  }
}

# collect nucleus shape and position features
nuc_shape_folder <- paste0(results_folder, "/", list.files(results_folder)[grep("nuc_shape_features",list.files(results_folder))])
nuc_shape_files <- list.files(nuc_shape_folder)

nucleus_features <- data.frame()
for (b in 1:length(nuc_shape_files)) {
  nucleus_features[b, 1:length(colnames(read.table(file = paste0(nuc_shape_folder,"/",nuc_shape_files[b]),header=T)))] <- read.table(file = paste0(nuc_shape_folder, "/", nuc_shape_files[b]), header = T,)
}

# re-order nucleus_features according to columns position and cell (i and j) and output to results file
nucleus_features <- nucleus_features[order(nucleus_features$i, nucleus_features$j),]
colnames(nucleus_features)[1] <- "position"
colnames(nucleus_features)[2] <- "cell"
write.table(file = paste0(results_folder,"/",date_time,"_nucleus_shape_features.csv"),nucleus_features, quote = F, sep = "\t", row.names = F)
if (file.exists(nuc_shape_folder)){unlink(nuc_shape_folder, recursive = T)}

# write position/cell indices,  intensity and area results into result file
res <- cbind(pos[1:(cnt-1)],cell_number_per_pos,mean_DNA_intensity[1:(cnt-1)],sd_DNA_intensity[1:(cnt-1)],mean_H3K9me3_intensity[1:(cnt-1)],sd_H3K9me3_intensity[1:(cnt-1)],mean_dCas9_intensity[1:(cnt-1)],sd_dCas9_intensity[1:(cnt-1)],mean_HP1_intensity[1:(cnt-1)],sd_HP1_intensity[1:(cnt-1)],area[1:(cnt-1)],mean_DNA_cc_intensity[1:(cnt-1)],sd_DNA_cc_intensity[1:(cnt-1)],mean_DNA_np_intensity[1:(cnt-1)],sd_DNA_np_intensity[1:(cnt-1)],mean_H3K9me3_cc_intensity[1:(cnt-1)],sd_H3K9me3_cc_intensity[1:(cnt-1)],mean_H3K9me3_np_intensity[1:(cnt-1)],sd_H3K9me3_np_intensity[1:(cnt-1)],mean_dCas9_cc_intensity[1:(cnt-1)],sd_dCas9_cc_intensity[1:(cnt-1)],mean_dCas9_np_intensity[1:(cnt-1)],sd_dCas9_np_intensity[1:(cnt-1)],mean_HP1_cc_intensity[1:(cnt-1)],sd_HP1_cc_intensity[1:(cnt-1)],mean_HP1_np_intensity[1:(cnt-1)],sd_HP1_np_intensity[1:(cnt-1)],area_cc[1:(cnt-1)],area_np[1:(cnt-1)],num_cc[1:(cnt-1)])
colnames(res) <- c("position","cell", "mean_DNA", "sd_DNA", "mean_H3K9me3", "sd_H3K9me3", "mean_dCas9", "sd_dCas9", "mean_HP1", "sd_HP1", "area", "mean_DNA_cc", "sd_DNA_cc", "mean_DNA_np", "sd_DNA_np", "mean_H3K9me3_cc", "sd_H3K9me3_cc", "mean_H3K9me3_np", "sd_H3K9me3_np", "mean_dCas9_cc", "sd_dCas9_cc", "mean_dCas9_np", "sd_dCas9_np", "mean_HP1_cc", "sd_HP1_cc", "mean_HP1_np", "sd_HP1_np", "area_cc", "area_np", "num_cc")
write.table(file=paste0(results_folder,"/",date_time, "_mean_sd_projection.csv"), res, sep = "\t", append = F, row.names = F, quote = F)
