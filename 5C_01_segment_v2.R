####################################################################################################################
############################################### Overview ###########################################################
####################################################################################################################

# thresholding-based segmentation strategies to generate nucleus (Nuc), chromocenter (CC) and nucleoplasm masks. Intensities and pixel area are quantified for the masks to extract relevant image features.
# Use: initially created for analysis of iMEF VPR recruitment to CC + IF (H3K27ac + p300)

####################################################################################################################
############################################### ADJUSTABLE PARAMETERS ##############################################
####################################################################################################################

# load functions
source("insert file path/functions/makeCCMask.R")
source("insert file path/functions/makeNucMask.R")

# specify location of input (maximum-intensity projection .tif files)
folder <- "insert file path/dn/VP16/"

# adjustable adaptive thresholding offset for nucleus segmentation
# (default = 0.0008)
threshold <- 0.0008

# select the channel in which the nuclei will be segmented
# ("dapi", "dCas9", "rna", "ac", "p300")
nucleus_seg_channel <- "dapi"

# remove incomplete nuclei at the image borders?
# ("y" = yes, "n" = no)
remove_bordercells <- "y"

# thresholding parameter for CC segmentation in the dCas9 channel using the makeCCMask function
# intensity threshold above which a pixel is considered to belong to a chromocenter (based on dCas9 enrichment) is given by threshold = median(nuclear intensity) + sat_cutoff * (maximum(nuclear intensity) - median(nuclear intensity))
# (default = 0.1)
sat_cutoff <- 0.15

# CC segmentation will only be carried out on cells that express dCas9 ("positive cells")
# cell is considered "positive" if its mean nuclear dCas9 intensity is "mean_cutoff"-fold above the median (background) of the total image intensity
# this parameter heavily depends on the type of construct used and the signal-to-background ratio
# (default = 1.5)
mean_cutoff <- 1.5

# size-filter for nucleus segmentation (pixels); needs to be adjusted to cell size in sample
# (default = 8000)
nucSize_cutoff <- 8000

# maximum number of nuclei expected per image; required for variable initialization
# (default = 20)
mn <- 20

####################################################################################################################
#################################################### FIXED #########################################################
####################################################################################################################
# load libraries
library(EBImage)
library(abind)

# create timestamp for this run
timestamp <- format(Sys.time(),"%Y%m%d%H%M%S")

# create results folders
dir.create(paste0(folder,timestamp,"_segmentCC_results"))
results_folder <- paste0(folder,timestamp,"_segmentCC_results")
dir.create(paste0(results_folder,"/nuc_shape_features"))

# document adjustable segmentation parameters used in this run
parameters_df <- data.frame(threshold, nucleus_seg_channel, remove_bordercells, sat_cutoff, mean_cutoff, nucSize_cutoff, mn)
colnames(parameters_df) <- c("threshold","nucleus_seg_channel","remove_bordercells","sat_cutoff","mean_cutoff","nucSize_cutoff","mn")
write.table(parameters_df, file = paste0(results_folder,"/segmentation_parameters.csv"), quote = F, sep = "\t", row.names = F)

# generate and write a position list to document the files that were processed in this run
positions <- list.files(path = folder, pattern = "tif")
write.table(cbind(seq(1,length(positions), by = 1), positions), file = paste0(results_folder,"/positions.csv"), quote = F, sep = "\t", row.names = F)

# initialization of variables and empty vectors for storing quantification results
mean_dapi_intensity <- rep(0, mn*length(positions))
sd_dapi_intensity <- rep(0, mn*length(positions))
mean_dapi_cc_intensity <- rep(0, mn*length(positions))
sd_dapi_cc_intensity <- rep(0, mn*length(positions))
mean_dapi_np_intensity <- rep(0, mn*length(positions))
sd_dapi_np_intensity <- rep(0, mn*length(positions))
mean_dcas9_intensity <- rep(0, mn*length(positions))
sd_dcas9_intensity <- rep(0, mn*length(positions))
mean_dcas9_cc_intensity <- rep(0, mn*length(positions))
sd_dcas9_cc_intensity <- rep(0, mn*length(positions))
mean_dcas9_np_intensity <- rep(0, mn*length(positions))
sd_dcas9_np_intensity <- rep(0, mn*length(positions))
mean_rna_intensity <- rep(0, mn*length(positions))
sd_rna_intensity <- rep(0, mn*length(positions))
mean_rna_cc_intensity <- rep(0, mn*length(positions))
sd_rna_cc_intensity <- rep(0, mn*length(positions))
mean_rna_np_intensity <- rep(0, mn*length(positions))
sd_rna_np_intensity <- rep(0, mn*length(positions))
mean_ac_intensity <- rep(0, mn*length(positions))
sd_ac_intensity <- rep(0, mn*length(positions))
mean_ac_cc_intensity <- rep(0, mn*length(positions))
sd_ac_cc_intensity <- rep(0, mn*length(positions))
mean_ac_np_intensity <- rep(0, mn*length(positions))
sd_ac_np_intensity <- rep(0, mn*length(positions))
mean_p300_intensity <- rep(0, mn*length(positions))
sd_p300_intensity <- rep(0, mn*length(positions))
mean_p300_cc_intensity <- rep(0, mn*length(positions))
sd_p300_cc_intensity <- rep(0, mn*length(positions))
mean_p300_np_intensity <- rep(0, mn*length(positions))
sd_p300_np_intensity <- rep(0, mn*length(positions))
area <- rep(0, mn*length(positions))
area_cc <- rep(0, mn*length(positions))
area_np <- rep(0, mn*length(positions))
pos <- rep(0, mn*length(positions))
num_cc <- rep(0, mn*length(positions))
cnt <- 1

# analyze nuclei
i <- 1
for(i in 1:length(positions)) {
  data <- readImage(paste0(folder,positions[i]))
  
  # create image objects from data
  dapi_projection <- data[,,seq(4*numberOfFrames(data)/5+1, 5*numberOfFrames(data)/5, by = 1)]    # channel pos 5: DAPI
  dcas9_projection <- data[,,seq(3*numberOfFrames(data)/5+1, 4*numberOfFrames(data)/5, by = 1)]   # channel pos 4: dCas9-GFP-VPR
  rna_projection <- data[,,seq(2*numberOfFrames(data)/5+1, 3*numberOfFrames(data)/5, by = 1)]     # channel pos 3: RNAScope
  ac_projection <- data[,,seq(numberOfFrames(data)/5+1, 2*numberOfFrames(data)/5, by = 1)]        # channel pos 2: H3K27ac
  p300_projection <- data[,,seq(1, numberOfFrames(data)/5, by = 1)]                               # channel pos 1: p300
  
  # generate nucleus masks using the above-defined threshold offset and channel selection
  if(nucleus_seg_channel=="dapi"){nucmask <- makeNucMask(gblur(dapi_projection, sigma = 3), threshold, nucSize_cutoff, remove_bordercells)}
  if(nucleus_seg_channel=="dcas9"){nucmask <- makeNucMask(gblur(dcas9_projection, sigma = 3), threshold, nucSize_cutoff, remove_bordercells)}
  if(nucleus_seg_channel=="rna"){nucmask <- makeNucMask(gblur(rna_projection, sigma = 3), threshold, nucSize_cutoff, remove_bordercells)}
  if(nucleus_seg_channel=="ac"){nucmask <- makeNucMask(gblur(ac_projection, sigma = 3), threshold, nucSize_cutoff, remove_bordercells)}
  if(nucleus_seg_channel=="p300"){nucmask <- makeNucMask(gblur(p300_projection, sigma = 3), threshold, nucSize_cutoff, remove_bordercells)}
  
  # if no nucleus could be found for this position:
  # create a dummy image with one squared "nucleus" of value 1 in the top left corner (40000 px area) and background intensity 0
  # this dummy_image will be processed as all other images, thereby maintaining results for all positions and marking out-of-focus positions for later filtering
  if(max(nucmask)==0){
    dummy_image <- Image(0,dim(dapi_projection)[1:2])
    dummy_image[50:249,50:249] <- 1
    nucmask <- makeNucMask(dummy_image, threshold, nucSize_cutoff, remove_bordercells)
  }
  
  # initialize variables for CC segmentation
  allcc <- matrix(0, nrow = dim(nucmask)[1], ncol = dim(nucmask)[2])
  allnp <- matrix(0, nrow = dim(nucmask)[1], ncol = dim(nucmask)[2])
  
  # for the above-identified nuclear masks, quantify intensities, SD and area in each mask throughout all channels
  rois <- as.numeric(names(table(nucmask)[2:length(table(nucmask))]))
  
  j <- 1
  for(j in 1:length(rois)) {
    mean_dapi_intensity[cnt] <- mean(dapi_projection[nucmask==rois[j]])
    sd_dapi_intensity[cnt] <- sd(dapi_projection[nucmask==rois[j]])
    mean_dcas9_intensity[cnt] <- mean(dcas9_projection[nucmask==rois[j]])
    sd_dcas9_intensity[cnt] <- sd(dcas9_projection[nucmask==rois[j]])
    mean_rna_intensity[cnt] <- mean(rna_projection[nucmask==rois[j]])
    sd_rna_intensity[cnt] <- sd(rna_projection[nucmask==rois[j]])
    mean_ac_intensity[cnt] <- mean(ac_projection[nucmask==rois[j]])
    sd_ac_intensity[cnt] <- sd(ac_projection[nucmask==rois[j]])
    mean_p300_intensity[cnt] <- mean(p300_projection[nucmask==rois[j]])
    sd_p300_intensity[cnt] <- sd(p300_projection[nucmask==rois[j]])
    area[cnt] <- sum(nucmask==rois[j])
    pos[cnt] <- i
    
    # compute nucleus shape features and positions, export them to a file (they are required later for manual curation)
    nuc_shape_features <- cbind(i,j,computeFeatures.moment(nucmask==rois[j]),computeFeatures.shape(nucmask==rois[j]))
    write.table(nuc_shape_features, file = paste0(paste0(results_folder,"/nuc_shape_features"),"/nuc_shape_features_pos_",i,"_cell",j), quote = F, sep = "\t", row.names = F)
    
    # identify if current nucleus is dcas9-"positive" (if yes, proceed with CC segmentation for this nucleus)
    if(mean(gblur(dcas9_projection, sigma = 1)[nucmask==rois[j]]) > mean_cutoff * median(gblur(dcas9_projection,sigma = 1))) {
      
      # generate CC masks
      mask <- nucmask
      mask[mask!=rois[j]] <- 0
      nuc <- gblur(dcas9_projection, sigma = 1)
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
      mean_dapi_cc_intensity[cnt] <- mean(dapi_projection[ccmask>0])
      sd_dapi_cc_intensity[cnt] <- sd(dapi_projection[ccmask>0])
      mean_dapi_np_intensity[cnt] <- mean(dapi_projection[npmask>0])
      sd_dapi_np_intensity[cnt] <- sd(dapi_projection[npmask>0])
      mean_dcas9_cc_intensity[cnt] <- mean(dcas9_projection[ccmask>0])
      sd_dcas9_cc_intensity[cnt] <- sd(dcas9_projection[ccmask>0])
      mean_dcas9_np_intensity[cnt] <- mean(dcas9_projection[npmask>0])
      sd_dcas9_np_intensity[cnt] <- sd(dcas9_projection[npmask>0])
      mean_rna_cc_intensity[cnt] <- mean(rna_projection[ccmask>0])
      sd_rna_cc_intensity[cnt] <- sd(rna_projection[ccmask>0])
      mean_rna_np_intensity[cnt] <- mean(rna_projection[npmask>0])
      sd_rna_np_intensity[cnt] <- sd(rna_projection[npmask>0])
      mean_ac_cc_intensity[cnt] <- mean(ac_projection[ccmask>0])
      sd_ac_cc_intensity[cnt] <- sd(ac_projection[ccmask>0])
      mean_ac_np_intensity[cnt] <- mean(ac_projection[npmask>0])
      sd_ac_np_intensity[cnt] <- sd(ac_projection[npmask>0])
      mean_p300_cc_intensity[cnt] <- mean(p300_projection[ccmask>0])
      sd_p300_cc_intensity[cnt] <- sd(p300_projection[ccmask>0])
      mean_p300_np_intensity[cnt] <- mean(p300_projection[npmask>0])
      sd_p300_np_intensity[cnt] <- sd(p300_projection[npmask>0])
      area_cc[cnt] <- sum(ccmask>0)
      area_np[cnt] <- sum(npmask>0)
    }
    # increase total nucleus counter
    cnt <- cnt + 1
  }
  
  # write out images with nucleus and CC masks overlayed onto DAPI/dCas9/ac/p300 for visualization & later curation
  # nucleus mask: yellow outlines
  # cc mask: blue outlines
  # np mask: opaque red area
  dapi_projection_with_mask <- toRGB(dapi_projection/max(dapi_projection))
  dapi_projection_with_mask <- paintObjects(nucmask, dapi_projection_with_mask, opac=c(1,0), col=c("yellow","yellow"))
  dapi_projection_with_mask <- paintObjects(allcc, dapi_projection_with_mask, opac=c(1,0), col=c("blue","blue"))
  dapi_projection_with_mask <- paintObjects(allnp, dapi_projection_with_mask, opac=c(0,0.2), col=c("red","red"))
  writeImage(combine(dapi_projection_with_mask, toRGB(dapi_projection/max(dapi_projection))), bits.per.sample=8, files=paste0(results_folder,"/",gsub(".tif","_dapi_with_masks.tif",positions[i])))
  
  dcas9_projection_with_mask <- toRGB(dcas9_projection/max(dcas9_projection))
  dcas9_projection_with_mask <- paintObjects(nucmask, dcas9_projection_with_mask, opac=c(1,0), col=c("yellow","yellow"))
  dcas9_projection_with_mask <- paintObjects(allcc, dcas9_projection_with_mask, opac=c(1,0), col=c("blue","blue"))
  dcas9_projection_with_mask <- paintObjects(allnp, dcas9_projection_with_mask, opac=c(0,0.2), col=c("red","red"))
  writeImage(combine(dcas9_projection_with_mask, toRGB(dcas9_projection/max(dcas9_projection))), bits.per.sample=8, files=paste0(results_folder,"/",gsub(".tif","_dcas9_with_masks.tif",positions[i])))
  
  rna_projection_with_mask <- toRGB(rna_projection/max(rna_projection))
  rna_projection_with_mask <- paintObjects(nucmask, rna_projection_with_mask, opac=c(1,0), col=c("yellow","yellow"))
  rna_projection_with_mask <- paintObjects(allcc, rna_projection_with_mask, opac=c(1,0), col=c("blue","blue"))
  rna_projection_with_mask <- paintObjects(allnp, rna_projection_with_mask, opac=c(0,0.2), col=c("red","red"))
  writeImage(combine(rna_projection_with_mask, toRGB(rna_projection/max(rna_projection))), bits.per.sample=8, files=paste0(results_folder,"/",gsub(".tif","_rna_with_masks.tif",positions[i])))
  
  ac_projection_with_mask <- toRGB(ac_projection/max(ac_projection))
  ac_projection_with_mask <- paintObjects(nucmask, ac_projection_with_mask, opac=c(1,0), col=c("yellow","yellow"))
  ac_projection_with_mask <- paintObjects(allcc, ac_projection_with_mask, opac=c(1,0), col=c("blue","blue"))
  ac_projection_with_mask <- paintObjects(allnp, ac_projection_with_mask, opac=c(0,0.2), col=c("red","red"))
  writeImage(combine(ac_projection_with_mask, toRGB(ac_projection/max(ac_projection))), bits.per.sample=8, files=paste0(results_folder,"/",gsub(".tif","_ac_with_masks.tif",positions[i])))
  
  p300_projection_with_mask <- toRGB(p300_projection/max(p300_projection))
  p300_projection_with_mask <- paintObjects(nucmask, p300_projection_with_mask, opac=c(1,0), col=c("yellow","yellow"))
  p300_projection_with_mask <- paintObjects(allcc, p300_projection_with_mask, opac=c(1,0), col=c("blue","blue"))
  p300_projection_with_mask <- paintObjects(allnp, p300_projection_with_mask, opac=c(0,0.2), col=c("red","red"))
  writeImage(combine(p300_projection_with_mask, toRGB(p300_projection/max(p300_projection))), bits.per.sample=8, files=paste0(results_folder,"/",gsub(".tif","_p300_with_masks.tif",positions[i])))
  
  # write out nucmask images
  writeImage((nucmask/2^16)+0.001, bits.per.sample = 16, files=paste0(results_folder,"/",gsub(".tif","_nucmask.tif", positions[i]))) # bits.per.sample only works for tiff files so doesn't do anything here !
  
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
write.table(file = paste0(results_folder,"/nucleus_shape_features.csv"),nucleus_features, quote = F, sep = "\t", row.names = F)
if (file.exists(nuc_shape_folder)){unlink(nuc_shape_folder, recursive = T)}

# write position/cell indices,  intensity and area results into result file
res <- cbind(pos[1:(cnt-1)],cell_number_per_pos,mean_dapi_intensity[1:(cnt-1)],sd_dapi_intensity[1:(cnt-1)],mean_dcas9_intensity[1:(cnt-1)],sd_dcas9_intensity[1:(cnt-1)],mean_rna_intensity[1:(cnt-1)],sd_rna_intensity[1:(cnt-1)],mean_ac_intensity[1:(cnt-1)],sd_ac_intensity[1:(cnt-1)],mean_p300_intensity[1:(cnt-1)],sd_p300_intensity[1:(cnt-1)],area[1:(cnt-1)],mean_dapi_cc_intensity[1:(cnt-1)],sd_dapi_cc_intensity[1:(cnt-1)],mean_dapi_np_intensity[1:(cnt-1)],sd_dapi_np_intensity[1:(cnt-1)],mean_dcas9_cc_intensity[1:(cnt-1)],sd_dcas9_cc_intensity[1:(cnt-1)],mean_dcas9_np_intensity[1:(cnt-1)],sd_dcas9_np_intensity[1:(cnt-1)],mean_rna_cc_intensity[1:(cnt-1)],sd_rna_cc_intensity[1:(cnt-1)],mean_rna_np_intensity[1:(cnt-1)],sd_rna_np_intensity[1:(cnt-1)],mean_ac_cc_intensity[1:(cnt-1)],sd_ac_cc_intensity[1:(cnt-1)],mean_ac_np_intensity[1:(cnt-1)],sd_ac_np_intensity[1:(cnt-1)],mean_p300_cc_intensity[1:(cnt-1)],sd_p300_cc_intensity[1:(cnt-1)],mean_p300_np_intensity[1:(cnt-1)],sd_p300_np_intensity[1:(cnt-1)],area_cc[1:(cnt-1)],area_np[1:(cnt-1)],num_cc[1:(cnt-1)])
colnames(res) <- c("position","cell", "mean_dapi", "sd_dapi", "mean_dcas9", "sd_dcas9","mean_rna", "sd_rna","mean_ac", "sd_ac", "mean_p300", "sd_p300", "area", "mean_dapi_cc", "sd_dapi_cc", "mean_dapi_np", "sd_dapi_np", "mean_dcas9_cc", "sd_dcas9_cc", "mean_dcas9_np", "sd_dcas9_np", "mean_rna_cc", "sd_rna_cc", "mean_rna_np", "sd_rna_np", "mean_ac_cc", "sd_ac_cc", "mean_ac_np", "sd_ac_np", "mean_p300_cc", "sd_p300_cc", "mean_p300_np", "sd_p300_np", "area_cc", "area_np", "num_cc")
write.table(file=paste0(results_folder,"/mean_sd_projection.csv"), res, sep = "\t", append = F, row.names = F, quote = F)
