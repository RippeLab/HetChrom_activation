############################################ Curate CC #############################################################
# Displays segmentation output from segmentCC for manual curation of segmentation results
# Adds a column to the mean_sd_projection.csv results file with user annotation of single cells
# and creates the new results file mean_sd_projection_curated.csv
# Cells with satisfactory segmentation results are manually annotated as "g" (good) or "b"(bad)
# Cells without segmentation are annotated with "n" (NA)

# load libraries and functions
library(EBImage)
library(abind)

####################################################################################################################
############################################### ADJUSTABLE PARAMETERS ##############################################
####################################################################################################################

# define input folder which corresponds to the segmentCC results folder
input_folder <- "insert file path/wt/p300core/xxxx-xx-xx_xx-xx-xx_segmentCC_results/"

####################################################################################################################
############################################### FIXED ##############################################################
# create timestamp variable for this run
date_time <- gsub(" ","_",gsub(":","-",substr(Sys.time(),1,19)))

# retrieve processed image positions from file
positions <- read.table(file=paste0(input_folder,list.files(input_folder)[grep("positions",list.files(input_folder))]))

# retrieve corresponding segmented image overlays to be displayed
DNA_overlays <-list.files(input_folder)[grep("DNA_with_masks",list.files(input_folder))]
dCas9_overlays <-list.files(input_folder)[grep("dCas9_with_masks",list.files(input_folder))]

#check if segmentCC run was complete
if (length(DNA_overlays)!=length(rownames(positions))){print("WARNING: It seems that not all of your images have been segmented. This can happen if a run of segmentCC was aborted. Please check if every image has a corresponding _DNA_with_masks.tif processed file.")}
if (length(list.files(input_folder, pattern=".csv"))<4){print("ERROR: It seems that not all results files have been created during the run of segmentCC. Please make sure that all neccessary .csv files are available in the results folder.")}

# retrieve nucleus coordinate & shape features from file
nuc_shape_features <- read.table(file=paste0(input_folder,list.files(input_folder)[grep("nucleus_shape_features",list.files(input_folder))]), header=T)

# retrieve area & intensity quantitifaction data from file
mean_sd_projection <- read.table(file=paste0(input_folder,list.files(input_folder)[grep("mean_sd_projection",list.files(input_folder))]), header=T)

# identify cells for which chromocenters have been segmented and reduce data frame accordingly
# also removes data corresponding to dummy nuclei of 40000 px area (positions where no nuclei have been found)
mean_sd_projection_seg <- mean_sd_projection[which(mean_sd_projection$mean_DNA_cc!=0 & mean_sd_projection$area!=40000),]

# collect nuc_shape_features for segmented cells only
nuc_shape_features_seg <- data.frame()
for (l in 1: length(mean_sd_projection_seg$position)){
  nuc_shape_features_seg <- rbind(nuc_shape_features_seg,nuc_shape_features[which(nuc_shape_features$position==mean_sd_projection_seg$position[l] & nuc_shape_features$cell==mean_sd_projection_seg$cell[l]),])
}

# inititalize annotation variable
annotation<-character()

# loop over positions for annotation
for (i in 1:length(rownames(positions))){
  dispIm1 <- readImage(file=paste0(input_folder,DNA_overlays[i]))
  dispIm2 <- readImage(file=paste0(input_folder,dCas9_overlays[i]))
  dispIm<-combine(dispIm1^0.5,dispIm2^0.5)
  display(dispIm,method = 'raster',all = T)
  
  # get number of segmented cells for the current position
  # skip if position contains zero segmented cells
  j <- NULL
  if(length(table(nuc_shape_features_seg$position)[names(table(nuc_shape_features_seg$position))==i])==0){
    j <- 0
  }else{
    j<-table(nuc_shape_features_seg$position)[names(table(nuc_shape_features_seg$position))==i][[1]]
  }
  if(j==0){print(paste0("NOTE:position ( ",i, " ) was skipped because it contained no segmented cells.")) 
    next()
  } 
  
  # get user annotation of cells
  flawless_image <- NULL
  flawless_image <- readline('Take segmentation as is? yes: y, no:n    ')
  if (flawless_image!="y"){
    # loop over segmented cells of each position and highlight for annotation
    for(k in 1:j){
      nuc_coord <- c(as.integer(nuc_shape_features_seg[nuc_shape_features_seg$position==i,][k,]$m.cx),nuc_shape_features_seg[nuc_shape_features_seg$position==i,][k,]$m.cy)
      display(drawCircle(dispIm,nuc_coord[1],nuc_coord[2],10,"red",fill=T,z=1),method = 'raster',all = T)
      annotation <- c(annotation, readline('Evaluate segmentation for the highlighted cell as good(g) or bad(b) or good+ring shape(r)     '))
    }
  }
  else{
    annotation <- c(annotation,rep("g",j))
  }
} 

# complete annotation vector
annotation_full <- character()
annotation_full[1:length(mean_sd_projection$position)] <- c(rep("n",length(mean_sd_projection$position)))
annotation_ids <- c(as.numeric(rownames(mean_sd_projection_seg)))

for(p in 1:length(rownames(mean_sd_projection_seg))){
  annotation_full[annotation_ids[p]] <- annotation[p]
}
mean_sd_projection_curated <- cbind(mean_sd_projection,annotation_full)
colnames(mean_sd_projection_curated)[31] <- "annotation"

# write results into mean_sd_projections_curated.csv
dir.create(paste0(input_folder,date_time,"_curateCC_results"))
write.table(file=paste0(input_folder,date_time,"_curateCC_results","/","mean_sd_projection_curated.csv"), mean_sd_projection_curated, sep = "\t", append = F, row.names = F, quote = F)
