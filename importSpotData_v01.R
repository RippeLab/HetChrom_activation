#' This function uses results tables from RS-FISH (FIJI plugin) containing information on spot position (x,y) and intensity
#' to generate an image object with spots drawn as single pixels, wherein each pixel is assigned the spot intensity as value.
#' Background has the value 0.
#' @param spotData Results table from RS-FISH (imported by read.csv())
#' @param nuc_masks Corresponding image of nuclear masks (bg = 0, nuclei = 1 to n)
#' @returns Image object with drawn spots 

importSpotData <- function(spotData,nuc_masks,norm16bit=TRUE){
  # Create empty image 
  img <- Image(0,c(dim(nuc_masks)[1:2]))
  # Get xy coordinates from the table - rounded so the values are pixels
  spots_xy <- round(spotData[,1:2])
  if(norm16bit==FALSE){
    #set all the spots to the value of their intensity
    c <- 1
    while (c <= nrow(spots_xy)){
      img[spots_xy[c,"x"], spots_xy[c,"y"]] <- round(spotData[c,"intensity"])
      c = c+1}
    return(img)
  }else{
    # set all the spots to the value of their intensity
    # divide by 2^16 so that values in the final image have a range from 0 to 1 
    c <- 1
    while (c <= nrow(spots_xy)){
      img[spots_xy[c,"x"]+1, spots_xy[c,"y"]+1] <- round(spotData[c,"intensity"])/2^16
      c = c+1}
    # Note that the +1 has to be added because EBImage objects start with (1,1) top left
    # while the output from RS-FISH assumes is based on a (0,0) origin
    return(img)
  }
}
