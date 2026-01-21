makeCCMask <- function(im, threshold, nuc_mask, dilate_disc=3, erode_disc=3, ignore_boundary=3, minsize=30) {
  # find regions of interest
  mask <- im>threshold
  mask <- dilate(mask, makeBrush(dilate_disc, shape='disc'))
  mask <- erode(mask, makeBrush(erode_disc, shape='disc'))
  mask[erode(nuc_mask, makeBrush(ignore_boundary, shape='disc'))==0] <- 0
  mask <- bwlabel(mask) # find connected sets of pixels
  mask <- fillHull(mask) # fill holes in the mask
  
  # only keep large regions of interest
  regions <- table(mask)
  largest_region <- which(regions==max(regions))
  small_regions <- which(regions<minsize)-1
  mask[mask %in% small_regions] <- 0
  mask[mask==as.numeric(names(largest_region)[1])] <- 0
  #mask[mask!=0] <- 1
  
  # return mask
  return(mask)
}