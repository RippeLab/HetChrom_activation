###############################################################################################
################                          Introduction                           ##############
###############################################################################################

# This script combines the spotcounts from RS-FISH output with the output table of curate_CC

###############################################################################################

# specify parent directory of input files
inputroot <- "insert file path/dn/p65/"

# specify path to FISH results.csv
FISHinput <- paste0(inputroot,"bg_and_ff_corrected/results.csv")
# specify path to ..._curated.csv
CUREinput <- paste0(inputroot,"xxxxxxxxxxxxxx_segmentCC_results/xxxxxxxxxxxxxx_curateCC_results/mean_sd_projection_curated.csv")

# import csvs
FISH <- read.csv(FISHinput,sep="\t")
CURE <- read.csv(CUREinput,sep="\t")

# combine spotcounts with curate.csv if the first two columns are identical (position + cell)

# Check if the first two columns are identical
if(all(FISH[,1:2] == CURE[,1:2])) {
  # If they are identical, then combine the dataframes
  combined <- cbind(CURE,FISH[,3])
  # Change the name of the last column to "spotcount"
  colnames(combined)[ncol(combined)] <- "spotcount"
} else {
  print("The first two columns are not identical.")
}

write.table(combined,paste0(inputroot,"combined.csv"),sep = "\t",row.names = FALSE)

