####################################################################################################################
############################################### Overview ###########################################################
####################################################################################################################

# Creates violinplots from a curated segmentation dataset generated using curateCC and exports them as PDF in ...
# ... _plotCC within the curateCC_results folder (adapted to images with 5 channels/readout)

# User can define to filter for good/bad segmentation quality of cells for plotting

####################################################################################################################
############################################### ADJUSTABLE PARAMETERS ##############################################
####################################################################################################################

# specify folder location of the curateCC results
folder_curated <- "insert file path/combined_plots/"

# select cell quality filtering 
# "g" (include good cells only), "b" (include bad cells only), "gb" (include both)
cell_filter <- "g"

# list of all cell lines
celllines <- c("wt","dn")

# list of all conditions (in preferred order, P is the plasmid declaration according to the CloneX ID)
condsP <- c("1299","1502","1618","1473")
conds <- c("GFP","VP16","p65","VPR")

# list of all readouts (in preferred order)
reots <- c("dna","dcas9","rna","ac","p300")

####################################################################################################################
############################################### FIXED ##############################################################
####################################################################################################################
# set working directory to folder containing the individual .csv files
setwd(folder_curated)

# load libraries
library(ggplot2)
library(ggsignif)
library(dplyr)

# create timestamp for this run
timestamp <- format(Sys.time(),"%Y%m%d%H%M%S")

# create outout folder
dir.create(paste0(folder_curated,timestamp,"_plotCC_PDFs"))

# generate list of all combinations of cell lines + condition
CLP <- do.call(paste0,expand.grid(celllines,"_",condsP))

# generate list of all combinations of cell lines + condition
group_simple <- unlist(lapply(conds,function(cond) {paste0(celllines,"_",cond)}))

# generate list of all combinations of cell lines + condition + readouts
group_complex <- unlist(lapply(conds,function(cond) {unlist(lapply(reots,function(reot) {paste0(celllines,"_",cond,"_",reot)}))}))

# read the file containing cellline and condition in the wd folder
for (i in 1:length(CLP)) {
  file_name <- list.files(folder_curated,pattern=CLP[i])
  temp <- read.table(file=paste0(folder_curated,file_name),header=T,stringsAsFactors=F)
  # filter table according to cell filter
  if (cell_filter=="gb"){
    temp_filtered <- temp[temp[,"annotation"]=="g"|temp[,"annotation"]=="b",]
    if(length(rownames(temp_filtered))==0){stop("No cells with specified quality annotation were found. Check filter settings and whether the dataset contains cells with this type of annotation.")}
  }else{
    temp_filtered <- temp[temp[,"annotation"]==cell_filter,]
    if(length(rownames(temp_filtered))==0){stop("No cells with specified quality annotation were found. Check filter settings and whether the dataset contains cells with this type of annotation.")}}
  # replace dapi in the column name with dna (important for later + consistency)
  new_colnames <- gsub("dapi","dna",colnames(temp_filtered))
  colnames(temp_filtered) <- new_colnames
  # assign filtered table to a new name
  assign(group_simple[i],as.data.frame(temp_filtered))}

# cleanup
rm(condsP,CLP,i,file_name,temp,new_colnames,temp_filtered)


####################################################################################################################
############### create "simple" dataframe for plotting #############################################################
####################################################################################################################
# create the background column dynamically
counter <- 1
background <- c()
for (i in 1:length(conds)) {
  for (j in 1:length(celllines)) {
    background <- c(background,rep(celllines[j],nrow(get(group_simple[counter]))))
    counter <- counter+1}}

# create the condition column dynamically
condition <- c()
for (i in 1:length(conds)) {condition <- c(condition,rep(conds[i],nrow(get(group_simple[2*i-1]))+nrow(get(group_simple[2*i]))))}

# create the group column dynamically
group <- unlist(lapply(group_simple,function(g) rep(g,nrow(get(g)))))

# create readout columns dynamically
area <- unlist(sapply(group_simple,function(g) get(g)[,"area"]),use.names=F)
decond <- unlist(sapply(group_simple,function(g) get(g)[,"area_cc"] / get(g)[,"area"] * 100),use.names=F)
H3K27ac <- unlist(sapply(group_simple,function(g) get(g)[,"mean_ac_cc"] / get(g)[,"mean_ac_np"]),use.names=F)
spots <- unlist(sapply(group_simple,function(g) get(g)[,"spotcount"]),use.names=F)
mean_rna <- unlist(sapply(group_simple,function(g) get(g)[,"mean_rna"]),use.names=F)

# create "simple" dataframe
simple_df <- data.frame(background,condition,group,area,decond,H3K27ac,spots,mean_rna)

# create a custom factor variable with the desired order for "background"
background_order <- c("wt","dn")
simple_df$background <- factor(simple_df$background,levels=background_order)

# cleanup
rm(counter,background,j,condition,i,group,area,decond,H3K27ac,spots,mean_rna)

####################################################################################################################
############### create "complex" dataframe for plotting ############################################################
####################################################################################################################
# create the background column dynamically
temp <- c()
background <- c()
counter <- 0
for (i in 1:length(conds)) {
  for (j in 1:length(reots)) {
    for (k in 1:length(celllines)) {
      temp <- c(rep(celllines[k],nrow(get(group_simple[(2*i-1)+counter]))))
      background <- c(background,temp)
      counter <- counter+1}
    counter <- 0}}

# create the condition column dynamically
temp <- c()
condition <- c()
for (i in 1:length(conds)) {
  for (j in 1:length(reots)) {
    temp <- rep(conds[i],nrow(get(group_simple[2*i-1]))+nrow(get(group_simple[2*i])))
    condition <- c(condition,temp)}}

# create the readout column dynamically
temp <- c()
readout <- c()
for (i in 1:length(conds)) {
  for (j in 1:length(reots)) {
    temp <- c(rep(reots[j],nrow(get(group_simple[2*i-1]))+nrow(get(group_simple[2*i]))))
    readout <- c(readout,temp)}}

# create the group column dynamically
temp <- c()
group <- c()
for (i in 1:length(group_complex)) {
  temp <- rep(group_complex[i],nrow(get(gsub("^(.*?)_(.*?)_.*","\\1_\\2",group_complex[i]))))
  group <- c(group,temp)}

# create readout columns dynamically
area <- unlist(lapply(group_complex,function(g) get(gsub("^(.*?)_(.*?)_.*","\\1_\\2",g))[,"area"]))
nuc_int <- unlist(lapply(group_complex,function(g) get(gsub("^(.*?)_(.*?)_.*","\\1_\\2",g))[,paste0("mean_",gsub("^.*_(.*)$","\\1",g))]))
enrich <- unlist(lapply(group_complex,function(g) {
  cc_col <- paste0("mean_",gsub("^.*_(.*)$","\\1",g),"_cc")
  np_col <- paste0("mean_",gsub("^.*_(.*)$","\\1",g),"_np")
  get(gsub("^(.*?)_(.*?)_.*","\\1_\\2",g))[,cc_col] / get(gsub("^(.*?)_(.*?)_.*","\\1_\\2",g))[,np_col]}))
spots <- unlist(lapply(group_complex,function(g) get(gsub("^(.*?)_(.*?)_.*","\\1_\\2",g))[,"spotcount"]))

# create "complex" dataframe
complex_df <- data.frame(background,condition,readout,group,area,nuc_int,enrich,spots)

# create a custom factor variable with the desired order for "background"
#background_order <- c("wt","dn")
#complex_df$background <- factor(complex_df$background,levels=background_order)

# cleanup
rm(cell_filter,background,counter,k,condition,readout,j,temp,group,i,area,nuc_int,enrich,spots)


####################################################################################################################
############################### Plots ##############################################################################
####################################################################################################################

############################### nucleus size #######################################################################
# violinplot for nucleus size readout
pViolin_nucsize <- ggplot(simple_df,aes(x=group,y=area)) +
  # Y/X axis titles
  xlab("") +
  ylab("Nuclear area (a.u.)") +
  # Violin plot with overlaid boxplot
  geom_violin(scale="width",lwd=0.5,aes(fill=group),trim=TRUE) +
  #geom_jitter() +
  # manual fill colors (alphabetically ordered group), use sort(x.labels_cc_area) to find out
  scale_fill_manual(values=c("#D20000","#D20000","#D20000","#D20000","#707070","#707070","#707070","#707070")) +
  geom_boxplot(width=0.4,outlier.size=0) +
  # keep the order of x elements as specified originally in the dataframe
  # usually, ggplot re-orders it alphabetically
  scale_x_discrete(limits=as.character(unique(simple_df$group)),labels=unique(simple_df[,"group"])) +
  # define tick marks of y axis
  scale_y_continuous(breaks=seq(0,50,by=5)) +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=2),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    #axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=1),
    axis.ticks.length=unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    legend.position = "none"
  )
ggsave(plot=pViolin_nucsize,width=unit(8,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_nucsize.pdf"))
rm(pViolin_nucsize)


############################### Decondensation #####################################################################
# violinplot for decondensation readout
pViolin_decond <- ggplot(simple_df,aes(x=condition,y=decond)) +
  # Y/X axis titles
  xlab("") +
  ylab("Nuclear area occupied by CC (%)") +
  # Violin plot with overlaid boxplot
  geom_violin(scale="width",lwd=0.5,aes(fill=background),alpha=0.75,trim=TRUE) +
  #geom_jitter() +
  # manual fill colors (alphabetically ordered background), use sort(x.labels_cc_area) to find out
  scale_fill_manual(values=c("#707070","#D20000")) +
  geom_boxplot(width=0.4,outlier.size=0) +
  # keep the order of x elements as specified originally in the dataframe
  # usually, ggplot re-orders it alphabetically
  scale_x_discrete(limits=as.character(unique(simple_df$condition)),labels=unique(simple_df[,"condition"])) +
  # define tick marks of y axis
  scale_y_continuous(breaks=seq(0,50,by=5)) +
  # make into two different plots split by cell background
  facet_wrap(~background,scales="free") +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=2),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    #axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=1),
    axis.ticks.length=unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_line(size=0.25,linetype='solid',colour="grey"), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    legend.position = "none"
  )
ggsave(plot=pViolin_decond,width=unit(7,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_decond.pdf"))
rm(pViolin_decond)


############################### Transcription (log) ################################################################
# violinplot for transcription readout
pViolin_trx <- ggplot(simple_df,aes(x=group,y=spots+1)) +
  # Y/X axis titles
  xlab("") +
  ylab("Number of transcripts + 1") +
  # Violin plot with overlaid boxplot
  geom_violin(scale="width",lwd=0.5,aes(fill=group),trim=TRUE) +
  #geom_jitter() +
  # manual fill colors (alphabetically ordered group), use sort(x.labels_cc_area) to find out
  scale_fill_manual(values=c("#D20000","#D20000","#D20000","#D20000","#707070","#707070","#707070","#707070")) +
  geom_boxplot(width=0.4,outlier.size=0) +
  # keep the order of x elements as specified originally in the dataframe
  # usually, ggplot re-orders it alphabetically
  scale_x_discrete(limits=as.character(unique(simple_df$group)),labels=unique(simple_df[,"group"])) +
  # define tick marks of y axis
  scale_y_continuous(breaks=c(seq(1,10,by=1),seq(10,100,by=10),seq(100,500,by=100)),trans = "log10") +
  # Set y-axis limits to exclude values above 100
  coord_cartesian(ylim = c(1, 500)) +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=2),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    #axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=1),
    axis.ticks.length=unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_line(size=0.25,linetype='solid',colour="grey"), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    legend.position = "none"
  )
ggsave(plot=pViolin_trx,width=unit(5,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_transcription.pdf"))
rm(pViolin_trx)


############################### Transcription (log) ################################################################
# violinplot for transcription readout
pViolin_trx <- ggplot(simple_df,aes(x=condition,y=spots+1)) +
  # Y/X axis titles
  xlab("") +
  ylab("Number of transcripts + 1") +
  # Violin plot with overlaid boxplot
  geom_violin(scale="width",lwd=0.5,aes(fill=background),alpha=0.75,trim=TRUE) +
  #geom_jitter() +
  # manual fill colors (alphabetically ordered background), use sort(x.labels_cc_area) to find out
  scale_fill_manual(values=c("#707070","#D20000")) +
  geom_boxplot(width=0.4,outlier.size=0) +
  # keep the order of x elements as specified originally in the dataframe
  # usually, ggplot re-orders it alphabetically
  scale_x_discrete(limits=as.character(unique(simple_df$condition)),labels=unique(simple_df[,"condition"])) +
  # define tick marks of y axis
  scale_y_continuous(breaks=c(seq(1,10,by=1),seq(10,100,by=10),seq(100,500,by=100)),trans = "log10") +
  # Set y-axis limits to exclude values above 100
  coord_cartesian(ylim = c(1, 500)) +
  # make into two different plots split by cell background
  facet_wrap(vars(simple_df$background)) +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=2),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    #axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=1),
    axis.ticks.length=unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_line(size=0.25,linetype='solid',colour="grey"), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    legend.position = "none"
  )
ggsave(plot=pViolin_trx,width=unit(7,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_transcription.pdf"))
rm(pViolin_trx)


############################### Transcription (cut) ################################################################
# violinplot for transcription readout
pViolin_trx <- ggplot(simple_df,aes(x=group,y=spots)) +
  # Y/X axis titles
  xlab("") +
  ylab("number of spots") +
  # Violin plot with overlaid boxplot
  geom_violin(scale="width",lwd=0.5,aes(fill=group),trim=TRUE) +
  #geom_jitter() +
  # manual fill colors (alphabetically ordered group), use sort(x.labels_cc_area) to find out
  scale_fill_manual(values=c("#D20000","#D20000","#D20000","#D20000","#707070","#707070","#707070","#707070")) +
  geom_boxplot(width=0.4,outlier.size=0) +
  # keep the order of x elements as specified originally in the dataframe
  # usually, ggplot re-orders it alphabetically
  scale_x_discrete(limits=as.character(unique(simple_df$group)),labels=unique(simple_df[,"group"])) +
  # define tick marks of y axis
  scale_y_continuous(breaks=seq(0,500,by=10)) +
  # Set y-axis limits to exclude values above 100
  coord_cartesian(ylim = c(0, 100)) +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=2),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    #axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=1),
    axis.ticks.length=unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_line(size=0.25,linetype='solid',colour="grey"), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    legend.position = "none"
  )
ggsave(plot=pViolin_trx,width=unit(5,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_transcription_cut.pdf"))
rm(pViolin_trx)


############################### Transcription (total nuclear signal) ###############################################
# violinplot for transcription readout
pViolin_trx <- ggplot(simple_df,aes(x=group,y=mean_rna*10)) +
  # Y/X axis titles
  xlab("") +
  ylab("mean intensity (a.u.)") +
  # Violin plot with overlaid boxplot
  geom_violin(scale="width",lwd=0.5,aes(fill=group),trim=TRUE) +
  #geom_jitter() +
  # manual fill colors (alphabetically ordered group), use sort(x.labels_cc_area) to find out
  scale_fill_manual(values=c("#D20000","#D20000","#D20000","#D20000","#707070","#707070","#707070","#707070")) +
  geom_boxplot(width=0.4,outlier.size=0) +
  # keep the order of x elements as specified originally in the dataframe
  # usually, ggplot re-orders it alphabetically
  scale_x_discrete(limits=as.character(unique(simple_df$group)),labels=unique(simple_df[,"group"])) +
  # define tick marks of y axis
  scale_y_continuous(breaks=seq(0,500,by=0.1)) +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=2),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    #axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=1),
    axis.ticks.length=unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_line(size=0.25,linetype='solid',colour="grey"), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    legend.position = "none"
  )
ggsave(plot=pViolin_trx,width=unit(5,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_transcription_meanrna.pdf"))
rm(pViolin_trx)


############################### spotcount vs decondensation ########################################################
# violinplots for spotcount vs decondensation
pViolin_spots_vs_decond <- ggplot(simple_df,aes(x=decond,y=spots)) +
  # Y/X axis titles
  xlab("Nuclear area occupied by CC (%)") +
  ylab("number of spots") +
  # Scatter plot
  geom_point(aes(color=condition),size=0.75) +
  scale_color_manual(values=c("GFP"="#000000","VP16"="#00CD00","p65"="#EEC900","VPR"="#D20000")) +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=2),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    #axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=1),
    axis.ticks.length = unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    #legend.position = "none"
  )
ggsave(plot=pViolin_spots_vs_decond,width=unit(8,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_spots_vs_decond.pdf"))
rm(pViolin_spots_vs_decond)


############################### spotcount vs decondensation ########################################################
# Calculate the average spot count for the control (GFP) condition
average_spot_count_control <- simple_df %>%
  filter(condition == "GFP") %>%
  summarise(average_spots = mean(spots)) %>%
  pull(average_spots)

# Print the result
print(sprintf("The average spot count for the control (GFP) condition is: %.2f", average_spot_count_control))



# Calculate the average spot count for the VP16 condition, considering only values up to 25% nuclear area occupied
average_spot_count_control <- simple_df %>%
  filter(condition == "VPR" & decond <= 25) %>%
  summarise(average_spots = median(spots)) %>%
  pull(average_spots)

# Print the result
print(sprintf("The average spot count for the VP16 condition (up to 25%% nuclear area occupied) is: %.2f", average_spot_count_control))

# Calculate the average spot count for the VP16 condition, considering only values up to 25% nuclear area occupied
average_spot_count_control <- simple_df %>%
  filter(condition == "VPR" & decond >= 25) %>%
  summarise(average_spots = median(spots)) %>%
  pull(average_spots)

# Print the result
print(sprintf("The average spot count for the VP16 condition (above 25%% nuclear area occupied) is: %.2f", average_spot_count_control))




# violinplots for spotcount vs decondensation
pViolin_spots_vs_decond <- ggplot(simple_df,aes(x=decond,y=spots)) +
  # Y/X axis titles
  xlab("Nuclear area occupied by CC (%)") +
  ylab("number of spots") +
  # Scatter plot
  geom_point(aes(color=condition),size=0.75) +
  scale_color_manual(values=c("GFP"="#000000","VP16"="#3ECB76","p65"="#FF961D","VPR"="#D10000")) +
  # Set axis limits
  #scale_x_continuous(limits = c(0.75, 2)) +
  scale_y_continuous(limits = c(0, 100)) +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=2),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    #axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=1),
    axis.ticks.length = unit(3,"mm"),
    # Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    #legend.position = "none"
  )
ggsave(plot=pViolin_spots_vs_decond,width=unit(8,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_spots_vs_decond_limited.pdf"))
rm(pViolin_spots_vs_decond)




library(ggplot2)
library(dplyr)
library(patchwork)

# Calculate the average spot count for the control (GFP) condition
average_spot_count_control <- simple_df %>%
  filter(condition == "GFP") %>%
  summarise(average_spots = mean(spots)) %>%
  pull(average_spots)

# Print the result
print(sprintf("The average spot count for the control (GFP) condition is: %.2f", average_spot_count_control))

# Function to create plot for each condition
create_plot <- function(data, condition) {
  ggplot(data, aes(x=decond, y=spots)) +
    xlab("Nuclear area occupied by CC (%)") +
    ylab("number of spots") +
    geom_point(color=switch(condition,
                            "GFP" = "#000000",
                            "VP16" = "#3ECB76",
                            "p65" = "#FF961D",
                            "VPR" = "#D10000"),
               size=0.75) +
    scale_y_continuous(limits = c(0, 100)) +
    ggtitle(condition) +
    theme(
      panel.border = element_rect(color="black",fill=NA,size=2),
      plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
      axis.text = element_text(size=10,colour="black"),
      axis.ticks.y = element_line(colour="black",size=1),
      axis.ticks.length = unit(3,"mm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

# Create a list of plots
plot_list <- lapply(unique(simple_df$condition), function(cond) {
  create_plot(simple_df %>% filter(condition == cond), cond)
})

# Combine plots
combined_plot <- wrap_plots(plot_list, ncol = 2)

# Save combined plot
ggsave(plot=combined_plot, width=unit(16,"cm"), height=unit(12,"cm"), dpi=300,
       filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_spots_vs_decond_separate_conditions.pdf"))

# Cleanup
rm(plot_list, combined_plot, create_plot)




############################### spotcount vs decondensation (per cell type) ########################################
# violinplots for spotcount vs decondensation
pViolin_spots_vs_decond_ct <- ggplot(simple_df,aes(x=decond,y=spots)) +
  # Y/X axis titles
  xlab("Nuclear area occupied by CC (%)") +
  ylab("number of spots") +
  # Scatter plot
  geom_point(aes(color=background),size=0.75) +
  scale_color_manual(values=c("wt"="#000000","dn"="#D20000")) +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=2),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    #axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=1),
    axis.ticks.length = unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    #legend.position = "none"
  )
ggsave(plot=pViolin_spots_vs_decond_ct,width=unit(8,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_spots_vs_decond_ct.pdf"))
rm(pViolin_spots_vs_decond_ct)



############################### H3K27ac enrich vs spotcount ########################################################
# violinplots for H3K27ac enrichment vs spotcount
pViolin_spots_vs_H3K27acenrich <- ggplot(simple_df,aes(x=H3K27ac,y=spots)) +
  # Y/X axis titles
  xlab("H3K27ac enrichment at CC (a.u.)") +
  ylab("number of spots") +
  # Scatter plot
  geom_point(aes(color=condition),size=0.75) +
  scale_color_manual(values=c("GFP"="#000000","VP16"="#00CD00","p65"="#EEC900","VPR"="#D20000")) +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=2),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    #axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=1),
    axis.ticks.length = unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    #legend.position = "none"
  )
ggsave(plot=pViolin_spots_vs_H3K27acenrich,width=unit(8,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_spots_vs_H3K27acenrich.pdf"))
rm(pViolin_spots_vs_H3K27acenrich)


############################### H3K27ac enrich vs spotcount zoom-in ################################################
# violinplots for H3K27ac enrichment vs spotcount
pViolin_spots_vs_H3K27acenrich <- ggplot(simple_df, aes(x=H3K27ac, y=spots+1)) +
  # Y/X axis titles
  xlab("H3K27ac enrichment at CC (a.u.)") +
  ylab("number of transcripts + 1") +
  # Scatter plot
  geom_point(aes(color=condition), size=0.75) +
  scale_color_manual(values=c("GFP"="#000000","VP16"="#3ECB76","p65"="#FF961D","VPR"="#D10000")) +
  # Set axis limits
  scale_x_continuous(limits = c(0.75, 2)) +
  scale_y_continuous(limits = c(0, 200)) +
  # Theme/style elements
  theme(
    panel.border = element_rect(color="black",fill=NA,size=2),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    axis.text = element_text(size=10,colour="black"),
    axis.ticks.y = element_line(colour="black",size=1),
    axis.ticks.length = unit(3,"mm"),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )

ggsave(plot=pViolin_spots_vs_H3K27acenrich, width=unit(8,"cm"), height=unit(6,"cm"), dpi=300,
       filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_spots_vs_H3K27acenrich_limited.pdf"))
rm(pViolin_spots_vs_H3K27acenrich)


############################### H3K27ac enrich vs spotcount zoom-in  + residuals and 99% confidence interval #######
library(ggplot2)
library(gridExtra)

# Fit a linear model
lm_model <- lm(spots ~ H3K27ac, data = simple_df)

# Calculate R-squared
r_squared <- summary(lm_model)$r.squared

# Create the main plot with linear fit and 95% confidence interval
pViolin_spots_vs_H3K27acenrich <- ggplot(simple_df, aes(x=H3K27ac, y=spots+1)) +
  xlab("H3K27ac enrichment at CC (a.u.)") +
  ylab("number of transcripts + 1") +
  geom_point(aes(color=condition), size=0.75) +
  scale_color_manual(values=c("GFP"="#000000","VP16"="#3ECB76","p65"="#FF961D","VPR"="#D10000")) +
  geom_smooth(method = "lm", se = TRUE, level = 0.99, color = "blue", fill = "lightblue", alpha = 0.3) + # Add linear fit with 95% CI
  scale_x_continuous(limits = c(0.75, 2)) +
  scale_y_continuous(limits = c(0, 200)) +
  annotate("text", x = 1.9, y = 190, label = sprintf("R² = %.3f", r_squared), hjust = 1) +
  theme(
    panel.border = element_rect(color="black",fill=NA,size=2),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    axis.text = element_text(size=10,colour="black"),
    axis.ticks.y = element_line(colour="black",size=1),
    axis.ticks.length = unit(3,"mm"),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )

# Create residuals plot
residuals_plot <- ggplot(data.frame(fitted = fitted(lm_model), residuals = residuals(lm_model)), 
                         aes(x = fitted, y = residuals)) +
  geom_point(size = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("Fitted values") +
  ylab("Residuals") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1),
    axis.text = element_text(size=8, colour="black"),
    axis.ticks = element_line(colour="black", size=0.5),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

# Combine plots
combined_plot <- grid.arrange(pViolin_spots_vs_H3K27acenrich, residuals_plot, ncol = 2)

# Save combined plot
ggsave(plot = combined_plot, width = unit(12,"cm"), height = unit(6,"cm"), dpi = 300,
       filename = paste0(folder_curated, timestamp, "_plotCC_PDFs", "/pViolin_spots_vs_H3K27acenrich_linear_with_residuals.pdf"))

# Cleanup
rm(pViolin_spots_vs_H3K27acenrich, residuals_plot, combined_plot, lm_model, r_squared)


############################### H3K27ac enrich vs spotcount (per cell type) ########################################
# violinplots for H3K27ac enrichment vs spotcount
pViolin_spots_vs_H3K27acenrich_ct <- ggplot(simple_df,aes(x=H3K27ac,y=spots)) +
  # Y/X axis titles
  xlab("H3K27ac enrichment at CC (a.u.)") +
  ylab("number of spots") +
  # Scatter plot
  geom_point(aes(color=background),size=0.75) +
  scale_color_manual(values=c("wt"="#000000","dn"="#D20000")) +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=2),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    #axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=1),
    axis.ticks.length = unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    #legend.position = "none"
  )
ggsave(plot=pViolin_spots_vs_H3K27acenrich_ct,width=unit(8,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_spots_vs_H3K27acenrich_ct.pdf"))
rm(pViolin_spots_vs_H3K27acenrich_ct)


############################### H3K27ac enrich vs decondensation ###################################################
# violinplots for H3K27ac enrichment vs decondensation
pViolin_H3K27acenrich_vs_decond <- ggplot(simple_df,aes(x=decond,y=H3K27ac)) +
  # Y/X axis titles
  xlab("Nuclear area occupied by CC (%)") +
  ylab("H3K27ac enrichment at CC (a.u.)") +
  # Scatter plot
  geom_point(aes(color=condition),size=0.75) +
  scale_color_manual(values=c("GFP"="#000000", "VP16"="#00CD00", "p65"="#EEC900", "VPR"="#D20000")) +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=2),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    #axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=1),
    axis.ticks.length = unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    #legend.position = "none"
  )
ggsave(plot=pViolin_H3K27acenrich_vs_decond,width=unit(8,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_H3K27acenrich_vs_decond.pdf"))
rm(pViolin_H3K27acenrich_vs_decond)


############################### H3K27ac enrich vs decondensation GFP only + fit #############################################
# Calculate linear model and R-squared
lm_model <- lm(H3K27ac ~ decond, data = simple_df[simple_df$condition == "GFP", ])
r_squared <- summary(lm_model)$r.squared

# Create the plot with linear fit
pViolin_H3K27acenrich_vs_decond <- ggplot(simple_df %>% filter(condition == "GFP"), aes(x=decond, y=H3K27ac)) +
  geom_point(color = "#000000", size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 0.5) +
  xlab("Nuclear area occupied by CCs (%)") +
  ylab("H3K27ac enrichment at CCs (a.u.)") +
  geom_text(aes(x = max(decond), y = max(H3K27ac), 
                label = sprintf("R² = %.3f", r_squared)),
            hjust = 1, vjust = 1, size = 3) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    axis.text = element_text(size=8, colour="black"),
    axis.ticks.y = element_line(colour="black", size=0.5),
    axis.ticks.x = element_line(colour="black", size=0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none"
  )
ggsave(plot=pViolin_H3K27acenrich_vs_decond, width=unit(4,"cm"), height=unit(3,"cm"), dpi=300,
       filename=paste0(folder_curated, timestamp, "_plotCC_PDFs/pViolin_H3K27acenrich_vs_decond_GFP_lm.pdf"))

# Cleanup
rm(lm_model,r_squared,pViolin_H3K27acenrich_vs_decond)


############################### H3K27ac enrich vs decondensation (per cell type) ###################################
# violinplots for H3K27ac enrichment vs decondensation
pViolin_H3K27acenrich_vs_decond_ct <- ggplot(simple_df,aes(x=decond,y=H3K27ac)) +
  # Y/X axis titles
  xlab("Nuclear area occupied by CC (%)") +
  ylab("H3K27ac enrichment at CC (a.u.)") +
  # Scatter plot
  geom_point(aes(color=background),size=0.75) +
  scale_color_manual(values=c("wt"="#000000","dn"="#D20000")) +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=2),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    #axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=1),
    axis.ticks.length = unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    #legend.position = "none"
  )
ggsave(plot=pViolin_H3K27acenrich_vs_decond_ct,width=unit(8,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_H3K27acenrich_vs_decond_ct.pdf"))
rm(pViolin_H3K27acenrich_vs_decond_ct)


############################### Mean nuclear intensity #############################################################
# violinplots for average absolute nuclear intensity values
pViolin_nuc_int <- ggplot(complex_df,aes(x=group,y=nuc_int)) +
  # Y/X axis titles
  xlab("") +
  ylab("Mean nuclear intensity (a.u.)") +
  # Violin plot with overlaid boxplot
  geom_violin(scale="width",lwd=0.4,aes(fill=group)) +
  # manual fill colors (alphabetically ordered group), use sort(unique(complex_df$group)) to find out
  scale_fill_manual(values=c("#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000",
                             "#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000",
                             "#D20000","#D20000","#D20000","#D20000","#707070","#707070","#707070","#707070",
                             "#707070","#707070","#707070","#707070","#707070","#707070","#707070","#707070",
                             "#707070","#707070","#707070","#707070","#707070","#707070","#707070","#707070")) +
  geom_boxplot(width=0.4,outlier.size=0,linewidth=0.3) +
  # keep the order of x elements as specified originally in the dataframe
  # usually, ggplot re-orders it alphabetically
  scale_x_discrete(limits=as.character(unique(complex_df$group)),labels=unique(complex_df[,"group"])) +
  # define tick marks of y axis
  scale_y_continuous(breaks=seq(0,0.6,by=0.1)) +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=1),
    plot.margin = unit(c(1,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=0.5),
    axis.ticks.length = unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_line(size=0.25,linetype='solid',colour="grey"), 
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    legend.position = "none"
  )
ggsave(plot=pViolin_nuc_int,width=unit(9,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_nuc_int.pdf"))
rm(pViolin_nuc_int)


############################### Mean nuclear intensity (no RNA) ####################################################
# violinplots for average absolute nuclear intensity values
pViolin_nuc_int <- ggplot(complex_df[complex_df$readout!="rna",],aes(x=group,y=nuc_int)) +
  # Y/X axis titles
  xlab("") +
  ylab("Mean nuclear intensity (a.u.)") +
  # Violin plot with overlaid boxplot
  geom_violin(scale="width",lwd=0.4,aes(fill=group)) +
  # manual fill colors (alphabetically ordered group), use sort(unique(complex_df$group)) to find out
  scale_fill_manual(values=c("#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000",
                             "#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000",
                             "#707070","#707070","#707070","#707070","#707070","#707070","#707070","#707070",
                             "#707070","#707070","#707070","#707070","#707070","#707070","#707070","#707070")) +
  geom_boxplot(width=0.4,outlier.size=0,linewidth=0.3) +
  # keep the order of x elements as specified originally in the dataframe
  # usually, ggplot re-orders it alphabetically
  scale_x_discrete(limits=as.character(unique(complex_df[complex_df$readout!="rna",]$group)),labels=unique(complex_df[complex_df$readout!="rna","group"])) +
  # define tick marks of y axis
  scale_y_continuous(breaks=seq(0,0.6,by=0.1)) +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=1),
    plot.margin = unit(c(1,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=0.5),
    axis.ticks.length = unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_line(size=0.25,linetype='solid',colour="grey"), 
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    legend.position = "none"
  )
ggsave(plot=pViolin_nuc_int,width=unit(9,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_nuc_int_norna.pdf"))
rm(pViolin_nuc_int)


############################### Enrichment (cc/np) #################################################################
# violinplots for average absolute nuclear intensity values
pViolin_enrich <- ggplot(complex_df,aes(x=group,y=enrich)) +
  # Y/X axis titles
  xlab("") +
  ylab("Normalized intensity (cc/np) (a.u.)") +
  # Violin plot with overlaid boxplot
  geom_violin(scale="width",lwd=0.5,aes(fill=group)) +
  # manual fill colors (alphabetically ordered group), use sort(unique(complex_df$group)) to find out
  scale_fill_manual(values=c("#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000",
                             "#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000",
                             "#D20000","#D20000","#D20000","#D20000","#707070","#707070","#707070","#707070",
                             "#707070","#707070","#707070","#707070","#707070","#707070","#707070","#707070",
                             "#707070","#707070","#707070","#707070","#707070","#707070","#707070","#707070")) +
  geom_boxplot(width=0.4,outlier.size=0,linewidth=0.3) +
  # keep the order of x elements as specified originally in the dataframe
  # usually, ggplot re-orders it alphabetically
  scale_x_discrete(limits=as.character(unique(complex_df$group)), labels=unique(complex_df[,"group"])) +
  # define tick marks of y axis
  scale_y_continuous(breaks=seq(0,6,by=0.5)) +
  # facet into different plots for better visualization
  #facet_wrap("condition") +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=1),
    plot.margin = unit(c(1,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=0.5),
    axis.ticks.length = unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_line(size=0.25,linetype='solid',colour="grey"), 
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    legend.position = "none"
  )
ggsave(plot=pViolin_enrich,width=unit(8,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_enrich.pdf"))
rm(pViolin_enrich)


############################### Enrichment (cc/np) (no RNA) ########################################################
# violinplots for average absolute nuclear intensity values
pViolin_enrich <- ggplot(complex_df[complex_df$readout!="rna",],aes(x=group,y=enrich)) +
  # Y/X axis titles
  xlab("") +
  ylab("Normalized intensity (cc/np) (a.u.)") +
  # Violin plot with overlaid boxplot
  geom_violin(scale="width",lwd=0.5,aes(fill=group)) +
  # manual fill colors (alphabetically ordered group), use sort(unique(complex_df$group)) to find out
  scale_fill_manual(values=c("#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000",
                             "#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000","#D20000",
                             "#707070","#707070","#707070","#707070","#707070","#707070","#707070","#707070",
                             "#707070","#707070","#707070","#707070","#707070","#707070","#707070","#707070")) +
  geom_boxplot(width=0.4,outlier.size=0,linewidth=0.3) +
  # keep the order of x elements as specified originally in the dataframe
  # usually, ggplot re-orders it alphabetically
  scale_x_discrete(limits=as.character(unique(complex_df[complex_df$readout!="rna",]$group)),labels=unique(complex_df[complex_df$readout!="rna","group"])) +
  # define tick marks of y axis
  scale_y_continuous(breaks=seq(0,6,by=0.5)) +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=1),
    plot.margin = unit(c(1,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=0.5),
    axis.ticks.length = unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_line(size=0.25,linetype='solid',colour="grey"), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    legend.position = "none"
  )
ggsave(plot=pViolin_enrich,width=unit(8,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_enrich_norna.pdf"))
rm(pViolin_enrich)


############################### Enrichment (cc/np) (no RNA) only wt ################################################
# filter complex dataframe to only contain wt conditions
complex_df_wt <- subset(complex_df, !grepl("dn", background))

# violinplot for enrichment values
pViolin_enrich <- ggplot(complex_df_wt[complex_df_wt$readout!="rna",], aes(x=group, y=enrich)) +
  # Y/X axis titles
  xlab("") +
  ylab("Normalized intensity (cc/np) (a.u.)") +
  # Violin plot with overlaid boxplot
  geom_violin(scale="width", lwd=0.5, aes(fill=group)) +
  # manual fill colors (adjust based on the number of unique groups in complex_df_wt)
  scale_fill_manual(values=c(rep("#707070",8),rep("#f1f1f1",8))) +
  geom_boxplot(width=0.4, outlier.size=0, linewidth=0.3) +
  # keep the order of x elements as specified originally in the dataframe
  scale_x_discrete(limits=as.character(unique(complex_df_wt[complex_df_wt$readout!="rna",]$group)),
                   labels=unique(complex_df_wt[complex_df_wt$readout!="rna","group"])) +
  # define tick marks of y axis
  scale_y_continuous(breaks=seq(0,6,by=0.5)) +
  # Theme/style elements (unchanged)
  theme(
    panel.border = element_rect(color="black",fill=NA,size=1),
    plot.margin = unit(c(1,0.5,0.5,0.5),"cm"),
    axis.text = element_text(size=10,colour="black"),
    axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=0.5),
    axis.ticks.length = unit(3,"mm"),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_line(size=0.25,linetype='solid',colour="grey"), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none"
  )
ggsave(plot=pViolin_enrich,width=unit(8,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_enrich_norna_wt.pdf"))
rm(complex_df_wt,pViolin_enrich)


############################### Enrichment (cc/np) (RNA only) ######################################################
# Filter rows based on the "rna" phrase in the "group" column
filtered_complex_df <- complex_df[grep("rna",complex_df$group,ignore.case=T),]

# violinplots for average absolute nuclear intensity values
pViolin_enrich_rnao <- ggplot(complex_df[grep("rna",complex_df$group,ignore.case=T),],aes(x=group,y=enrich)) +
  # Y/X axis titles
  xlab("") +
  ylab("Normalized intensity (cc/np) (a.u.)") +
  # Violin plot with overlaid boxplot
  geom_violin(scale="width",lwd=0.5,aes(fill=group)) +
  # manual fill colors (alphabetically ordered group), use sort(unique(filtered_complex_df$group)) to find out
  scale_fill_manual(values=c("#D20000","#D20000","#D20000","#D20000","#707070","#707070","#707070","#707070")) +
  geom_boxplot(width=0.4,outlier.size=0,linewidth=0.3) +
  # keep the order of x elements as specified originally in the dataframe
  # usually, ggplot re-orders it alphabetically
  scale_x_discrete(limits=as.character(unique(filtered_complex_df$group)), labels=unique(filtered_complex_df[,"group"])) +
  # define tick marks of y axis
  scale_y_continuous(breaks=seq(0,4,by=0.5)) +
  # Set y-axis limits to exclude values above 3.5
  coord_cartesian(ylim = c(1, 3.5)) +
  # facet into different plots for better visualization
  #facet_wrap("condition") +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=1),
    plot.margin = unit(c(1,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=0.5),
    axis.ticks.length = unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_line(size=0.25,linetype='solid',colour="grey"), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    legend.position = "none"
  )
ggsave(plot=pViolin_enrich_rnao,width=unit(4,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_enrich_rnao.pdf"))
rm(filtered_complex_df,pViolin_enrich_rnao)


############################### Enrichment (cc/np) (RNA only, cut) #################################################
# Filter rows based on the "rna" phrase in the "group" column
filtered_complex_df <- complex_df[grep("rna",complex_df$group,ignore.case=T),]

# violinplots for average absolute nuclear intensity values
pViolin_enrich_rnao <- ggplot(complex_df[grep("rna",complex_df$group,ignore.case=T),],aes(x=group,y=enrich)) +
  # Y/X axis titles
  xlab("") +
  ylab("Normalized intensity (cc/np) (a.u.)") +
  # Violin plot with overlaid boxplot
  geom_violin(scale="width",lwd=0.5,aes(fill=group)) +
  # manual fill colors (alphabetically ordered group), use sort(unique(filtered_complex_df$group)) to find out
  scale_fill_manual(values=c("#D20000","#D20000","#D20000","#D20000","#707070","#707070","#707070","#707070")) +
  geom_boxplot(width=0.4,outlier.size=0,linewidth=0.3) +
  # keep the order of x elements as specified originally in the dataframe
  # usually, ggplot re-orders it alphabetically
  scale_x_discrete(limits=as.character(unique(filtered_complex_df$group)), labels=unique(filtered_complex_df[,"group"])) +
  # define tick marks of y axis
  scale_y_continuous(breaks=seq(0,4,by=0.05)) +
  # Set y-axis limits to exclude values above 1.5
  coord_cartesian(ylim = c(1, 1.5)) +
  # facet into different plots for better visualization
  #facet_wrap("condition") +
  # Theme/style elements
  theme(
    # give the plot a border and margin
    panel.border = element_rect(color="black",fill=NA,size=1),
    plot.margin = unit(c(1,0.5,0.5,0.5),"cm"),
    # define axis text color, fontsize and tick properties
    axis.text = element_text(size=10,colour="black"),
    axis.text.x = element_text(angle=90),
    axis.ticks.y = element_line(colour="black",size=0.5),
    axis.ticks.length = unit(3,"mm"),
    # remove x axis ticks
    axis.ticks.x = element_blank(),
    # Remove grid lines
    panel.grid.major = element_line(size=0.25,linetype='solid',colour="grey"), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    # Change panel background
    panel.background = element_blank(),
    #remove legend
    legend.position = "none"
  )
ggsave(plot=pViolin_enrich_rnao,width=unit(4,"cm"),height=unit(6,"cm"),dpi=300,filename=paste0(folder_curated,timestamp,"_plotCC_PDFs","/pViolin_enrich_rnao_cut.pdf"))
rm(filtered_complex_df,pViolin_enrich_rnao)

