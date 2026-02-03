############### Description ##########################################################################################


############### User-defined parameters ##############################################################################
folder_curated <- "insert file path/_analysis_VPR/"
cell_filter <- "g"  # "g" for good cells, "b" for bad cells, "gb" for both
conditions <- c("ctrl", "VPR")
channels <- c("DNA", "dCas9", "H3K9me3", "HP1")

############### Setup and data reading ###############################################################################
# Load required libraries
library(ggplot2)
library(ggsignif)

# Create timestamp variable for this run
date_time <- format(Sys.time(), "analysis_%Y%m%d_%H%M%S")

# Create output folder
output_folder <- file.path(folder_curated, format(Sys.time(), "analysis_%Y%m%d_%H%M%S"))
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Read data
read_data <- function(channel) {
  file <- list.files(folder_curated, pattern = paste0(channel, "_mean_sd_projection_curated"), full.names = TRUE)
  read.table(file, header = TRUE, stringsAsFactors = FALSE)
}
data_conditions <- lapply(conditions, read_data)

# Filter data
filter_data <- function(data) {
  if (cell_filter == "gb") {
    data[data$annotation %in% c("g", "b"), ]
  } else {
    data[data$annotation == cell_filter, ]
  }
}
data_filtered <- lapply(data_conditions, filter_data)

# Check if filtered data is empty
if (any(sapply(data_filtered, nrow) == 0)) {
  stop("No cells with specified quality annotation were found. Check filter settings and dataset annotations.")
}

# Cleanup
rm(cell_filter,read_data,data_conditions,filter_data)


########## Mean nuclear size #####################################################################
# Define function for data preparation
prepare_area_data <- function(data_list, conditions) {
  plot_data <- data.frame()
  for (i in seq_along(data_list)) {
    condition_data <- data.frame(
      group = paste0("area_", conditions[i]),
      value = data_list[[i]]$area,
      status = conditions[i]
    )
    plot_data <- rbind(plot_data, condition_data)
  }
  plot_data
}
# Use the function
area_data <- prepare_area_data(data_filtered, conditions)

# Define function for violin plot
create_area_violinplot <- function(data, y_label, y_breaks, filename) {
  x_labels <- c("area_ctrl", "area_VPR")
  p <- ggplot(data, aes(x = group, y = value)) +
    geom_violin(scale = "width", lwd = 0.5, aes(fill = group)) +
    geom_boxplot(width = 0.4, outlier.size = 0) +
    scale_fill_manual(values = c("#707070", "#f1f1f1")) +
    scale_x_discrete(limits = as.character(unique(data$group)), labels = x_labels) +
    scale_y_continuous(breaks = y_breaks) +
    labs(x = "", y = y_label) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 2),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      axis.text = element_text(size = 16, colour = "black"),
      axis.ticks.y = element_line(colour = "black", size = 1),
      axis.ticks.length = unit(3, "mm"),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  ggsave(plot = p, width = 8, height = 6, dpi = 300, filename = filename)
}
# Generate plot using the function
create_area_violinplot(area_data, "Nuclear area (pixels)", seq(0, max(area_data$value), by = 5000), file.path(output_folder, "pViolin_area_ViolinBoxplot.pdf"))

# Calculate mean, median, sd, and SEM for each dataframe in the list
area_stats_list <- lapply(seq_along(data_filtered), function(i) {
  df <- data_filtered[[i]]
  data.frame(
    condition = conditions[i],
    mean = mean(df$area),
    median = median(df$area),
    sd = sd(df$area),
    sem = sd(df$area) / sqrt(length(df$area))
  )
})
# Combine the results into a single dataframe
area_stats <- do.call(rbind, area_stats_list)

# Cleanup
rm(prepare_area_data,create_area_violinplot,area_stats_list)


########## Mean nuclear intensity ################################################################
# Define function for data preparation
prepare_nucint_data <- function(data_list, conditions, channels) {
  plot_data <- data.frame()
  for (channel in channels) {
    channel_data <- data.frame(
      group = rep(paste0(channel, "_", conditions), sapply(data_list, nrow)),
      value = unlist(lapply(data_list, function(x) x[[paste0("mean_", channel)]])),
      status = rep(conditions, sapply(data_list, nrow))
    )
    plot_data <- rbind(plot_data, channel_data)
  }
  plot_data
}
# Use the function
nucint_data <- prepare_nucint_data(data_filtered, conditions, channels)

# Define function for violin plot
create_nucint_violinplot <- function(data, y_label, y_breaks, filename) {
  x_labels <- c("dna-", "dnavpr", "dcas9-", "dcas9vpr", "H3K9me3-", "H3K9me3vpr", "HP1-", "HP1vpr")
  p <- ggplot(data, aes(x = group, y = value)) +
    geom_violin(scale = "width", lwd = 0.5, aes(fill = group)) +
    geom_boxplot(width = 0.4, outlier.size = 0) +
    scale_fill_manual(values = rep(c("#707070", "#f1f1f1"), 4)) +
    scale_x_discrete(limits = as.character(unique(data$group)), labels = x_labels) +
    scale_y_continuous(breaks = y_breaks) +
    labs(x = "", y = y_label) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 2),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      axis.text = element_text(size = 16, colour = "black"),
      axis.ticks.y = element_line(colour = "black", size = 1),
      axis.ticks.length = unit(3, "mm"),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  ggsave(plot = p, width = 8, height = 6, dpi = 300, filename = filename)
}
# Generate plot using the function
create_nucint_violinplot(nucint_data, "Mean nuclear intensity (a.u.)", seq(0, 0.6, by = 0.1), file.path(output_folder, "pViolin_nucint_ViolinBoxplot.pdf"))

# Define function to calculate statistics for each condition and channel
calculate_nucint_stats <- function(data_list, conditions, channels) {
  stats_list <- lapply(seq_along(data_list), function(i) {
    df <- data_list[[i]]
    condition <- conditions[i]
    
    channel_stats <- lapply(channels, function(channel) {
      values <- df[[paste0("mean_", channel)]]
      data.frame(
        condition = condition,
        channel = channel,
        median = median(values),
        mean = mean(values),
        sd = sd(values),
        sem = sd(values) / sqrt(length(values))
      )
    })
    do.call(rbind, channel_stats)
  })
  # Combine all statistics into a single dataframe
  do.call(rbind, stats_list)
}

# Calculate statistics for each condition and channel
nucint_stats <- calculate_nucint_stats(data_filtered, conditions, channels)

# Cleanup
rm(prepare_nucint_data,create_nucint_violinplot,calculate_nucint_stats)


########## Mean integrated nuclear intensity #####################################################
# Define function for data preparation
prepare_nucintegint_data <- function(data_list, conditions, channels) {
  plot_data <- data.frame()
  for (channel in channels) {
    channel_data <- data.frame(
      group = rep(paste0(channel, "_", conditions), sapply(data_list, nrow)),
      value = unlist(lapply(data_list, function(x) x[[paste0("mean_", channel)]] * x$area)),
      status = rep(conditions, sapply(data_list, nrow))
    )
    plot_data <- rbind(plot_data, channel_data)
  }
  plot_data
}
# Use the function
nucintegint_data <- prepare_nucintegint_data(data_filtered, conditions, channels)

# Define function for violin plot
create_nucintegint_violinplot <- function(data, y_label, y_breaks, filename) {
  x_labels <- c("dna-", "dnavpr", "dcas9-", "dcas9vpr", "H3K9me3-", "H3K9me3vpr", "HP1-", "HP1vpr")
  p <- ggplot(data, aes(x = group, y = value)) +
    geom_violin(scale = "width", lwd = 0.5, aes(fill = group)) +
    geom_boxplot(width = 0.4, outlier.size = 0) +
    scale_fill_manual(values = rep(c("#707070", "#f1f1f1"), 4)) +
    scale_x_discrete(limits = as.character(unique(data$group)), labels = x_labels) +
    scale_y_continuous(breaks = y_breaks, limits = c(0, max(y_breaks))) +
    labs(x = "", y = y_label) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 2),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      axis.text = element_text(size = 16, colour = "black"),
      axis.ticks.y = element_line(colour = "black", size = 1),
      axis.ticks.length = unit(3, "mm"),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  ggsave(plot = p, width = 8, height = 6, dpi = 300, filename = filename)
}
# Generate plot using the function
create_nucintegint_violinplot(nucintegint_data, "Mean integrated nuclear intensity (a.u.)", seq(0, 6000, by = 1000), file.path(output_folder, "pViolin_nucintegint_ViolinBoxplot.pdf"))

# Define function to calculate statistics for each condition and channel
calculate_nucintegint_stats <- function(data_list, conditions, channels) {
  stats_list <- lapply(seq_along(data_list), function(i) {
    df <- data_list[[i]]
    condition <- conditions[i]
    channel_stats <- lapply(channels, function(channel) {
      values <- df[[paste0("mean_", channel)]] * df$area
      data.frame(
        condition = condition,
        channel = channel,
        median = median(values),
        mean = mean(values),
        sd = sd(values),
        sem = sd(values) / sqrt(length(values))
      )
    })
    do.call(rbind, channel_stats)
  })
  # Combine all statistics into a single dataframe
  do.call(rbind, stats_list)
}

# Calculate statistics for each condition and channel
nucintegint_stats <- calculate_nucintegint_stats(data_filtered, conditions, channels)

# Cleanup
rm(prepare_nucintegint_data,create_nucintegint_violinplot,calculate_nucintegint_stats)


########## Mean nucleoplasm intensity / Mean chromocenter intensity ##############################
# Define function for data preparation
prepare_maskint_data <- function(data_list, conditions, channels) {
  plot_data <- data.frame()
  for (channel in channels) {
    for (compartment in c("np", "cc")) {
      for (condition in conditions) {
        data <- data_list[[which(conditions == condition)]]
        new_data <- data.frame(
          group = paste0(channel, "_", compartment, "_", condition),
          value = data[[paste0("mean_", channel, "_", compartment)]],
          status = condition,
          compartment = compartment
        )
        plot_data <- rbind(plot_data, new_data)
      }
    }
  }
  # Define the desired order
  desired_order <- c()
  for (channel in channels) {
    for (compartment in c("np", "cc")) {
      for (condition in conditions) {
        desired_order <- c(desired_order, paste0(channel, "_", compartment, "_", condition))
      }
    }
  }
  # Create a factor to ensure the desired order
  plot_data$group <- factor(plot_data$group, levels = desired_order)
  return(plot_data)
}
# Use the function
maskint_data <- prepare_maskint_data(data_filtered, conditions, channels)

# Define function for violin plot
create_maskint_violinplot <- function(data, y_label, y_breaks, filename) {
  x_labels <- levels(data$group)
  x_labels <- gsub("_", "-", x_labels)  # Replace underscores with hyphens in labels
  p <- ggplot(data, aes(x = group, y = value)) +
    geom_violin(scale = "width", lwd = 0.5, aes(fill = group)) +
    geom_boxplot(width = 0.4, outlier.size = 0) +
    scale_fill_manual(values = rep(c("#707070", "#f1f1f1"), length(unique(data$group))/2)) +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(breaks = y_breaks) +
    labs(x = "", y = y_label) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 2),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      axis.text.y = element_text(size = 16, colour = "black"),
      axis.text.x = element_text(size = 12, colour = "black", angle = 90, vjust = 0.5, hjust = 1),
      axis.ticks.y = element_line(colour = "black", size = 1),
      axis.ticks.length = unit(3, "mm"),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  ggsave(plot = p, width = 12, height = 6, dpi = 300, filename = filename)
}
# Generate plot using the function
create_maskint_violinplot(maskint_data, "Mean mask intensities (a.u.)", seq(0, 1, by = 0.1), file.path(output_folder, "pViolin_maskint_ViolinBoxplot.pdf"))

# Define function to calculate statistics for each condition, channel, and compartment
calculate_maskint_stats <- function(data_list, conditions, channels) {
  stats_list <- lapply(seq_along(data_list), function(i) {
    df <- data_list[[i]]
    condition <- conditions[i]
    
    channel_stats <- lapply(channels, function(channel) {
      compartment_stats <- lapply(c("np", "cc"), function(compartment) {
        values <- df[[paste0("mean_", channel, "_", compartment)]]
        data.frame(
          condition = condition,
          channel = paste0(channel, "_", compartment),
          median = median(values),
          mean = mean(values),
          sd = sd(values),
          sem = sd(values) / sqrt(length(values))
        )
      })
      do.call(rbind, compartment_stats)
    })
    do.call(rbind, channel_stats)
  })
  # Combine all statistics into a single dataframe
  do.call(rbind, stats_list)
}
# Calculate statistics for each condition, channel, and compartment
maskint_stats <- calculate_maskint_stats(data_filtered, conditions, channels)

# Cleanup
rm(prepare_maskint_data,create_maskint_violinplot,calculate_maskint_stats)


########## Mean integrated nucleoplasm intensity / Mean integrated chromocenter intensity ########
# Define function for data preparation
prepare_maskintegint_data <- function(data_list, conditions, channels) {
  plot_data <- data.frame()
  for (channel in channels) {
    for (compartment in c("np", "cc")) {
      for (condition in conditions) {
        data <- data_list[[which(conditions == condition)]]
        new_data <- data.frame(
          group = paste0(channel, "_", compartment, "_", condition),
          value = data[[paste0("mean_", channel, "_", compartment)]] * data[[paste0("area_", compartment)]],
          status = condition,
          compartment = compartment
        )
        plot_data <- rbind(plot_data, new_data)
      }
    }
  }
  # Define the desired order
  desired_order <- c()
  for (channel in channels) {
    for (compartment in c("np", "cc")) {
      for (condition in conditions) {
        desired_order <- c(desired_order, paste0(channel, "_", compartment, "_", condition))
      }
    }
  }
  # Create a factor to ensure the desired order
  plot_data$group <- factor(plot_data$group, levels = desired_order)
  return(plot_data)
}
# Use the function
maskintegint_data <- prepare_maskintegint_data(data_filtered, conditions, channels)

# Define function for violin plot
create_maskintegint_violinplot <- function(data, y_label, y_breaks, filename) {
  x_labels <- levels(data$group)
  x_labels <- gsub("_", "-", x_labels)  # Replace underscores with hyphens in labels
  p <- ggplot(data, aes(x = group, y = value)) +
    geom_violin(scale = "width", lwd = 0.5, aes(fill = group)) +
    geom_boxplot(width = 0.4, outlier.size = 0) +
    scale_fill_manual(values = rep(c("#707070", "#f1f1f1"), length(unique(data$group))/2)) +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(breaks = y_breaks) +
    labs(x = "", y = y_label) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 2),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      axis.text.y = element_text(size = 16, colour = "black"),
      axis.text.x = element_text(size = 12, colour = "black", angle = 90, vjust = 0.5, hjust = 1),
      axis.ticks.y = element_line(colour = "black", size = 1),
      axis.ticks.length = unit(3, "mm"),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  ggsave(plot = p, width = 12, height = 6, dpi = 300, filename = filename)
}
# Generate plot using the function
create_maskintegint_violinplot(maskintegint_data, "Mean integrated mask intensities (a.u.)", seq(0, 6000, by = 1000), file.path(output_folder, "pViolin_maskintegint_ViolinBoxplot.pdf"))

# Define function to calculate statistics for each condition, channel, and compartment
calculate_maskintegint_stats <- function(data_list, conditions, channels) {
  stats_list <- lapply(seq_along(data_list), function(i) {
    df <- data_list[[i]]
    condition <- conditions[i]
    channel_stats <- lapply(channels, function(channel) {
      compartment_stats <- lapply(c("np", "cc"), function(compartment) {
        values <- df[[paste0("mean_", channel, "_", compartment)]] * df[[paste0("area_", compartment)]]
        data.frame(
          condition = condition,
          channel = channel,
          compartment = compartment,
          median = median(values),
          mean = mean(values),
          sd = sd(values),
          sem = sd(values) / sqrt(length(values)),
          n = length(values)
        )
      })
      do.call(rbind, compartment_stats)
    })
    do.call(rbind, channel_stats)
  })
  # Combine all statistics into a single dataframe
  do.call(rbind, stats_list)
}
# Calculate statistics for each condition, channel, and compartment
maskintegint_stats <- calculate_maskintegint_stats(data_filtered, conditions, channels)

# Cleanup
rm(prepare_maskintegint_data,create_maskintegint_violinplot,calculate_maskintegint_stats)


########## Mean enrichment #######################################################################
# Define function for data preparation
prepare_enrich_data <- function(data_list, conditions, channels) {
  plot_data <- data.frame()
  for (channel in channels) {
    channel_data <- data.frame(
      group = rep(paste0(channel, "_", conditions), sapply(data_list, nrow)),
      value = unlist(lapply(data_list, function(x) x[[paste0("mean_", channel, "_cc")]] / x[[paste0("mean_", channel, "_np")]])),
      status = rep(conditions, sapply(data_list, nrow))
    )
    plot_data <- rbind(plot_data, channel_data)
  }
  plot_data
}
# Use the function
enrich_data <- prepare_enrich_data(data_filtered, conditions, channels)

# Define function for enrichment violin plot
create_enrich_violinplot <- function(data, y_label, y_breaks, filename) {
  x_labels <- c("dna-", "dnavpr", "dcas9-", "dcas9vpr", "H3K9me3-", "H3K9me3vpr", "HP1-", "HP1vpr")
  p <- ggplot(data, aes(x = group, y = value)) +
    geom_violin(scale = "width", lwd = 0.5, aes(fill = group)) +
    geom_boxplot(width = 0.4, outlier.size = 0) +
    scale_fill_manual(values = rep(c("#707070", "#f1f1f1"), 4)) +
    scale_x_discrete(limits = as.character(unique(data$group)), labels = x_labels) +
    scale_y_continuous(breaks = y_breaks) +
    labs(x = "", y = y_label) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 2),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      axis.text = element_text(size = 16, colour = "black"),
      axis.ticks.y = element_line(colour = "black", size = 1),
      axis.ticks.length = unit(3, "mm"),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  ggsave(plot = p, width = 8, height = 6, dpi = 300, filename = filename)
}
# Generate plot using the function
create_enrich_violinplot(enrich_data, "Normalized intensity (cc/np)", seq(0, 6, by = 0.5), file.path(output_folder, "pViolin_enrich_ViolinBoxplot.pdf"))

# Define function to calculate enrichment statistics for each condition and channel
calculate_enrich_stats <- function(data_list, conditions, channels) {
  stats_list <- lapply(seq_along(data_list), function(i) {
    df <- data_list[[i]]
    condition <- conditions[i]
    channel_stats <- lapply(channels, function(channel) {
      cc_values <- df[[paste0("mean_", channel, "_cc")]]
      np_values <- df[[paste0("mean_", channel, "_np")]]
      enrichment_values <- cc_values / np_values
      data.frame(
        condition = condition,
        channel = channel,
        median = median(enrichment_values),
        mean = mean(enrichment_values),
        sd = sd(enrichment_values),
        sem = sd(enrichment_values) / sqrt(length(enrichment_values))
      )
    })
    do.call(rbind, channel_stats)
  })
  # Combine all statistics into a single dataframe
  do.call(rbind, stats_list)
}
# Calculate statistics for each condition and channel
enrich_stats <- calculate_enrich_stats(data_filtered, conditions, channels)

# Calculate relative enrichment (VPR / ctrl) and its standard error
calculate_relative_enrichment <- function(enrich_stats) {
  channels <- unique(enrich_stats$channel)
  relative_enrichment <- sapply(channels, function(channel) {
    ctrl_stats <- enrich_stats[enrich_stats$condition == "ctrl" & enrich_stats$channel == channel, ]
    vpr_stats <- enrich_stats[enrich_stats$condition == "VPR" & enrich_stats$channel == channel, ]
    rel_enrich <- vpr_stats$mean / ctrl_stats$mean
    rel_enrich_sem <- rel_enrich * sqrt(
      (vpr_stats$sem / vpr_stats$mean)^2 + 
        (ctrl_stats$sem / ctrl_stats$mean)^2
    )
    c(relative_enrichment = rel_enrich, sem = rel_enrich_sem)
  })
  
  t(relative_enrichment)
}
# Calculate relative enrichment
relative_enrichment <- calculate_relative_enrichment(enrich_stats)

# Create data frame for plotting
relative_enrichment_df <- data.frame(
  channel = rownames(relative_enrichment),
  relative_enrichment = relative_enrichment[, "relative_enrichment"],
  sem = relative_enrichment[, "sem"]
)

# Create bar plot function
create_bar_plot <- function(data, y_label, filename) {
  # Reorder the channels
  data$channel <- factor(data$channel, levels = c("DNA", "dCas9", "H3K9me3", "HP1"))
  
  p <- ggplot(data, aes(x = channel, y = relative_enrichment)) +
    geom_bar(stat = "identity", width = 0.6, fill = "grey50") +
    geom_errorbar(aes(ymin = relative_enrichment - sem, ymax = relative_enrichment + sem), 
                  width = 0.2, color = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    labs(x = "", y = y_label) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      axis.text = element_text(size = 12, colour = "black"),
      axis.title = element_text(size = 14),
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.ticks.length = unit(2, "mm"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  ggsave(plot = p, width = 6, height = 5, dpi = 300, filename = filename)
}

# Generate relative enrichment bar plot
create_bar_plot(relative_enrichment_df, "Relative enrichment (VPR / ctrl)", file.path(output_folder, "pBar_relative_enrichment.pdf"))

# Cleanup
rm(prepare_enrich_data,create_enrich_violinplot,calculate_enrich_stats,calculate_relative_enrichment,create_bar_plot)

