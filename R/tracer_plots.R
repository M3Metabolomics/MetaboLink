
# Normalization functions

normalizeRowSums <- function(data) {
    normalized <- data / rowSums(data, na.rm = TRUE)
    return(normalized)
}

normalizeReference <- function(data) {
    normalized <- data / data[,1]
    return(normalized)
}



# FRACTIONAL CONTRIBUTION

selectFCtable <- function(data, sequence, plot_settings) {
    fc_data <- data[, c("Analyte", "Analysis", "FC")]

    # Filter data for selected metabolite and group
    fc_data <- fc_data[fc_data$Analyte == plot_settings$metabolite_group & fc_data$Analysis %in% plot_settings$sample, ]

    fc_data <- merge(fc_data, sequence[, c('sample', 'time')], by.x = 'Analysis', by.y = 'sample')

    return(fc_data)
}

plotFractionalContribution <- function(data, sequence, plot_settings) {
    # Filter data for selected metabolite and group
    group_data <- data[data$Analyte == plot_settings$metabolite_group & data$Analysis %in% plot_settings$sample, c('Analyte', 'Analysis', 'FC')]
    
    # Merge with sequence time column 
    group_data <- merge(group_data, sequence[, c('sample', 'time')], by.x = 'Analysis', by.y = 'sample')

    # Remove missing rows with NA in FC or time columns
    group_data <- group_data[!is.na(group_data$FC) & !is.na(group_data$time), ]
    
    # Summarize FC data
    melted_data <- group_data %>%
        group_by(time) %>%
        summarise(mean_abundance = mean(FC, na.rm = TRUE),
                  se_abundance = sd(FC, na.rm = TRUE) / sqrt(n()),
                  .groups = 'drop')

    # Plot
    if (plot_settings$plot_type == "errorbar") {
        plot <- ggplot(melted_data, aes(x = as.factor(time), y = mean_abundance)) +
            geom_col(fill = "steelblue", alpha = 0.7) +
            geom_errorbar(aes(ymin = mean_abundance - se_abundance, ymax = mean_abundance + se_abundance), width = 0.2) +
            labs(x = "Time", y = "Fractional Contribution ± SE", title = paste("FC ± SE —", plot_settings$metabolite_group, "(Group:", plot_settings$group, ")")) +
            theme_minimal()
    } else {
        plot <- ggplot(melted_data, aes(x = as.factor(time), y = mean_abundance)) +
            geom_col(fill = "steelblue", alpha = 0.7) +
            labs(x = "Time", y = "Fractional Contribution", title = paste("FC —", plot_settings$metabolite_group, "(Group:", plot_settings$group, ")")) +
            theme_minimal()
    }

    return(plot)
}


# TRACER PLOTS (1)
# Plot metabolite & sample "Isotopologue distribution — [metabolite] (Normalized to total%)"
plotIsotopologueDist <- function(data, plot_settings) {
    metabolite_data <- data[data$Analyte == plot_settings$metabolite & data$Analysis == plot_settings$sample, ]

    melted_data <- reshape2::melt(metabolite_data, id.vars = c("Analyte", "Analysis"), variable.name = "isotopologue", value.name = "abundance")

    plot <- ggplot(melted_data, aes(x = isotopologue, y = abundance, fill = isotopologue)) +
        geom_bar(stat = "identity") +
        labs(x = "Isotopologues", y = "Abundance", title = paste("Isotopologue distribution —", plot_settings$metabolite, "(Normalized to total%)")) +
        theme_minimal() +
        theme(legend.position = "none")

    return(plot)
}


# Plot metabolite & sample (A.0 normalized): "Isotopologue distribution — [metabolite] (Normalized to A.0)"
plotIsotopologueDistA0 <- function(data, plot_settings) {
    metabolite_data <- data[data$Analyte == plot_settings$metabolite & data$Analysis == plot_settings$sample, ]
    melted_data <- reshape2::melt(metabolite_data, id.vars = c("Analyte", "Analysis"), variable.name = "isotopologue", value.name = "abundance")

    plot <- ggplot(melted_data, aes(x = isotopologue, y = abundance, fill = isotopologue)) +
        geom_bar(stat = "identity") +
        labs(x = "Isotopologues", y = "Abundance", title = paste("Isotopologue distribution —", plot_settings$metabolite, "(Normalized to A.0)")) +
        theme_minimal() +
        theme(legend.position = "none")

    return(plot)
}


# TRACER PLOTS (2)
# Plot one metabolite "Isotopologue profile over time — [metabolite] (Group: [group])" and errorbar variant "Mean isotopologue abundance ± SE — [metabolite] (Group: [group])"
plotIsotopologueProfile <- function(data, sequence, plot_settings) {
    # Filter data for selected metabolite and group
    group_data <- sequence[sequence[,'group'] %in% plot_settings$group, 'sample']
    filtered_data <- data[data$Analyte == plot_settings$metabolite & data$Analysis == group_data, ]
    
    # Melt data to long format
    melted_data <- reshape2::melt(filtered_data, id.vars = c("Analyte", "Analysis", "time"), variable.name = "isotopologue", value.name = "abundance")
    
    # Plot
    plot <- ggplot(melted_data, aes(x = time, y = abundance, color = isotopologue)) +
        geom_line() +
        geom_point() +
        labs(x = "Time", y = "Abundance", title = paste("Isotopologue profile over time —", plot_settings$metabolite, "(Group:", plot_settings$group_time, ")")) +
        theme_minimal()

    return(plot)
}

# Plot isotopologues stacked (per group_time): "Stacked isotopologue abundances — [metabolite] ([Raw/Normalized])" and errorbar variant "Stacked isotopologue abundances (Mean ± SE) — [metabolite] ([Raw/Normalized])"
plotStackedIsotopologues <- function(data, plot_settings) {
    
    # Melt data to long format
    melted_data <- reshape2::melt(data, id.vars = c("Analyte", "Analysis"), variable.name = "isotopologue", value.name = "abundance")

    # ensure isotopologue is ordered so A+0 is at the bottom of stacked bars
    iso_levels <- unique(melted_data$isotopologue)
    # try to extract numeric part after '+' to sort numerically; fallback to lexical order
    nums <- suppressWarnings(as.numeric(gsub(".*\\+(\\d+).*", "\\1", iso_levels)))
    if (all(!is.na(nums))) {
        iso_levels <- iso_levels[order(nums, decreasing = TRUE)]
    } else {
        iso_levels <- sort(iso_levels, decreasing = TRUE)
    }
    melted_data$isotopologue <- factor(melted_data$isotopologue, levels = iso_levels)

    # Plot
    plot <- ggplot(melted_data, aes(x = Analysis, y = abundance, fill = isotopologue)) +
        geom_col(position = "stack") +
        labs(x = plot_settings$group_time, y = "Abundance", title = paste("Stacked isotopologue abundances —", plot_settings$metabolite_iso)) +
        theme_minimal()

    return(plot)
}


# Plot multiple groups, one time: "Isotopologue distribution across groups — [metabolite] at time [time]" and errorbar variant add "(Mean ± SE)"
plotIsotopologueDistAcrossGroups <- function(data, plot_settings) {
    # Filter data for selected metabolite and time point
    filtered_data <- data[data$Analyte == plot_settings$metabolite_time_table & data$time == plot_settings$time_point, ]
    
    # Melt data to long format
    melted_data <- reshape2::melt(filtered_data, id.vars = c("Analyte", "Analysis", "group"), variable.name = "isotopologue", value.name = "abundance")
    
    plot <- ggplot(melted_data, aes(x = group, y = abundance, fill = isotopologue)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(x = "Group", y = "Abundance", title = paste("Isotopologue distribution across groups —", plot_settings$metabolite_time_table, "at time", plot_settings$time_point)) +
        theme_minimal()

    return(plot)
}
