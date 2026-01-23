
# Normalization functions

normalizeRowSums <- function(data) {
    normalized <- data / rowSums(data, na.rm = TRUE)
    return(normalized)
}

normalizeReference <- function(data) {
    normalized <- data / data[,1]
    return(normalized)
}

format4ggplot <- function(data, sequence, isotopologues) {
    colnames(sequence)[which(names(sequence) == "sample")] <- "Analysis"
    merged_df <- left_join(data, sequence, by = "Analysis")

    pivoted <- merged_df %>%
      pivot_longer(
        cols = all_of(isotopologues),
        names_to = "Isotopologue",
        values_to = "Abundance"
      ) %>%
      select(Analyte, Isotopologue, Abundance, Analysis, Time = time, Group = group, groupTime = group_time) %>%
      filter(!is.na(Abundance))
  
    # Order isotopologues by the numeric part after '+' (A+0, A+1, A+2, A+10 -> correct numeric order)
    iso_levels <- unique(pivoted$Isotopologue)
    nums <- suppressWarnings(as.numeric(gsub(".*\\+(\\d+).*", "\\1", iso_levels)))
    if (all(!is.na(nums))) {
      iso_levels <- iso_levels[order(nums)]
    } else {
      iso_levels <- sort(iso_levels)
    }
    pivoted$Isotopologue <- factor(pivoted$Isotopologue, levels = iso_levels)

    # Ensure time is numeric and order groupTime by Group then Time
    pivoted$Time <- as.numeric(as.character(pivoted$Time))
    group_levels <- pivoted %>%
      dplyr::distinct(Group, Time, groupTime) %>%
      arrange(Group, Time) %>%
      pull(groupTime)

    group_levels <- unique(group_levels)
    pivoted$groupTime <- factor(pivoted$groupTime, levels = group_levels)
    pivoted$Group <- factor(pivoted$Group, levels = unique(pivoted$Group))
  

    return(pivoted)
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

    plot <- ggplot(metabolite_data, aes(x = Isotopologue, y = Abundance, fill = Isotopologue)) +
        geom_bar(stat = "identity") +
        labs(x = "Isotopologues", y = "Abundance", title = paste("Isotopologue distribution —", plot_settings$metabolite, "(Normalized to total%)")) +
        theme_minimal() +
        theme(legend.position = "none")

    return(plot)
}


# TRACER PLOTS (2)
# Plot one metabolite "Isotopologue profile over time — [metabolite] (Group: [group])" and errorbar variant "Mean isotopologue abundance ± SE — [metabolite] (Group: [group])"
plotIsotopologueProfile <- function(data, sequence, plot_settings) {
    # Filter data for selected metabolite and group
    filtered_data <- data[data$Analyte == plot_settings$metabolite & data$Group %in% plot_settings$group, ]
     
    plot <- ggplot(filtered_data, aes(x = Time, y = Abundance, color = Isotopologue)) +
        geom_line() +
        geom_point() +
        labs(x = "Time", y = "Abundance", title = paste("Isotopologue profile over time —", plot_settings$metabolite, "(Group:", plot_settings$group_time, ")")) +
        theme_minimal()

    return(plot)
}

# Plot isotopologues stacked (per group_time): "Stacked isotopologue abundances — [metabolite] ([Raw/Normalized])" and errorbar variant "Stacked isotopologue abundances (Mean ± SE) — [metabolite] ([Raw/Normalized])"
plotStackedIsotopologues <- function(data, plot_settings) {
    
    melted_data <- data
    # ensure isotopologue is ordered so A+0 is at the bottom of stacked bars
    iso_levels <- unique(melted_data$Isotopologue)
    # try to extract numeric part after '+' to sort numerically; fallback to lexical order
    nums <- suppressWarnings(as.numeric(gsub(".*\\+(\\d+).*", "\\1", iso_levels)))
    if (all(!is.na(nums))) {
        iso_levels <- iso_levels[order(nums, decreasing = TRUE)]
    } else {
        iso_levels <- sort(iso_levels, decreasing = TRUE)
    }
    melted_data$Isotopologue <- factor(melted_data$Isotopologue, levels = iso_levels)
    

    if (plot_settings$plot_type == "errorbar") {
        plot <- ggplot(melted_data, aes(x = groupTime, y = mean_abundance, fill = Isotopologue)) +
          geom_col(position = "stack") +
          geom_errorbar(aes(ymin = error_lower, ymax = error_upper), width = 0.2, position = "identity") +
          labs(x = "Group / Time", y = "Abundance", title = paste("Stacked isotopologue abundances (Mean ± SE) —", plot_settings$metabolite_iso)) +
          theme_minimal()
    }
    else {
       plot <- ggplot(melted_data, aes(x = groupTime, y = mean_abundance, fill = Isotopologue)) +
        geom_col(position = "stack") +
        labs(x = "Group / Time", y = "Abundance", title = paste("Stacked isotopologue abundances —", plot_settings$metabolite_iso)) +
        theme_minimal()
    }

    plot <- ggplotly(plot)

    return(plot)
}


# Plot group over time: "Isotopologue distribution over time — [metabolite] in group [group]" and errorbar variant add "(Mean ± SE)"
plotIsotopologueDistAcrossGroups <- function(data, plot_settings) {
    # Filter data for selected metabolite and group
    filtered_data <- data[data$Analyte == plot_settings$metabolite & data$Group == plot_settings$group, ]

    # Summarize --> mean and se for errorbar plot
    filtered_data <- filtered_data %>%
        group_by(Time, Isotopologue) %>%
        summarise(mean_abundance = mean(Abundance, na.rm = TRUE),
                  se_abundance = sd(Abundance, na.rm = TRUE) / sqrt(n()),
                  error_lower = mean_abundance - se_abundance,
                  error_upper = mean_abundance + se_abundance,
                  .groups = 'drop')
    

    plot <- ggplot(filtered_data, aes(x = Time, y = mean_abundance, fill = Isotopologue)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(x = "Group", y = "Abundance", title = paste("Isotopologue distribution over time —", plot_settings$metabolite, ", group", plot_settings$group)) +
        theme_minimal()
    plot <- ggplotly(plot)
    return(plot)
}




# FROM MFA3 APP
selectGTtable <- function(data, sequence, settings) {
  filtered <- data %>%
    filter(Analyte == settings$metabolite & Group %in% settings$group) %>%
    select(Isotopologue, Analysis, Abundance, Time) %>%
    group_by(Isotopologue, Time) %>%
    dplyr::summarise(
        mean_x = mean(Abundance, na.rm = TRUE),
        se = sd(Abundance, na.rm = TRUE) / sqrt(n()),
        .groups = 'drop'
    ) %>%
    ungroup()

    filtered <- filtered[filtered$mean_x > 0, ]

    iso_levels <- unique(filtered$Isotopologue)
    # try to extract numeric part after '+' to sort numerically; fallback to lexical order
    nums <- suppressWarnings(as.numeric(gsub(".*\\+(\\d+).*", "\\1", iso_levels)))
    if (all(!is.na(nums))) {
        iso_levels <- iso_levels[order(nums, decreasing = TRUE)]
    } else {
        iso_levels <- sort(iso_levels, decreasing = TRUE)
    }
    filtered$Isotopologue <- factor(filtered$Isotopologue, levels = iso_levels)

  return(filtered)
}

plotGroup <- function(data, sequence, settings) {
    
    filtered <- selectGTtable(data, sequence, settings)
  #filtered <- summarySE(filtered, measurevar = "X1", groupvars = c("Isotopologue", "Time"), na.rm=TRUE)

  plot <- ggplotly(ggplot(filtered, aes(fill = Isotopologue, x = as.factor(Time), y = mean_x)) +
    geom_bar(position = "stack", stat = "identity") +
    # geom_errorbar(aes(ymin=X1-se, ymax=X1+se),
    #           linewidth=.3,    # Thinner lines
    #           width=.2,
    #           position=position_dodge(.9)) +
    labs(x = "Time", y = "Abundance", title = "Group distribution") +
    theme_bw())

    return(plot)
}


# Plot 4 - relative abundance of isotopologues for mulitple metabolites, one sample
# ALWAYS NORMALIZED DATA
plotMultipleAnalytes <- function(data, sample, metabolites) {
  #TODO when selecting different sized (?) metabolites, it doesn't look good

    filtered <- data %>%
        filter(Analysis == sample & Analyte %in% metabolites) %>%
        select(Isotopologue, X1, Analyte) %>%
        mutate(Analyte = factor(Analyte, levels = metabolites))

    plot <- ggplot(filtered, aes(fill = Isotopologue, x = Analyte, y = X1)) +
        geom_bar(position = "stack", stat = "identity") +
        labs(x = "Metabolites", y = "Abundance", title = sample) +
        theme_bw()

    plot <- ggplotly(plot)
    return(plot)
}

#Metabolite per group
plotIsotopologue <- function(data, sequence, plot_settings) {
  if (plot_settings$plot_type == "barplot") {
    # Bar plot logic
    metabolite_data <- data[data$Analyte == plot_settings$metabolite, ]
    
    plot <- ggplot(metabolite_data, aes(x = as.factor(Time), y = Abundance, fill = Isotopologue)) +
      geom_bar(position = "dodge", stat = "identity") +
      labs(
        x = "Time", 
        y = "Relative Abundance", 
        title = paste("Isotopologue distribution across time -", plot_settings$metabolite),
        fill = "Isotopologue"
      ) +
      theme_minimal()
    
  } else if (plot_settings$plot_type == "errorbar") {
    # Error bar plot logic
    filtered_data <- data[data$Analyte == plot_settings$metabolite, ]
    
    # Calculate summary statistics
    summary_data <- filtered_data %>%
      group_by(Isotopologue, Time) %>%
      summarise(
        mean_abundance = mean(Abundance, na.rm = TRUE),
        se_abundance = sd(Abundance, na.rm = TRUE) / sqrt(n()),
        .groups = 'drop'
      )
    
    plot <- ggplot(summary_data, aes(x = as.factor(Time), y = mean_abundance, 
                                     fill = Isotopologue)) +
      geom_bar(position = position_dodge(width = 0.9), stat = "identity", width = 0.7) +
      geom_errorbar(
        aes(ymin = mean_abundance - se_abundance, ymax = mean_abundance + se_abundance),
        position = position_dodge(width = 0.9),
        width = 0.25
      ) +
      labs(
        x = "Time", 
        y = "Relative Abundance", 
        title = paste("Mean isotopologue abundance ± SE -", plot_settings$metabolite),
        fill = "Isotopologue"
      ) +
      theme_minimal()
  }
  
  return(ggplotly(plot))
}

#Group X Time Plot
# Time Plot Function
plotGroupTime <- function(data, sequence, plot_settings) {
  if (plot_settings$plot_type == "errorbar") {
    # Error bar plot logic
    summary_data <- data %>%
      filter(Analyte == plot_settings$metabolite, Time %in% plot_settings$time_points) %>%
      group_by(Group, Isotopologue, Time) %>%
      summarise(
        mean_abundance = mean(Abundance, na.rm = TRUE),
        se_abundance = sd(Abundance, na.rm = TRUE) / sqrt(n()),
        .groups = 'drop'
      )
    
    plot <- ggplot(summary_data, aes(x = Group, y = mean_abundance, 
                                     fill = Isotopologue)) +
      geom_bar(position = position_dodge(width = 0.9), stat = "identity", width = 0.7) +
      geom_errorbar(
        aes(ymin = mean_abundance - se_abundance, ymax = mean_abundance + se_abundance),
        position = position_dodge(width = 0.9),
        width = 0.25
      ) +
      labs(
        x = "Group", 
        y = "Relative Abundance", 
        title = paste(plot_settings$metabolite, "- Time:", paste(plot_settings$time_points, collapse = ", "), "(Error Bars)"),
        fill = "Isotopologue"
      ) +
      theme_minimal()
    
  } else {
    # Regular bar plot logic
    filtered_data <- data %>%
      filter(Analyte == plot_settings$metabolite, Time %in% plot_settings$time_points)
    
    plot <- ggplot(filtered_data, aes(x = Group, y = Abundance, fill = Isotopologue)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(
        x = "Group", 
        y = "Relative Abundance", 
        title = paste(plot_settings$metabolite, "- Time:", paste(plot_settings$time_points, collapse = ", ")),
        fill = "Isotopologue"
      ) +
      theme_minimal()
  }
  
  return(ggplotly(plot))
}

#Normalized to A.0 plot
plotIsotopologueDist2 <- function(normalizeReference, plot_settings) {
  metabolite_data <- normalizeReference[normalizeReference$Analyte == plot_settings$metabolite & normalizeReference$Analysis == plot_settings$sample, ]

    print(metabolite_data)
  print(nrow(metabolite_data))
  
  plot <- ggplot(metabolite_data, aes(x = Isotopologue, y = Abundance, fill = Isotopologue)) +
    geom_bar(stat = "identity") +
    labs(x = "Isotopologues", y = "Abundance", title = paste("Isotopologue distribution —", plot_settings$metabolite, "(Normalized to A.0)")) +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(ggplotly(plot))
}