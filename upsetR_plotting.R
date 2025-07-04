library(ComplexUpset)
library(tidyr)
library(dplyr)

reach_scores <- read.csv("/data/UCSF_ComeBACK/REACH_scores_20250423_clean.csv")
reach_scores <- reach_scores[!is.na(reach_scores$X50_50_probs), ]
reach_scores <- reach_scores[reach_scores$visit == "baseline", ]
reach_scores$participant_id <- gsub("CMBK-010", "", reach_scores$participant_id)
reach_scores$participant_id <- gsub("-", "", reach_scores$participant_id)

file <- "tmcs/outputs/subtype_stage_results_PSMatched.csv" # use this file if doing the lkm for subtype since it includes all the participants

df <- read.csv(file)

df$stage <- df$stage * 21

df <- merge(df, reach_scores, by = "participant_id")

df <- df[, -which(names(df) == "C_DMGENDER")]
df <- df[, -which(names(df) == "C_ANBMISRCL")]
df <- df[, -which(names(df) == "smoking")]
df <- df[, -which(names(df) == "C_DMAGE")]

#########################################################################################################
load("/tmcs/data/REACH_data.rda")
#########################################################################################################

threshold <- 4
# binarise the variables of interest if variables containing ddd are >= threshold make them one else make them zero
df$disc_ddd_12_binary <- ifelse(df$disc_ddd_12_nl >= threshold, 1, 0)
df$disc_ddd_23_binary <- ifelse(df$disc_ddd_23_nl >= threshold, 1, 0)
df$disc_ddd_34_binary <- ifelse(df$disc_ddd_34_nl >= threshold, 1, 0)
df$disc_ddd_45_binary <- ifelse(df$disc_ddd_45_nl >= threshold, 1, 0)
df$disc_ddd_51_binary <- ifelse(df$disc_ddd_51_nl >= threshold, 1, 0)

threshold <- 2
# repeat the above but for facet_fj_oa_12_nl, etc
df$facet_fj_oa_12_binary <- ifelse(df$facet_fj_oa_12_nl >= threshold, 1, 0)
df$facet_fj_oa_23_binary <- ifelse(df$facet_fj_oa_23_nl >= threshold, 1, 0)
df$facet_fj_oa_34_binary <- ifelse(df$facet_fj_oa_34_nl >= threshold, 1, 0)
df$facet_fj_oa_45_binary <- ifelse(df$facet_fj_oa_45_nl >= threshold, 1, 0)
df$facet_fj_oa_51_binary <- ifelse(df$facet_fj_oa_51_nl >= threshold, 1, 0)

# change the column names to be more descriptive
colnames(df)[colnames(df) == "disc_ddd_12_binary"] <- "L1-L2 IVD"
colnames(df)[colnames(df) == "disc_ddd_23_binary"] <- "L2-L3 IVD"
colnames(df)[colnames(df) == "disc_ddd_34_binary"] <- "L3-L4 IVD"      
colnames(df)[colnames(df) == "disc_ddd_45_binary"] <- "L4-L5 IVD"
colnames(df)[colnames(df) == "disc_ddd_51_binary"] <- "L5-S1 IVD"

colnames(df)[colnames(df) == "facet_fj_oa_12_binary"] <- "L1-L2 FJ"
colnames(df)[colnames(df) == "facet_fj_oa_23_binary"] <- "L2-L3 FJ"
colnames(df)[colnames(df) == "facet_fj_oa_34_binary"] <- "L3-L4 FJ"
colnames(df)[colnames(df) == "facet_fj_oa_45_binary"] <- "L4-L5 FJ"
colnames(df)[colnames(df) == "facet_fj_oa_51_binary"] <- "L5-S1 FJ"


df$subtype <- as.character(df$subtype)

subset_data <- df

subset_data$subtype <- ifelse(subset_data$subtype == 1, "Facet-first", subset_data$subtype)
subset_data$subtype <- ifelse(subset_data$subtype == 0, "Disc-first", subset_data$subtype)

# intersect_columns <- rev(c("L1-L2 IVD", "L2-L3 IVD", "L3-L4 IVD", "L4-L5 IVD", "L5-S1 IVD"))
# plot_name <- "IVD degeneration intersections"

intersect_columns <- rev(c("L1-L2 FJ", "L2-L3 FJ", "L3-L4 FJ", "L4-L5 FJ", "L5-S1 FJ"))
plot_name <- "Facet joint OA intersections"


subset_data$intersection_id <- interaction(subset_data[, intersect_columns], drop = TRUE)
median_stage_by_intersection <- subset_data %>%
  group_by(intersection_id) %>%
  summarise(median_stage = median(stage, na.rm = TRUE)) %>%
  arrange(median_stage)  # Or `asc()` if you want ascending order
subset_data$intersection_id <- factor(
  subset_data$intersection_id,
  levels = median_stage_by_intersection$intersection_id
)
sorted_intersections <- median_stage_by_intersection$intersection_id
levels_vec <- levels(sorted_intersections)  # Assuming this is your factor

convert_to_labels <- function(code_str, labels) {
  bits <- strsplit(code_str, "\\.")[[1]]         # Split by '.'
  active <- labels[as.logical(as.numeric(bits))] # Keep only labels where bit == 1
  if (length(active) == 0) "None" else paste(active, collapse = " & ")
}
pretty_levels <- sapply(levels_vec, convert_to_labels, labels = intersect_columns)
names(pretty_levels) <- levels_vec
pretty_levels[pretty_levels == "None"] <- "Outside of known sets"
pretty_sorted <- unname(pretty_levels[sorted_intersections])
sorted_intersections_list <- strsplit(pretty_sorted, " & ")
sorted_intersections_list <- lapply(sorted_intersections_list, trimws)

selected_columns <- c(
  "L1-L2 IVD", "L2-L3 IVD", "L3-L4 IVD", "L4-L5 IVD", "L5-S1 IVD", 
  "L1-L2 FJ", "L2-L3 FJ", "L3-L4 FJ", "L4-L5 FJ", "L5-S1 FJ",
  "Age", "BMI", "stage", "subtype"
)
rating_scale = scale_fill_manual(values=c(
    'Disc-first'='#377EB8',
    'Facet-first'='#FF7F00'
))
show_hide_scale = scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide='none')
text_size <- 15
# Subset the data frame
subset_data <- subset_data[selected_columns]
size = get_size_mode("exclusive_intersection")
common_ymax <- 100
upset(
  subset_data,
  intersect = intersect_columns,
  sort_sets = FALSE,
  # min_size = 5         # Only show intersections with at least 3 items
  sort_intersections = FALSE,
  intersections = sorted_intersections_list,
  n_intersections = 10,
  ,name = plot_name
  ,width_ratio = 0.2,     # Adjust layout  
  keep_empty_groups = TRUE,
  set_sizes=(
    upset_set_size(
        geom=geom_bar(
            aes(fill=subtype, x=group),
            width=0.8
        ),
        position='right'
    )
  ),
  guides='over',
  annotations = list(
    "Age" = upset_annotate("Age", geom_boxplot(na.rm = TRUE)),
    "EBM stage" = upset_annotate("stage", geom_boxplot(na.rm = TRUE)),
    'Subtype'=list(
    aes=aes(x=intersection, fill=subtype),
    geom=list(
        geom_bar(stat='count', position='fill', na.rm=TRUE),
        # geom_text(
        #     aes(
        #         label=!!aes_percentage(relative_to='intersection'),
        #         color=ifelse(subtype == '0', 'show', 'hide')
        #     ),
        #     stat='count',
        #     position=position_fill(vjust = .5)
        # ),
        scale_y_continuous(labels=scales::percent_format()),
        # show_hide_scale,
        rating_scale
    ))
  ),
  base_annotations = list(
    "Intersection Size" = intersection_size(
      text = list(size = 5)  
    ) 
    # + scale_y_continuous(limits = c(0, common_ymax)) 
  ),

  
  themes = upset_modify_themes(list(
    default = theme(
      text = element_text(size = text_size),
      axis.text = element_text(size = text_size),
      axis.title = element_text(size = text_size),
      plot.title = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      legend.title = element_text(size = text_size),
      strip.text = element_text(size = text_size),  
      plot.margin = margin(10, 10, 10, 10)
    ),
    intersections_matrix = theme(
      text = element_text(size = text_size),
      axis.text.y = element_text(size = text_size)
    ),
    overall_sizes = theme(
      axis.text = element_text(size = text_size),
      axis.title = element_text(size = text_size)
    ),
    intersection_size = theme(
      axis.text.x = element_text(size = text_size),  
      axis.text.y = element_text(size = text_size),
      axis.title.x = element_text(size = text_size),
      strip.text = element_text(size = text_size)  
    )
  ))
)+ 
scale_fill_manual(values = c('Disc-first' = '#377EB8', 'Facet-first' = '#FF7F00'), guide = "legend")  # Add color scale here

