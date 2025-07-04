# Required packages
library(dplyr)
library(purrr)
library(ggplot2)
library(rcompanion)  
library(tidyr) 
library(knitr)
library(ggcorrplot)
library(lme4)
library(broom.mixed)  # For tidy model summaries
library(jtools)
library(tidyverse)


file <- "/tmcs/outputs/strict_ivd_fjoa_controls_included_trained_model_ALL_DATA.csv"
df <- read.csv(file)

all_binary_vars <- c(
  "disc_af_12_nl", "disc_af_23_nl", "disc_af_34_nl", "disc_af_45_nl", "disc_af_51_nl", 
  "ep_defect_l2_sup_shape", "ep_defect_l3_sup_shape", "ep_defect_l4_sup_shape", "ep_defect_l5_sup_shape", "ep_defect_s1_sup_shape", 
  "ep_defect_l1_inf_shape", "ep_defect_l2_inf_shape", "ep_defect_l3_inf_shape", "ep_defect_l4_inf_shape", "ep_defect_l5_inf_shape",
  "stenosis_congenital_narrowing_spinal_canal_12_nl", "stenosis_congenital_narrowing_spinal_canal_23_nl", "stenosis_congenital_narrowing_spinal_canal_34_nl", "stenosis_congenital_narrowing_spinal_canal_45_nl", "stenosis_congenital_narrowing_spinal_canal_51_nl", 
  "stenosis_lumbo_sacral_seg_abnormality"
)    
all_ordinal_vars <- c(
  "disc_d_herniation_12_nl", "disc_d_herniation_23_nl", "disc_d_herniation_34_nl", "disc_d_herniation_45_nl", "disc_d_herniation_51_nl",
  "disc_nr_involvement_12_nl", "disc_nr_involvement_23_nl", "disc_nr_involvement_34_nl", "disc_nr_involvement_45_nl", "disc_nr_involvement_51_nl",
  "stenosis_ccs_12_nl", "stenosis_ccs_23_nl", "stenosis_ccs_34_nl", "stenosis_ccs_45_nl", "stenosis_ccs_51_nl", 
  "stenosis_lrs_12_nl", "stenosis_lrs_23_nl", "stenosis_lrs_34_nl", "stenosis_lrs_45_nl", "stenosis_lrs_51_nl",
  "stenosis_fs_12_nl", "stenosis_fs_23_nl", "stenosis_fs_34_nl", "stenosis_fs_45_nl", "stenosis_fs_51_nl" 
)  

#### fix some variables: convert nerve root involvement to 4 levels:
df <- df %>%
  mutate(across(
    c(disc_nr_involvement_12_nl, disc_nr_involvement_23_nl, disc_nr_involvement_34_nl, 
      disc_nr_involvement_45_nl, disc_nr_involvement_51_nl),
    ~ case_when(
      .x == 3 ~ 1L,
      .x == 4 ~ 2L,
      .x == 5 ~ 3L,
      TRUE ~ .x
    )
  ))

#### sum nerve root involvements:
df$summed_nerve_root_involvement <- rowSums(df[, c(
  "disc_nr_involvement_12_nl", 
  "disc_nr_involvement_23_nl", 
  "disc_nr_involvement_34_nl", 
  "disc_nr_involvement_45_nl", 
  "disc_nr_involvement_51_nl"
)], na.rm = TRUE)

#### binarise LSTV:
df <- df %>%
  mutate(stenosis_lumbo_sacral_seg_abnormality = ifelse(stenosis_lumbo_sacral_seg_abnormality == 0, 0, 1)) 

#### select only schmorl"s nodes and sum them:
schmorls_nodes <- c(
  "ep_defect_l2_sup_shape", 
  "ep_defect_l3_sup_shape", 
  "ep_defect_l4_sup_shape", 
  "ep_defect_l5_sup_shape", 
  "ep_defect_s1_sup_shape", 
  "ep_defect_l1_inf_shape", 
  "ep_defect_l2_inf_shape", 
  "ep_defect_l3_inf_shape", 
  "ep_defect_l4_inf_shape", 
  "ep_defect_l5_inf_shape"
)

df <- df %>%
  mutate(across(all_of(schmorls_nodes), ~ ifelse(.x == 2, 1, 0)))

df$summed_schmorls_nodes <- rowSums(df[, schmorls_nodes], na.rm = TRUE)

#### combine the protrusion and extrusion categories ( "disc_d_herniation_12_nl" = 3 gets converted to 2):
df <- df %>%
  mutate(across(
    c(disc_d_herniation_12_nl, disc_d_herniation_23_nl, disc_d_herniation_34_nl, disc_d_herniation_45_nl, disc_d_herniation_51_nl),
    ~ case_when(
      .x == 3 ~ 2L,
      TRUE ~ .x
    )
  ))

#### sum disc bulges:
df$summed_disc_herniation <- rowSums(df[, c(
  "disc_d_herniation_12_nl", 
  "disc_d_herniation_23_nl", 
  "disc_d_herniation_34_nl", 
  "disc_d_herniation_45_nl", 
  "disc_d_herniation_51_nl"
)], na.rm = TRUE)

#### Sum annular fissures:
df$summed_annular_fissures <- rowSums(df[, c(
  "disc_af_12_nl",
  "disc_af_23_nl",
  "disc_af_34_nl", 
  "disc_af_45_nl",
  "disc_af_51_nl"
)], na.rm = TRUE)

#### sum stenosis   "stenosis_congenital_narrowing_spinal_canal_12_nl", "stenosis_congenital_narrowing_spinal_canal_23_nl", "stenosis_congenital_narrowing_spinal_canal_34_nl", "stenosis_congenital_narrowing_spinal_canal_45_nl", "stenosis_congenital_narrowing_spinal_canal_51_nl", 
df$summed_stenosis_congenital_narrowing_spinal_canal <- rowSums(df[, c(
  "stenosis_congenital_narrowing_spinal_canal_12_nl",
  "stenosis_congenital_narrowing_spinal_canal_23_nl",
  "stenosis_congenital_narrowing_spinal_canal_34_nl", 
  "stenosis_congenital_narrowing_spinal_canal_45_nl",
  "stenosis_congenital_narrowing_spinal_canal_51_nl"
)], na.rm = TRUE)

#### sum stenosis "stenosis_ccs_12_nl", "stenosis_ccs_23_nl", "stenosis_ccs_34_nl", "stenosis_ccs_45_nl", "stenosis_ccs_51_nl",
df$summed_stenosis_ccs <- rowSums(df[, c(
  "stenosis_ccs_12_nl",
  "stenosis_ccs_23_nl",
  "stenosis_ccs_34_nl", 
  "stenosis_ccs_45_nl",
  "stenosis_ccs_51_nl"
)], na.rm = TRUE)

#### convert    "stenosis_lrs_lrb_12_nl" = 0 to 1
df <- df %>%
  mutate(across(
    c(stenosis_lrs_lrb_12_nl, stenosis_lrs_lrb_23_nl, stenosis_lrs_lrb_34_nl, stenosis_lrs_lrb_45_nl, stenosis_lrs_lrb_51_nl),
    ~ case_when(
      .x == 0 ~ 1L,
      TRUE ~ .x
    )
  ))                               

#### multiply stenosis_lrs_lrb_12_nl by stenosis_lrs_12_nl
lrb_cols <- c(
  "stenosis_lrs_lrb_12_nl",
  "stenosis_lrs_lrb_23_nl", 
  "stenosis_lrs_lrb_34_nl", 
  "stenosis_lrs_lrb_45_nl", 
  "stenosis_lrs_lrb_51_nl"
)
lrs_cols <- paste0("stenosis_lrs_", c("12_nl", "23_nl", "34_nl", "45_nl", "51_nl"))

# Multiply matching columns
df <- df %>%
  mutate(across(all_of(lrb_cols), ~ .x * df[[lrs_cols[which(lrb_cols == cur_column())]]]))

#### sum stenosis "stenosis_lrs_12_nl", "stenosis_lrs_23_nl", "stenosis_lrs_34_nl", "stenosis_lrs_45_nl", "stenosis_lrs_51_nl",
df$summed_lateral_recess_stenosis <- rowSums(df[, c(
  "stenosis_lrs_12_nl", 
  "stenosis_lrs_23_nl", 
  "stenosis_lrs_34_nl", 
  "stenosis_lrs_45_nl",
  "stenosis_lrs_51_nl"
)], na.rm = TRUE)

#### multiply stenosis_fs_lrb_12_nl by stenosis_fs_12_nl
lrb_cols <- c(
  "stenosis_fs_lrb_12_nl", 
  "stenosis_fs_lrb_23_nl", 
  "stenosis_fs_lrb_34_nl", 
  "stenosis_fs_lrb_45_nl", 
  "stenosis_fs_lrb_51_nl"
)
fs_cols <- paste0("stenosis_fs_", c("12_nl", "23_nl", "34_nl", "45_nl", "51_nl"))

# Multiply matching columns
df <- df %>%
  mutate(across(all_of(lrb_cols), ~ .x * df[[lrs_cols[which(lrb_cols == cur_column())]]]))

#### sum stenosis "stenosis_fs_12_nl", "stenosis_fs_23_nl", "stenosis_fs_34_nl", "stenosis_fs_45_nl", "stenosis_fs_51_nl",
df$summed_facet_stenosis <- rowSums(df[, c(
  "stenosis_fs_12_nl", 
  "stenosis_fs_23_nl", 
  "stenosis_fs_34_nl", 
  "stenosis_fs_45_nl", 
  "stenosis_fs_51_nl"
)], na.rm = TRUE)

#### sum summed_lateral_recess_stenosis and summed_stenosis_ccs
df$summed_intervertebral_foramen_stenosis <- df$summed_lateral_recess_stenosis + df$summed_stenosis_ccs

all_var_names <- c(all_binary_vars, all_ordinal_vars)

composite_var_names <- c(
  "summed_nerve_root_involvement", 
  "summed_disc_herniation", 
  "summed_annular_fissures", 
  "summed_schmorls_nodes", 
  # summed_stenosis_congenital_narrowing_spinal_canal + 
  "summed_stenosis_ccs",
  # summed_lateral_recess_stenosis + 
  # summed_facet_stenosis +
  "summed_intervertebral_foramen_stenosis"
  # "stenosis_lumbo_sacral_seg_abnormality" 
)

######### merge with reach server data
# rename participant_ID to C_ID
colnames(df)[1] <- "C_ID"

load("/tmcs/data/REACH_data.rda")

 variables <- c(
    "C_ID", 
    "C_DMAGE", # Age
    "C_DMGENDER", # Gender
    "C_ANBMISRCL", # BMI
    "C_SUCIG" # How would you describe your tobacco use? 1, Never smoker | 2, Current smoker | 3, Used to smoke, but have now quit
    )

# merge the data on ID
df <- merge(df, df_all, by = "C_ID")
df$C_SUCIG <- ifelse(df$C_SUCIG == 3, 2, df$C_SUCIG)
df$C_SUCIG <- factor(df$C_SUCIG)
df$C_DMGENDER <- factor(df$C_DMGENDER)
df <- df[df$C_DMGENDER != 3, ]
df <- df[df$C_DMGENDER != 5, ]

colnames(df)[colnames(df) == "C_DMAGE"] <- "Age"
colnames(df)[colnames(df) == "C_DMGENDER"] <- "Gender"
colnames(df)[colnames(df) == "C_ANBMISRCL"] <- "BMI"
colnames(df)[colnames(df) == "C_SUCIG"] <- "Smoking"

df$Smoking <- as.character(df$Smoking)  
df$Smoking[is.na(df$Smoking)] <- "1"    
df$Smoking <- as.factor(df$Smoking)  

df <- df %>% filter(!is.na(match_id))

model <- glmer(
  subtype ~ 
    Age +
      Gender +
      BMI +
      Smoking +
      summed_nerve_root_involvement + 
      summed_disc_herniation + 
      summed_annular_fissures + 
      summed_schmorls_nodes + 
      # summed_stenosis_congenital_narrowing_spinal_canal + 
      summed_stenosis_ccs + 
      # summed_lateral_recess_stenosis + 
      # summed_facet_stenosis +
      # summed_intervertebral_foramen_stenosis +
      stenosis_lumbo_sacral_seg_abnormality +
      (1 | match_id),
  data = df,
  family = binomial
)

summary(model)

#### feature correlation matrix:
vars <- model.matrix(model)[, -1] 
cor_matrix <- cor(vars, use = "pairwise.complete.obs")

ggcorrplot(cor_matrix, lab = TRUE, type = "upper")

tidy_model <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)
print(tidy_model)

summ(model, exp = TRUE, scale = TRUE, confint = TRUE)

model_summary <- tidy(model, conf.int = TRUE)
model_summary %>%
  mutate(estimate = round(estimate, 3),
         std.error = round(std.error, 3),
         statistic = round(statistic, 2),
         p.value = round(p.value, 3),
         conf.low = round(conf.low, 3),
         conf.high = round(conf.high, 3)) %>%
  kable(caption = "Logistic Regression Model Results", align = c("l", "r", "r", "r", "r", "r", "r"))
model_summary <- tidy(model, conf.int = TRUE) %>%
  mutate(OR = exp(estimate),  
         OR_lower = exp(conf.low),  
         OR_upper = exp(conf.high))  
write.csv(model_summary, file = "tmcs/data/model_results/subtype_stage_inferences_glmer_other_features.csv", row.names = FALSE)


df <- df %>% filter(!is.na(summed_nerve_root_involvement))

ggplot(df, aes(x = as.factor(summed_nerve_root_involvement), fill = as.factor(subtype))) +
  geom_bar(position = "dodge") +
  # scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_size = 14) +
  labs(x = "Your Variable", y = "Count", fill = "subtype")


ggplot(df, aes(x = as.factor(subtype), y = summed_disc_herniation)) +
  geom_boxplot() +
  theme_minimal()


### check LSTV differences:
data_clean <- df %>% filter(match_id != 2)

# Reshape the data for matched pairs
data_matched <- data_clean %>%
  group_by(match_id) %>%
  reframe(
    facet_first = .data[["stenosis_lumbo_sacral_seg_abnormality"]][subtype == 1],
    disc_first = .data[["stenosis_lumbo_sacral_seg_abnormality"]][subtype == 0]
  ) 

contingency_table <- table(data_matched$disc_first, data_matched$facet_first)

# # # Perform McNemar"s Test
test_result <- mcnemar.test(contingency_table)


#### check for summed / composite variable differences:
results_ordinal <- data.frame(variable = character(),
                      p_value = numeric(),
                      effect_size = numeric(),
                      method = character(),
                      type = character(),
                      stringsAsFactors = FALSE)

for (var_name in composite_var_names) {
  # Reshape the data for matched pairs
  data_matched <- data_clean %>%
    group_by(match_id) %>%
    summarise(
      var_treatment = .data[[var_name]][subtype == 1],
      var_control = .data[[var_name]][subtype == 0]
    ) %>%
    filter(!is.na(var_treatment) & !is.na(var_control))  

  print(length(data_matched$var_treatment))
  print(length(data_matched$var_control))

  r_rc <- wilcox.test(data_matched$var_treatment, data_matched$var_control, paired = TRUE, exact = FALSE)

  differences <- data_matched$var_treatment - data_matched$var_control
  ranks <- rank(differences)

  # Effect size calculation (Rank-Biserial Correlation)
  n <- length(differences)  
  W <- r_rc$statistic 
  effect_size <- (W - n*(n + 1)/2) / (n * (n - 1) / 2)

  result <- (data.frame(
    variable = var_name,
    p_value = r_rc$p.value,
    effect_size = effect_size,  # rank-biserial correlation coefficient
    method = "Wilcoxon Paired Test (Rank Biserial Correlation)",
    type = "ordinal"
  ))

  results_ordinal <- rbind(results_ordinal, result)
}
print(results_ordinal)









