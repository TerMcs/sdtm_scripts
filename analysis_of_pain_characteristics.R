
library(lme4)
library(lmerTest)  # for p-values
library(marginaleffects)
library(nnet)
library(car)
library(broom)
library(knitr)
library(lme4)
library(survival)
library(ggcorrplot)
library(performance)
library(ComplexUpset)
library(tidyr)
library(dplyr)
library(broom.mixed)  # For tidy model summaries
library(jtools)
library(ggplot2)
library(broom.mixed)
library(betareg)
library(MatchIt)


psm_matched <- FALSE

if (psm_matched == TRUE) {
  file <- "tmcs/data/subtype_stage_results_PSMatched.csv"
  df <- read.csv(file)
  df <- df[, c("participant_id", "match_id", "subtype", "stage")]  # "group", "stage"
} else {
  file <- "tmcs/outputs/subtype_stage_results.csv" 
  df <- read.csv(file)
  df <- df[, c("participant_id", "subtype", "stage", "group")]
  df$participant_id <- gsub("CMBK-010", "", df$participant_id)
  df$participant_id <- gsub("-", "", df$participant_id)
}

# rename participant_ID to C_ID
colnames(df)[1] <- "C_ID"

load("/tmcs/data/REACH_data.rda")

 variables <- c(
    "C_ID", 
    "C_DMAGE", # Age
    "C_DMGENDER", # Gender
    "C_ANBMISRCL", # BMI
    "C_SUTOBAC12M", # Smoking in the last 12 months
    "C_SUCIG", # How would you describe your tobacco use? 1, Never smoker | 2, Current smoker | 3, Used to smoke, but have now quit

    "C_LBLBPDURM", # Length of low back pain (LBP) (self-reported) 1= < 3 months (ineligible for study) 2= 3 - 6 months 3= 6.1 months - 1 year 4= 1.1 - 5 years 5= > 5 years

    "C_QLPROMPFT", # PROMIS physical function T-score
    "C_QLPROMPFS", # PROMIS Physical Function 6b raw score (range: 6-30, higher=more function)
    "C_PIPROMPIS", # PROMIS Pain Interference SF 4a Score (range: 4-20, higher=more interference)
    "C_PIPROMPIT", # PROMIS Pain Interference SF 4a T-SCORE (range: 41.6-75.6, higher=more interference)
    "C_PIPNIMPSC", # Pain impact score (range: 8=least impact to 50=greatest impact)

    "C_PMNEUROSUM", # neuropathic pain (pain detect) - I don't seem to have this variable...

    "C_PIPEGSCORE", # PEG Score-3 Item (Pain, Enjoyment of Life and General Activity scale) (range: 0-10, higher=more pain interference)

    "C_PMMCGSENPS", # McGill Sensory Pain score (range: 0-10, higher=more pain)[Short Form McGill Pain Questionnaire (SF-MPQ)]
    "C_PMMCGAFFPS", # McGill Affective Pain score (range: 0-10, higher=more pain)[Short Form McGill Pain Questionnaire (SF-MPQ)]
 
    "C_PMWSCHPNST", # Has widespread chronic pain 0 or 1

    "C_BEFABQS", # fear avoidance beliefs questionnaire

    "C_PMPNMECCAT", # Pain mechanism category (1 =  nociceptive, 2 = neuropathic, 3 = nociplastic, 4 = mixed)

    "C_PMCNPPT", # pressure pain threshold (PPT) on upper trapezius, a control site on the body
    "C_PMBPPPT", # pressure pain threshold (PPT) on the low back pain (LBP) site on the body
    "C_PMPPTDIF", # the difference between LBP site PPT and control site PPT; more negative, more sensitive to pressure on the back
    "C_PMTSDIF", # the difference between LBP site TS and control site TS; more positive, more summation on the back
    "C_PMAFTERDIF", # the difference between "pressure pain threshold (PPT) on upper trapezius, a control site on the body" & "prepain pain threshold on upper trapezius after hand immersion", more positive, more conditioned pain modulation

    "C_IMFFAL51MU", # multifidus  L5S1 
    "C_IMFFAL51PS", # psoas  L5S1 
    "C_IMFFAL45ES", # erector spinae  L4L5 
    "C_IMFFAL45MU", # multifidus  L4L5 
    "C_IMFFAL45PS", # psoas  L4L5 
    "C_IMFFAL34ES", # erector spinae  L3L4 
    "C_IMFFAL34MU", # multifidus  L3L4 
    "C_IMFFAL34PS", # psoas  L3L4 
    "C_IMFFAL23ES", # erector spinae  L2L3 
    "C_IMFFAL23MU", # multifidus  L2L3 
    "C_IMFFAL23PS", # psoas  L2L3 
    "C_IMFFAL12ES", # erector spinae  L1L2 
    "C_IMFFAL12MU", # multifidus  L1L2 
    "C_IMFFAL12PS", # psoas  L1L2 

    "C_PFBALRT", # Standing Balance Test: single limb: right leg, time held (seconds)
    "C_PFBALLT", # Standing Balance Test: single limb: left leg, time held (seconds)
    "C_PFPRONEINS" # Prone Instability Test
    )


df <- merge(df, df_all, by = "C_ID")
df$C_SUCIG <- factor(df$C_SUCIG)
df$C_DMGENDER <- factor(df$C_DMGENDER)
df <- df[df$C_DMGENDER != 3, ]
df <- df[df$C_DMGENDER != 5, ]
levels_to_drop <- c("3", "5")
df$C_DMGENDER <- droplevels(df$C_DMGENDER, levels = setdiff(levels(df$C_DMGENDER), levels_to_drop))

#### Rename columns ####
colnames(df)[colnames(df) == "C_DMAGE"] <- "Age"
colnames(df)[colnames(df) == "C_DMGENDER"] <- "Gender"
colnames(df)[colnames(df) == "C_ANBMISRCL"] <- "BMI"
colnames(df)[colnames(df) == "C_SUCIG"] <- "Smoking"
colnames(df)[colnames(df) == "stage"] <- "Stage"
colnames(df)[colnames(df) == "C_LBLBPDURM"] <- "LBP duration"
colnames(df)[colnames(df) == "C_QLPROMPFT"] <- "PROMIS physical function T-score"
colnames(df)[colnames(df) == "C_PIPROMPIT"] <- "PROMIS pain interference T-score"
colnames(df)[colnames(df) == "C_PIPEGSCORE"] <- "PEG-3"
colnames(df)[colnames(df) == "C_PIPNIMPSC"] <- "Pain impact score"
colnames(df)[colnames(df) == "C_PMMCGSENPS"] <- "McGill sensory pain"
colnames(df)[colnames(df) == "C_PMMCGAFFPS"] <- "McGill affective pain"
colnames(df)[colnames(df) == "C_PMWSCHPNST"] <- "Chronic widespread pain"
colnames(df)[colnames(df) == "C_BEFABQS"] <- "FABQ"
colnames(df)[colnames(df) == "C_PMPNMECCAT"] <- "Pain mechanism category"
colnames(df)[colnames(df) == "C_PMCNPPT"] <- "PPT control (trapezius)"
colnames(df)[colnames(df) == "C_PMBPPPT"] <- "PPT LBP site"
colnames(df)[colnames(df) == "C_PMPPTDIF"] <- "PPT difference LBP site - control"
colnames(df)[colnames(df) == "C_PMTSDIF"] <- "TS difference LBP site - control"
colnames(df)[colnames(df) == "C_PMAFTERDIF"] <- "PPT trapezius difference after immersion"
colnames(df)[colnames(df) == "C_PFPRONEINS"] <- "Prone instability test"

df$`Nociceptive pain` <- ifelse(df$`Pain mechanism category` == 1, 1, 0)
df$`Nociceptive pain` <- as.factor(df$`Nociceptive pain`)
df$`Neuropathic pain` <- ifelse(df$`Pain mechanism category` == 2, 1, 0)
df$`Neuropathic pain` <- as.factor(df$`Neuropathic pain`)
df$`Nociplastic pain` <- ifelse(df$`Pain mechanism category` == 3, 1, 0)
df$`Nociplastic pain` <- as.factor(df$`Nociplastic pain`)
df$`Mixed pain` <- ifelse(df$`Pain mechanism category` == 4, 1, 0)
df$`Mixed pain` <- as.factor(df$`Mixed pain`)
df$`Neuropathic and mixed pain` <- ifelse(df$`Pain mechanism category` == 2 | df$`Pain mechanism category` == 4, 1, 0)
df$`Chronic widespread pain` <- as.factor(df$`Chronic widespread pain`)

##### PSM analysis #####
df <- df[complete.cases(df[, c(
  "FABQ", "Chronic widespread pain", "Neuropathic and mixed pain", "Pain impact score",
  "McGill sensory pain", "McGill affective pain",
    "PEG-3",
  "PROMIS physical function T-score",
  "PROMIS pain interference T-score" 
                                )]), ]

glm_model <- glm(subtype ~ 
                  Age 
                  + Gender 
                  + BMI 
                  + Smoking , data=df, family=binomial())
ps <- predict(glm_model, type="response")
logit_ps <- log(ps / (1 - ps))
caliper_bound <- 0.2 * sd(logit_ps)

m.out <- matchit(subtype ~ 
                  Age 
                  + Gender 
                  + BMI 
                  + Smoking 
                  ,
                 data = df,
                 method = "nearest",       
                 ratio = 2,                
                 caliper = caliper_bound,            
                 replace = FALSE)           

matched_data <- match.data(m.out)

summary(m.out)

################################## test PCA of subset of variables ########################################################
predictors <- matched_data[, c(
  "FABQ", 
  # "Chronic widespread pain", 
  "Pain impact score",
  "PEG-3",
  "PROMIS physical function T-score",
  "PROMIS pain interference T-score", 
  "McGill sensory pain", "McGill affective pain"
  )]

predictors <- predictors %>%
  mutate(across(everything(), as.numeric))

predictors_scaled <- scale(predictors)

pca_results <- prcomp(predictors_scaled, center = TRUE, scale. = TRUE)

summary(pca_results)

explained_var <- summary(pca_results)$importance[2,]
cum_var <- cumsum(explained_var)
num_components <- which(cum_var >= 0.80)[1]
print(paste("Number of components to keep:", num_components))

pca_scores <- pca_results$x[, 1:num_components]

matched_data$PC1 <- pca_scores[, 1]
if (num_components >= 2) matched_data$PC2 <- pca_scores[, 2]
if (num_components >= 3) matched_data$PC3 <- pca_scores[, 3]
if (num_components >= 4) matched_data$PC4 <- pca_scores[, 4]

######## Logistic regression with random effects for subtype ####################################################

model1 <- glmer(subtype ~ 
                  Age 
                  + Gender 
                  + BMI 
                  + Smoking 
                  + Stage
                 + `PROMIS physical function T-score`
                 + `PROMIS pain interference T-score`
                  + `PEG-3` 
                  + `McGill sensory pain`
                  + `McGill affective pain`
                  + FABQ
                  + `Chronic widespread pain`
                 + `Nociceptive pain`
                  + `Neuropathic and mixed pain`
                 + `Nociplastic pain`
                  + `Pain impact score`
                  + (1 | subclass)
            , data = matched_data, family = binomial)

model2 <- glmer(subtype ~ 
                  Age 
                  + Gender 
                  + BMI 
                  + Smoking 
                  + Stage
                 + `PROMIS physical function T-score`
                 + `PROMIS pain interference T-score`
                  + `PEG-3` 
                  + `McGill sensory pain`
                  + `McGill affective pain`
                  + FABQ
                  + `Chronic widespread pain`
                  + `Neuropathic and mixed pain`
                  + `Pain impact score`
                  + (1 | subclass)
            , data = matched_data, family = binomial)

modelPCA <- glmer(subtype ~ 
                  Age 
                  + Gender 
                  + BMI 
                  + Smoking 
                  + Stage
                  + PC1
                 + PC2
                 + PC3
                  + `Chronic widespread pain`
                  + `Neuropathic and mixed pain`
                  + `Pain impact score`
                  + (1 | subclass)
            , data = matched_data, family = binomial)

model3 <- glmer(subtype ~ 
                  Age 
                  + Gender 
                  + BMI 
                  + Smoking 
                  + FABQ
                  + `Chronic widespread pain`
                  + `Neuropathic and mixed pain`
                  + `Pain impact score`
                  + (1 | subclass)
                  
            , data = matched_data, family = binomial)

# compare AIC/BIC for best model:
model_comparison <- bind_rows(
  broom.mixed::glance(model1),
  broom.mixed::glance(model2),
  broom.mixed::glance(modelPCA),
  broom.mixed::glance(model3),
)
print(model_comparison)

# compare nested models using ANOVA
anova(model2, model1, test="Chisq")
anova(model3, model2, test="Chisq")

summary(model3)

shapiro.test(matched_data$`PPT difference LBP site - control`)
shapiro.test(matched_data$`TS difference LBP site - control`)

wilcox.test(matched_data$`PPT difference LBP site - control` ~ matched_data$subtype)
wilcox.test(matched_data$`TS difference LBP site - control` ~ matched_data$subtype)

model <- model3

tidy_model <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)
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
  mutate(OR = exp(estimate),  # Calculate Odds Ratio
        OR_lower = exp(conf.low),  # Lower CI
        OR_upper = exp(conf.high))  # Upper CI 

  write.csv(model_summary, file = "tmcs/outputs/results_PSMatched_glmer.csv", row.names = FALSE)


######## Beta regression for stage ####################################################

epsilon <- 1e-6
df$stage_scaled <- (df$Stage / 21) * (1 - 2 * epsilon) + epsilon
df$stage_scaled[df$stage_scaled == 1] <- df$stage_scaled[df$stage_scaled == 1] - epsilon

df$Age <- as.numeric(as.character(df$Age))
str(df$Age)
unique(df$Age)
df$subtype <- as.factor(df$subtype)

beta_model_betareg <- betareg(stage_scaled ~ 
                    Age
                  + Gender 
                  + BMI 
                  + Smoking 
                   + (FABQ
                   + `Chronic widespread pain`
                   + `Neuropathic and mixed pain`
                   + `Pain impact score`
                   ) * subtype,
                  data = df)
summary(beta_model_betareg)
plot(residuals(beta_model_betareg) ~ fitted(beta_model_betareg))

par(mfrow = c(3, 2))
plot(beta_model_betareg, which = 1:4, type = "pearson")
plot(beta_model_betareg, which = 5, type = "deviance")
plot(beta_model_betareg, which = 1, type = "deviance")

# Diagnostic Plots
plot(beta_model_betareg, which = 1:3)  # Examine residual plots

lrtest(model, model_phi)

plot(beta_model_betareg, which = 1:3) 
plot(residuals(beta_model_betareg) ~ fitted(beta_model_betareg))


ggplot(df, aes(x = `Pain impact score`, y = stage_scaled, color = subtype)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x) +  
    facet_wrap(~ subtype) 

plot_stage_scaled_by_binary_and_subtype <- function(data, binary_variable, subtype_variable = "subtype", stage_scaled_variable = "stage_scaled") {
  ggplot(data, aes(x = factor(.data[[binary_variable]]), y = .data[[stage_scaled_variable]], fill = .data[[binary_variable]])) +  # Use factor()
    geom_boxplot() +
    facet_wrap(~ .data[[subtype_variable]]) +
    labs(
      title = paste("Boxplot of", stage_scaled_variable, "by", binary_variable, "and", subtype_variable),
      x = binary_variable,
      y = stage_scaled_variable,
      fill = binary_variable
    ) +
    scale_x_discrete(labels = c("0", "1")) + 
    theme_bw() +
    theme(legend.position = "bottom")


model <- beta_model_betareg
estimates <- coef(model, model = "mean")
std_errors <- summary(model)$coefficients$mean[, "Std. Error"]
z_values <- summary(model)$coefficients$mean[, "z value"]
p_values <- summary(model)$coefficients$mean[, "Pr(>|z|)"]

odds_ratios <- exp(estimates)

z_critical <- qnorm(0.975)  # For 95% CI
lower_ci <- estimates - (z_critical * std_errors)
upper_ci <- estimates + (z_critical * std_errors)
odds_lower_ci <- exp(lower_ci)
odds_upper_ci <- exp(upper_ci)

# Extract Phi parameters
phi_estimates <- coef(model, model = "precision")
phi_std_errors <- summary(model)$coefficients$precision[, "Std. Error"]
phi_z_values <- summary(model)$coefficients$precision[, "z value"]
phi_p_values <- summary(model)$coefficients$precision[, "Pr(>|z|)"]
phi_lower_ci <- phi_estimates - (z_critical * phi_std_errors)
phi_upper_ci <- phi_estimates + (z_critical * phi_std_errors)

results_df <- data.frame(
    Predictor = names(estimates),
    Estimate = estimates,
    Std.Error = std_errors,
    z.value = z_values,
    p.value = p_values,
    Odds.Ratio = odds_ratios,
    CI.Lower = lower_ci,
    CI.Upper = upper_ci,
    Odds.CI.Lower = odds_lower_ci,
    Odds.CI.Upper = odds_upper_ci
)

write.csv(results_df, file = "tmcs/outputs/results_betareg_stage_scaled.csv", row.names = FALSE)
}

#### Plot BMI by age group ###########################################################################################################################

df_subtype <- df[df$subtype == 1, ]

df_subtype$age_group <- cut_number(df_subtype$Age, n = 2)
# # Check the updated labels
levels(df_subtype$age_group)
levels(df_subtype$age_group) <- c("20-63", "63-84")
levels(df_subtype$age_group)



df_no_controls$age_group <- cut_number(df_no_controls$Age, n = 2)
# # Check the updated labels
levels(df_no_controls$age_group)
levels(df_no_controls$age_group) <- c("20-59", "60-91")
levels(df_no_controls$age_group)

df_subtype <- df_no_controls[df_no_controls$subtype == 0, ]

# Fit model per age_group and extract R²
stats_df <- df_subtype %>%
  group_by(age_group) %>%
  do({
    mod <- lm(BMI ~ stage, data = .)
    glance_mod <- glance(mod)
    tidy_mod <- tidy(mod)

    data.frame(
      n = nrow(.),
      r_squared = round(glance_mod$r.squared, 3),
      p_value = round(tidy_mod$p.value[tidy_mod$term == "stage"], 3)
    )
  }) %>%
  mutate(
    facet_label = paste0("Age group: ", age_group, ", N = ", n,
                        "\nR² = ", r_squared, ", p = ", p_value)
  )

stats_df$facet_label <- factor(stats_df$facet_label,
                              levels = stats_df$facet_label[order(stats_df$age_group)])

df_plot <- df_subtype %>%
  left_join(stats_df, by = "age_group")

ggplot(df_plot, aes(x = stage, y = BMI)) + # , color = subtype
  geom_jitter(alpha = 0.6, width = 2, height = 0) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) + # aes(group = subtype), 1
  facet_wrap(~ facet_label) +
  labs(
    title = "BMI vs stage by age group in facet-first subtype",
    x = "Stage", y = "BMI", color = "Subtype"
  ) +
  scale_y_continuous(limits = c(15, NA)) +
  theme_minimal()
