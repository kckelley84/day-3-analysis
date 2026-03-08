# Library -----------------------------------------------------------------

library(tidyverse)
library(glmmTMB)
library(marginaleffects)

# Data Read ---------------------------------------------------------------

main_data <- read.csv("clean_data.csv", header = TRUE, na.strings = c("", "NA"))

glimpse(main_data)

# Create/modify new variables as needed

main_data <- main_data %>%
  mutate(
    eq_score_beta = eq_score_beta / 10,
    tEB_reached = ifelse(!is.na(tEB), 1L, 0L),
    tB_reached = ifelse(!is.na(tB), 1L, 0L),
    icm_num = case_when(
      gard_icm == "A" ~ 3L,
      gard_icm == "B" ~ 2L,
      gard_icm == "C" ~ 1L,
      TRUE ~ NA_integer_
    ),
    te_num = case_when(
      gard_te == "A" ~ 3L,
      gard_te == "B" ~ 2L,
      gard_te == "C" ~ 1L,
      TRUE ~ NA_integer_
    ),
    gardner_score = if_else(
      is.na(gard_exp) | is.na(icm_num) | is.na(te_num),
      NA_real_,
      gard_exp * 100 + icm_num * 10 + te_num
    )
  )

summary(main_data$eq_score_beta)


# Demographics ------------------------------------------------------------

main_data %>%
  filter(!is.na(pgt_status)) %>%
  group_by(patient_id) %>%
  summarise(
    n_embryos = n(),
    .groups = "drop"
  ) %>%
  summarise(
    mean_embryos_per_patient = mean(n_embryos),
    sd_embryos_per_patient   = sd(n_embryos),
    median_embryos           = median(n_embryos),
    min_embryos              = min(n_embryos),
    max_embryos              = max(n_embryos)
  )

main_data %>%
  filter(!is.na(gard_exp)) %>%
  group_by(d3_hatching) %>%
  summarize(
    n = n(),
    mean_age = mean(mat_age, na.rm = TRUE),
    mean_eq = mean(eq_score, na.rm = TRUE),
    mean_blast = mean(blast_score, na.rm = TRUE),
    
  )
# EQ Score - d3 AH --------------------------------------------------

eq_ah_train <- main_data %>%
  filter(
    split == "train",
    !is.na(eq_score_beta),
    !is.na(d3_hatching),
    !is.na(t8)
  )

nrow(eq_ah_train)

eq_ah_train %>%
  group_by(d3_hatching) %>%
  summarise(
    n = n(),
    mean_eq = mean(eq_score_beta, na.rm = TRUE),
    sd_eq   = sd(eq_score_beta, na.rm = TRUE)
  )

hist(
  eq_ah_train$eq_score_beta,
  breaks = 40,
  main = "EQ score distribution",
  xlab = "EQ score"
)

eq_ah_model <- glmmTMB(
  eq_score_beta ~ d3_hatching + mat_age + (1 | patient_id),
  data = eq_ah_train,
  family = beta_family(link = "logit")
)

summary(eq_ah_model)

beta <- summary(eq_ah_model)$coefficients$cond["d3_hatching", "Estimate"]
se   <- summary(eq_ah_model)$coefficients$cond["d3_hatching", "Std. Error"]

or <- exp(beta)
ci <- exp(beta + c(-1.96, 1.96) * se)

or
ci

eq_ah_test <- main_data %>%
  filter(
    split == "test",
    !is.na(eq_score_beta),
    !is.na(d3_hatching),
    !is.na(tSB)
  ) %>%
  mutate(
    pred_eq = predict(eq_ah_model, newdata = ., type = "response")
  )

rmse_eq <- sqrt(mean((eq_ah_test$eq_score_beta - eq_ah_test$pred_eq)^2))
cor_eq  <- cor(eq_ah_test$eq_score_beta, eq_ah_test$pred_eq)

rmse_eq
cor_eq


# tSB and D3 AH -----------------------------------------------------

tsb_ah_train <- main_data %>%
  filter(
    split == "train",
    !is.na(tSB),
    !is.na(d3_hatching)
  )

nrow(tsb_ah_train)

tsb_ah_means <- tsb_ah_train %>%
  group_by(d3_hatching) %>%
  summarise(
    n = n(),
    mean_tsb = round(mean(tSB, na.rm = TRUE), 2),
    sd_tsb   = round(sd(tSB, na.rm = TRUE), 2)
  )

tsb_ah_means

hist(
  tsb_ah_train$tSB,
  breaks = 40,
  main = "tSB distribution",
  xlab = "tSB score"
)

tSB_ah_model <- glmmTMB(
  tSB ~ d3_hatching + mat_age + (1 | patient_id),
  data = eq_ah_train
)

summary(tSB_ah_model)

# tB and D3 AH -----------------------------------------------------

tB_ah_train <- main_data %>%
  filter(
    split == "train",
    !is.na(tB),
    !is.na(d3_hatching)
  )

nrow(tB_ah_train)

tB_ah_means <- tB_ah_train %>%
  group_by(d3_hatching) %>%
  summarise(
    n = n(),
    mean_tB = round(mean(tB, na.rm = TRUE), 2),
    sd_tB   = round(sd(tB, na.rm = TRUE), 2)
  )

tB_ah_means

hist(
  tB_ah_train$tB,
  breaks = 40,
  main = "tB distribution",
  xlab = "tB score"
)

tB_ah_model <- glmmTMB(
  tB ~ d3_hatching + mat_age + (1 | patient_id),
  data = eq_ah_train
)

summary(tB_ah_model)

# tEB and D3 AH -----------------------------------------------------

tEB_ah_train <- main_data %>%
  filter(
    split == "train",
    !is.na(tEB),
    !is.na(d3_hatching)
  )

nrow(tEB_ah_train)

tEB_ah_means <- tEB_ah_train %>%
  group_by(d3_hatching) %>%
  summarise(
    n = n(),
    mean_tEB = round(mean(tEB, na.rm = TRUE), 2),
    sd_tEB   = round(sd(tEB, na.rm = TRUE), 2)
  )

tEB_ah_means

hist(
  tEB_ah_train$tEB,
  breaks = 40,
  main = "tEB distribution",
  xlab = "tEB score"
)

tEB_ah_model <- glmmTMB(
  tEB ~ d3_hatching + mat_age + (1 | patient_id),
  data = eq_ah_train
)

summary(tEB_ah_model)

# Ci

beta  <- fixef(tEB_ah_model)$cond["d3_hatching"]
se    <- sqrt(vcov(tEB_ah_model)$cond["d3_hatching","d3_hatching"])

ci <- beta + c(-1.96, 1.96) * se

tEB_ah_test <- main_data %>%
  filter(
    split == "test",
    !is.na(d3_hatching),
    !is.na(tEB)
  ) %>%
  mutate(
    pred_tEB = predict(tEB_ah_model, newdata = ., type = "response")
  )

rmse_tEB <- sqrt(mean((tEB_ah_test$tEB - tEB_ah_test$pred_tEB)^2))
cor_tEB  <- cor(tEB_ah_test$tEB, tEB_ah_test$pred_tEB)

rmse_tEB
cor_tEB


# tEB Reached and AH ------------------------------------------------

tb_no_teb <- main_data %>%
  filter(
    tB_reached == 1,
    is.na(tEB)
  )

tb_no_teb %>%
  count(emb_fate)

tb_no_teb %>%
  count(emb_fate, d3_hatching)

main_data <- main_data %>%
  mutate(
    tEB_viable = case_when(
      tEB_reached == 1 ~ 1,
      tB_reached == 1 & is.na(tEB) & emb_fate == "FROZEN" ~ 1,
      TRUE ~ 0
    )
  )

table(main_data$tEB_reached, main_data$tEB_viable, useNA = "ifany")

tEB_viable_model <- glmmTMB(
  tEB_viable ~ d3_hatching + mat_age + (1 | patient_id),
  data = main_data %>%
    filter(!is.na(d3_hatching), !is.na(t8)),
  family = binomial()
)

summary(tEB_viable_model)

beta <- summary(tEB_viable_model)$coefficients$cond["d3_hatching", "Estimate"]
se   <- summary(tEB_viable_model)$coefficients$cond["d3_hatching", "Std. Error"]

or <- exp(beta)
ci <- exp(beta + c(-1.96, 1.96) * se)

or
ci

main_data %>%
  filter(!is.na(d3_hatching), !is.na(t8)) %>%
  group_by(d3_hatching) %>%
  summarise(
    n = n(),
    viable = sum(tEB_viable == 1),
    rate = viable / n
  )

# AH and PGT --------------------------------------------------------

main_data %>%
  group_by(pgt_group)%>%
  summarise(
    n = n()
  )

pgt_ah_train <- main_data %>%
  filter(
    split == "train",
    !is.na(d3_hatching),
    pgt_group %in% c("euploid", "aneuploid")
  ) %>%
  mutate(
    pgt_euploid = ifelse(pgt_group == "euploid", 1L, 0L)
    )

pgt_ah_model <- glmmTMB(
  pgt_euploid ~ d3_hatching + mat_age + (1 | patient_id),
  data = pgt_ah_train,
  family = binomial()
)

summary(pgt_ah_model)

pgt_ah_train %>%
  group_by(d3_hatching) %>%
  summarise(
    n = n(),
    euploid_rate = mean(pgt_euploid),
    euploid_n = sum(pgt_euploid),
    aneuploid_n = n() - euploid_n
  )

beta <- summary(pgt_ah_model)$coefficients$cond["d3_hatching", "Estimate"]
se   <- summary(pgt_ah_model)$coefficients$cond["d3_hatching", "Std. Error"]

or <- exp(beta)
ci <- exp(beta + c(-1.96, 1.96) * se)

or
ci

# BLAST Score - d3 AH --------------------------------------------------

blast_ah_train <- main_data %>%
  filter(
    split == "train",
    !is.na(blast_score_beta),
    !is.na(d3_hatching),
    !is.na(t8)
  )

nrow(blast_ah_train)

blast_ah_train %>%
  group_by(d3_hatching) %>%
  summarise(
    n = n(),
    mean_blast = mean(blast_score_beta, na.rm = TRUE),
    sd_blast   = sd(blast_score_beta, na.rm = TRUE)
  )

hist(
  blast_ah_train$blast_score_beta,
  breaks = 40,
  main = "blast score distribution",
  xlab = "blast score"
)

blast_ah_model <- glmmTMB(
  blast_score_beta ~ d3_hatching + mat_age + (1 | patient_id),
  data = blast_ah_train,
  family = gaussian()
)

summary(blast_ah_model)

beta <- summary(blast_ah_model)$coefficients$cond["d3_hatching", "Estimate"]
se   <- summary(blast_ah_model)$coefficients$cond["d3_hatching", "Std. Error"]

or <- exp(beta)
ci <- exp(beta + c(-1.96, 1.96) * se)

or
ci

# Gardner Score - d3 AH --------------------------------------------------

gard_ah_train <- main_data %>%
  filter(
    split == "train",
    !is.na(gardner_score),
    !is.na(d3_hatching),
    !is.na(t8)
  )

nrow(gard_ah_train)

gard_ah_train %>%
  group_by(d3_hatching) %>%
  summarise(
    n = n(),
    mean_gard = mean(gardner_score, na.rm = TRUE),
    sd_gard   = sd(gardner_score, na.rm = TRUE)
  )

hist(
  gard_ah_train$gardner_score,
  breaks = 40,
  main = "gardner score distribution",
  xlab = "gardner score"
)

gard_ah_model <- glmmTMB(
  gardner_score ~ d3_hatching + mat_age + (1 | patient_id),
  data = gard_ah_train,
  family = gaussian()
)

summary(gard_ah_model)

# Gardner Component Analysis

m_exp <- glmmTMB(
  gard_exp ~ d3_hatching + mat_age + (1 | patient_id),
  data = gard_ah_train,
  family = gaussian()
)
summary(m_exp)

m_icm <- glmmTMB(
  icm_num ~ d3_hatching + mat_age + (1 | patient_id),
  data = gard_ah_train,
  family = gaussian()
)
summary(m_icm)

m_te <- glmmTMB(
  te_num ~ d3_hatching + mat_age + (1 | patient_id),
  data = gard_ah_train,
  family = gaussian()
)
summary(m_te)

