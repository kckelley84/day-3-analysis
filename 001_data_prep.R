
# Library -----------------------------------------------------------------

library(tidyverse)


# Data Read ---------------------------------------------------------------

main_data <- read.csv("raw_data.csv", header = TRUE, na.strings = c("", "NA"))

glimpse(main_data)

# Data Prep ---------------------------------------------------------------

# Split

set.seed(41284)  # for reproducibility

patient_split <- main_data %>%
  distinct(patient.esId) %>%
  mutate(
    split = ifelse(runif(n()) < 0.7, "train", "test")
  )

main_data_split <- main_data %>%
  left_join(patient_split, by = "patient.esId") %>%
  mutate(
    blastulation = ifelse(!is.na(tSB), 1L, 0L)
  )

# Trim Down Dataset

selected_data <- main_data_split %>%
  select(patient.oocyteAge.years, patient.esId, inseminationMethod, status, clinicGardner.expansion,
         clinicGardner.icm, clinicGardner.te, pgta, kidScore, t8, tSB, tB, tEB, blastAtHr, assistedHatching.Day.5..key,
         kidRank, assistedHatching.Day3..key, manualAnnotation.Day.6.assistedHatching.0..key, blastulation, split)

# Standardize Variable Names

rename_data <- selected_data %>%
  rename(
    mat_age = patient.oocyteAge.years,
    patient_id = patient.esId,
    insem_method = inseminationMethod,
    emb_fate = status,
    gard_exp = clinicGardner.expansion,
    gard_icm = clinicGardner.icm,
    gard_te = clinicGardner.te,
    pgt_status = pgta,
    eq_score = kidScore,
    blast_score = blastAtHr,
    d5_ah = assistedHatching.Day.5..key,
    eq_rank = kidRank,
    d3_ah = assistedHatching.Day3..key,
    d6_ah = manualAnnotation.Day.6.assistedHatching.0..key
  )

# Adjust EQ Score for Beta Regression

n_eq <- sum(!is.na(rename_data$eq_score))

rename_data <- rename_data %>%
  mutate(
    eq_score_beta = ifelse(
      is.na(eq_score),
      NA_real_,
      (eq_score * (n_eq - 1) + 0.5) / n_eq
    )
  )

# Adjust blast score for beta regression

n_blast <- sum(!is.na(rename_data$blast_score))

rename_data <- rename_data %>%
  mutate(
    blast_score_beta = ifelse(
      is.na(blast_score),
      NA_real_,
      (blast_score * (n_blast - 1) + 0.5) / n_blast
    )
  )

# Clean Variable: d3 hatching

rename_data %>%
  filter(!is.na(tB)) %>%
  group_by(d5_ah, d3_ah) %>%
  summarise(
    n = n()
  )

clean_data <- rename_data %>%
  mutate(
    d3_hatching = case_when(
      d3_ah == "ASSISTED_HATCHING_YES" ~ 1L,
      TRUE       ~ 0L
    )
  )

clean_data %>%
  group_by(d3_hatching) %>%
  summarise(
    n = n()
  )

# Clean Variable: ploidy

clean_data %>%
  group_by(pgt_status) %>%
  summarise(
    n = n()
  )

clean_pgt_data <- clean_data %>%
  mutate(
    pgt_euploid = case_when(
      str_detect(pgt_status, "^euploid") ~ 1L,
      is.na(pgt_status)                  ~ NA_integer_,
      TRUE                               ~ 0L
    ),
    pgt_group = case_when(
      str_detect(pgt_status, "^euploid")               ~ "euploid",
      str_detect(pgt_status, "aneuploid|segmental")   ~ "aneuploid",
      str_detect(pgt_status, "^mosaic")                ~ "mosaic",
      str_detect(pgt_status, "^no_result")             ~ "no_result",
      is.na(pgt_status)                                ~ NA_character_,
      TRUE                                             ~ NA_character_
    )
  )

clean_pgt_data %>%
  group_by(pgt_euploid) %>%
  summarise(
    n = n()
  )


# Data Write --------------------------------------------------------------

write.csv(clean_pgt_data, "clean_data.csv", row.names = FALSE)

