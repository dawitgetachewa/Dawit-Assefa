#===========================================================
# Title	:Epidemiological shifts and trends in visceral leishmaniasis incidence, relapse, and mortality in Brazil: 
# Analysis using National Notifiable Diseases Information System (SINAN) system, 2007-2023
#===========================================================

# Load necessary packages--------

library(sf)
library(gganimate)
library(SpatialEpi)
library(INLAjoint)
library(plotly)
library(ggplot2)
library(readxl)
library(janitor)
library(dplyr)
library(readr)
library(gghighlight)
library(spdep)
library(RColorBrewer)
library(tidyverse)
library(tmap)
library(here)
library(Hmisc)  # for label function
library(table1)  # for table1 function
library(cli)
library(ggthemes)
library(lme4)
library(splines)
library(geobr)
library(ggforce)  # For better flow lines
library(viridis)  # For color-blind friendly palettes
library(ggraph)
library(ggrepel)
library(gridExtra)
library(cowplot)
library(forcats)
library(ggplot2)
library(dplyr)
library(jpeg)
library(grid)
library(ggspatial)
library(patchwork) # For combining plots
library(cowplot)
library(png)    # For reading PNG images
library(grid)   # For rendering images as grobs
library(classInt) # For creating natural breaks
library(MASS)
library(broom.mixed) 
library(broom)

# Load your dataset---------------------------------------------------------

combined_vl_data <- read.csv("~/Desktop/from PC/IDDO Work/VL in brazil/scripts/combined_dataset_after_drop.csv") %>%
    clean_names()
# data is available at : https://datasus.saude.gov.br/transferencia-de-arquivos/
#===========================================================

# Population data --------
# source (https://www.ibge.gov.br/estatisticas/sociais/populacao/9109-projecao-da-populacao.html)

Pop2007_2023 <-  read_csv("~/Desktop/from PC/IDDO Work/VL in brazil/scripts/new script/population data 2007_2023.csv") %>%
    clean_names() %>%
    group_by(federative_unit, year_of_notification) %>%
    summarise(total_pop_sum = sum(total_pop, na.rm = TRUE), .groups = 'drop')


#===========================================================
# Data cleaning---------------------------------------------
#===========================================================

# Recoding and renaming variables ----------------------------------------
combined_vl_data <- combined_vl_data %>%
    mutate(
        # Recoding `evolution_of_cases`
        evolution_of_cases = case_when(
            evolution_of_cases == 1 ~ "Cured",
            evolution_of_cases == 2 ~ "Abandonment",
            evolution_of_cases == 3 ~ "Death from VL",
            evolution_of_cases == 4 ~ "Death from other causes",
            evolution_of_cases == 5 ~ "Transfer",
            TRUE ~ NA_character_
        ),
        
        # Recoding `type_of_rx`
        type_of_rx = case_when(
            type_of_rx == 1 ~ "Pentavalent antimony",
            type_of_rx == 2 ~ "Amphotericin B deoxycholate",
            type_of_rx == 3 ~ "Pentamidine",
            type_of_rx == 4 ~ "Amphotericin B liposomal",
            type_of_rx == 5 ~ "Treated but not documented properly",
            type_of_rx == 6 ~ "Treated but not documented properly",
            is.na(type_of_rx) ~ "Missing"
        ),
        
        # Recoding `type_of_cases_at_baseline`
        type_of_cases_at_baseline = case_when(
            type_of_cases_at_baseline == 1 ~ "New cases",
            type_of_cases_at_baseline == 2 ~ "Relapse",
            type_of_cases_at_baseline == 3 ~ "Transferred cases",
            TRUE ~ "Other"
        ),
        
        # Recoding `treatment_failure`
        treatment_failure = case_when(
            treatment_failure == 1 ~ "Amphotericin B deoxycholate",
            treatment_failure == 2 ~ "Amphotericin B liposomal",
            treatment_failure == 3 ~ "Other",
            treatment_failure == 6 ~ "Doesn't apply",
            TRUE ~ NA_character_
        ),
        
        # Recoding `race`
        race = case_when(
            race == 1 ~ "White",
            race == 2 ~ "Black",
            race == 3 ~ "Yellow",
            race == 4 ~ "Brown",
            race == 5 ~ "Indigenous",
            race == 9 ~ "Ignored",
            TRUE ~ NA_character_
        )
    )

# Creating age group -----------------------------------------------------
combined_vl_data$Age_Group <- cut(
  as.numeric(combined_vl_data$age_in_years),
  breaks = c(-Inf, 1, 5, 15, 50, Inf),
  labels = c("<1", "1–<5", "5–<15", "15–50", "50+"),
  right = FALSE
)


table(combined_vl_data$`Age-Group`)

# Cleaning symptom and status variables ----------------------------------
combined_vl_data <- combined_vl_data %>%
    mutate(
        fever = case_when(fever == "1" ~ "Yes", fever == "2" ~ "No", TRUE ~ "Ignored"),
        weakness = case_when(weakness == "1" ~ "Yes", weakness == "2" ~ "No", TRUE ~ "Ignored"),
        edema = case_when(edema == "1" ~ "Yes", edema == "2" ~ "No", TRUE ~ "Ignored"),
        weight_loss = case_when(weight_loss == "1" ~ "Yes", weight_loss == "2" ~ "No", TRUE ~ "Ignored"),
        cough_and_dyspnea = case_when(cough_and_dyspnea == "1" ~ "Yes", cough_and_dyspnea == "2" ~ "No", TRUE ~ "Ignored"),
        pallor = case_when(pallor == "1" ~ "Yes", pallor == "2" ~ "No", TRUE ~ "Ignored"),
        spleen_enlargement = case_when(spleen_enlargement == "1" ~ "Yes", spleen_enlargement == "2" ~ "No", TRUE ~ "Ignored"),
        increase_infection = case_when(increase_infection == "1" ~ "Yes", increase_infection == "2" ~ "No", TRUE ~ "Ignored"),
        hemorrhagic_phenomena = case_when(hemorrhagic_phenomena == "1" ~ "Yes", hemorrhagic_phenomena == "2" ~ "No", TRUE ~ "Ignored"),
        liver_enlargement = case_when(liver_enlargement == "1" ~ "Yes", liver_enlargement == "2" ~ "No", TRUE ~ "Ignored"),
        jaundice = case_when(jaundice == "1" ~ "Yes", jaundice == "2" ~ "No", TRUE ~ "Ignored"),
        other = case_when(other == "1" ~ "Yes", other == "2" ~ "No", TRUE ~ "Ignored"),
        
        # Modify Age_Group labels
        `Age Group` = case_when(
            Age_Group == "0-1" ~ "Age Group <1 Year",
            Age_Group == "1-5" ~ "Age Group 1-5 Years",
            Age_Group == "5-15" ~ "Age Group 5-15 Years",
            Age_Group == "15-50" ~ "Age Group 15-50 Years",
            Age_Group == "50+" ~ "Age Group 50+ Years",
            TRUE ~ as.character(Age_Group)
        ),
        
        # Recode Indigenous and HIV status
        `Indigenous Case` = case_when(indigenous_case == "1" ~ "Yes", indigenous_case == "2" ~ "No", TRUE ~ NA_character_),
        `HIV Status` = case_when(
            hiv == "1" ~ "Positive",
            hiv == "2" ~ "Negative",
            is.na(hiv) ~ "Unknown",      # handle NA values
            TRUE ~ "Unknown"             # catch any other unexpected values
        )
    )

# Recode missing outcome
combined_vl_data <- combined_vl_data %>%
    mutate(evolution_of_cases = ifelse(is.na(evolution_of_cases), "Missing", evolution_of_cases))

# Renaming columns -------------------------------------------------------
combined_vl_data <- combined_vl_data %>%
    rename(
        "Type of cases at baseline" = type_of_cases_at_baseline,
        "Evolution of cases" = evolution_of_cases,
        "Type of treatments" = type_of_rx,
        "Treatment failure" = treatment_failure,
        "Gender" = gender,
        "Race" = race,
        "Age-Group" = `Age Group`
    )

# Additional renaming of symptoms ----------------------------------------
combined_vl_data <- combined_vl_data %>%
    rename(
        "Increase Infection" = increase_infection,
        "Bleeding" = hemorrhagic_phenomena,
        "Cough and Diarrhea" = cough_and_dyspnea
    )

# Addressing empty or NA values ------------------------------------------
combined_vl_data <- combined_vl_data %>%
    mutate(
        `Increase Infection` = ifelse(is.na(`Increase Infection`) | `Increase Infection` == "", "Ignored", `Increase Infection`),
        `Bleeding` = ifelse(is.na(`Bleeding`) | `Bleeding` == "", "Ignored", `Bleeding`),
        `Cough and Diarrhea` = ifelse(is.na(`Cough and Diarrhea`) | `Cough and Diarrhea` == "", "Ignored", `Cough and Diarrhea`)
    )

#create year range------------- 
combined_vl_data<- combined_vl_data %>%
  mutate(
    year_range = cut(
      year_of_notification,
      breaks = c(2006, 2010, 2014, 2018, 2023),   # cut points
      labels = c("2007–2009", "2010–2014", "2015–2018", "2019–2023"),
      right = TRUE   # include upper endpoint
    )
  )
#===========================================================
# converting age as patient reported to year----------------------------------------------
# it is a four digit code descriped as 
# "The composition of the variable obeys the following criterion: 1stdigit:
# Hour
# Day
# Month
# Year
# Ex: 3009 – nine months, 4018 – eighteen years
#===========================================================

combined_vl_data$age_in_years <- with(combined_vl_data, {
    # Extract first digit (unit)
    unit <- as.numeric(substr(age_as_pt_repored, 1, 1))
    
    # Extract remaining digits (value)
    value <- as.numeric(substr(age_as_pt_repored, 2, 4))
    
    # Convert all to age in years based on unit
    ifelse(unit == 1, value / (24 * 365),  # hours to years
           ifelse(unit == 2, value / 365,         # days to years
                  ifelse(unit == 3, value / 12,          # months to years
                         ifelse(unit == 4, value, NA))))        # years
})

#============================================================

# filter confirmed cases only ---------------------------------------
combined_vl_data <- combined_vl_data %>%
    filter( final_calssification == 1)

#===========================================================
# Flow diagram----------------------------------------------
#===========================================================

library(DiagrammeR)

grViz("
digraph flowchart {
  graph [layout = dot, rankdir = TB]

  # Define node styles
  node [shape = rectangle, style = filled, color = black, fillcolor = white, fontname = Arial, fontsize = 12]

  # Nodes
  A [label = '102,220 VL cases from 2007–2017']
  B [label = '49,384 VL cases from 2018–2023']
  C [label = '151,609 VL cases from 2007–2023']
  D [label = '84,613 cases discarded']
  E [label = '5,443 inconclusive']
  F [label = '55,723 confirmed VL cases']
  G [label = '36,098 VL cases among patients age >5 years']
  H [label = '18,220 VL cases in under-five children']

  # Edges
  A -> C
  B -> C
  C -> D
  C -> E
  C -> F
  F -> G
  F -> H
}
")



#===========================================================
# Baseline description of patients characteristics -------------
#===========================================================

# Filter dataset with age  and assign unique ID---------------------------

filtered_age <- combined_vl_data %>%
    filter(!is.na(age_in_years))

# summerize the age------------------------------------------------------

cleaned_data <- filtered_age %>%
    summarise(
        # Age statistics
        mean_age = mean(age_in_years, na.rm = TRUE),
        median_age = median(age_in_years, na.rm = TRUE),
        sd_age = sd(age_in_years, na.rm = TRUE),
        iqr25_age = quantile(age_in_years, 0.25, na.rm = TRUE),
        iqr75_age = quantile(age_in_years, 0.75, na.rm = TRUE),
        min_age = min(age_in_years, na.rm = TRUE),
        max_age = max(age_in_years, na.rm = TRUE))


# View the result

cleaned_data

#median age by year range-----------------------------------
cleaned_data_year <- filtered_age %>%
    group_by(year_range) %>%     # group by your year categories
    summarise(
        mean_age   = mean(age_in_years, na.rm = TRUE),
        median_age = median(age_in_years, na.rm = TRUE),
        sd_age     = sd(age_in_years, na.rm = TRUE),
        iqr25_age  = quantile(age_in_years, 0.25, na.rm = TRUE),
        iqr75_age  = quantile(age_in_years, 0.75, na.rm = TRUE),
        min_age    = min(age_in_years, na.rm = TRUE),
        max_age    = max(age_in_years, na.rm = TRUE),
        n          = n()   # number of patients in each period
    ) %>%
    ungroup()

cleaned_data_year

# Age distribution---------------------------------------------

combined_vl_data$`Age-Group` <- factor(combined_vl_data$`Age-Group`, 
                                       levels = c("Age Group <1 Year", 
                                                  "Age Group 1-5 Years", 
                                                  "Age Group 5-15 Years", 
                                                  "Age Group 15-50 Years", 
                                                  "Age Group 50+ Years"))
age_tab <- table(combined_vl_data$`Age-Group`)

age_df <- data.frame(
  Age_Group = names(age_tab),
  Count = as.vector(age_tab),
  Percent = round(100 * age_tab / sum(age_tab), 1)  # rounded %
)

age_df

#Linear regression for age by year range-------------------------------------------
#Fit linear mixed-effects model with federative_unit as random intercept
model_age_mixed <- lmer(age_in_years ~ year_range + (1 | federative_unit), data = filtered_age)

# Summarise fixed effects with 95% CIs
results <- broom.mixed::tidy(model_age_mixed, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    estimate = round(estimate, 1),
    conf.low = round(conf.low, 1),
    conf.high = round(conf.high, 1),
    report = paste0(
      "Compared to 2007–2010, patients were older by ", estimate, " years [95% CI: ",
      conf.low, " to ", conf.high, "]"
    )
  )

print(results$report)


#Table one-------------------------- --------

# -----------------------
# Handle missing values
# -----------------------
combined_vl_data$Race[is.na(combined_vl_data$Race)] <- "Ignored"
combined_vl_data$`HIV Status`[is.na(combined_vl_data$`HIV Status`)] <- "Unknown"
combined_vl_data$`Indigenous Case`[is.na(combined_vl_data$`Indigenous Case`)] <- "Unknown"
combined_vl_data$Gender[is.na(combined_vl_data$Gender)] <- "Unknown"

# -----------------------
# Helper function
# -----------------------
summarize_tabyl <- function(df, var1, var2) {
  tabyl(df, !!sym(var1), !!sym(var2)) %>%
    adorn_percentages("col") %>%
    adorn_pct_formatting(digits = 1) %>%
    adorn_ns()
}

# -----------------------
# Summaries
# -----------------------
# Age sex distribution---------------------------------------- 

summarize_tabyl(combined_vl_data, "Gender", "Age-Group")

# Age group × Race
summarize_tabyl(combined_vl_data, "Race", "Age-Group")

# Age group × HIV status
summarize_tabyl(combined_vl_data, "HIV Status", "Age-Group")

# Age group × Indigenous case
summarize_tabyl(combined_vl_data, "Indigenous Case", "Age-Group")

# Age group × Type of cases at baseline
summarize_tabyl(combined_vl_data, "Type of cases at baseline", "Age-Group")

# Age group × Final outcome
summarize_tabyl(combined_vl_data, "Age-Group", "Evolution of cases")

# Age group × Treatment
summarize_tabyl(combined_vl_data, "Age-Group", "Type of treatments")

# Type of cases × Final outcome
summarize_tabyl(combined_vl_data, "Evolution of cases", "Type of cases at baseline")

# Treatment type × Final outcome
summarize_tabyl(combined_vl_data, "Evolution of cases", "Type of treatments")

# HIV status × Final outcome
summarize_tabyl(combined_vl_data, "Evolution of cases", "HIV Status")

# Treatment failure × Final outcome
summarize_tabyl(combined_vl_data, "Evolution of cases", "Treatment failure")

# Age group × Diagnostics
summarize_tabyl(combined_vl_data, "Age-Group", "parasitological_dx")
summarize_tabyl(combined_vl_data, "Age-Group", "immunological_dx")
summarize_tabyl(combined_vl_data, "Age-Group", "confirmation_criteria")

# Year × Diagnostics
summarize_tabyl(combined_vl_data, "year_of_notification", "parasitological_dx")
summarize_tabyl(combined_vl_data, "year_of_notification", "immunological_dx")

#clinical symptom at presentation----------------------------------

# Define symptom columns

symptoms <- c(
  "fever", "weakness", "edema", "weight_loss", "Cough and Diarrhea",
  "pallor", "spleen_enlargement", "Increase Infection",
  "Bleeding", "liver_enlargement", "jaundice", "other"
)

# Calculate counts and percentages
symptom_summary <- combined_vl_data %>%
  dplyr::select(all_of(symptoms)) %>%
  summarise(across(everything(), ~sum(. == "Yes", na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "Symptom", values_to = "Count") %>%
  mutate(
    Total = nrow(combined_vl_data),
    Percentage = round((Count / Total) * 100, 1)
  ) %>%
  dplyr::select(Symptom, Count, Percentage)

print(symptom_summary)
#==========================================================
#Table 2 Clinical presentation at baseline-------------------------------
#==========================================================
# List of symptom column names
symptom_cols <- c(
  "fever", "weakness", "edema", "weight_loss", "Cough and Diarrhea",
  "pallor", "spleen_enlargement", "Increase Infection",
  "Bleeding", "liver_enlargement", "jaundice", "other"
)

# Create table
symptom_table <- combined_vl_data %>%
  dplyr::select(Age_Group, all_of(symptom_cols)) %>%
  pivot_longer(
    cols = -Age_Group,
    names_to = "Symptom",
    values_to = "Response"
  ) %>%
  filter(Response == "Yes") %>%
  mutate(
    Symptom = case_when(
      Symptom == "fever" ~ "Fever",
      Symptom == "weakness" ~ "Weakness",
      Symptom == "edema" ~ "Edema",
      Symptom == "weight_loss" ~ "Weight loss",
      Symptom == "Cough and Diarrhea" ~ "Cough & Diarrhea",
      Symptom == "pallor" ~ "Pallor",
      Symptom == "spleen_enlargement" ~ "Spleen enlargement",
      Symptom == "Increase Infection" ~ "Increase Infection",
      Symptom == "Bleeding" ~ "Bleeding",
      Symptom == "liver_enlargement" ~ "Liver enlargement",
      Symptom == "jaundice" ~ "Jaundice",
      Symptom == "other" ~ "Other symptoms",
      TRUE ~ Symptom
    )
  ) %>%
  group_by(Age_Group, Symptom) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Age_Group, values_from = Count, values_fill = 0) %>%
  rowwise() %>%
  mutate(Total = sum(c_across(-Symptom))) %>%
  mutate(across(-c(Symptom, Total), ~ round(.x / Total * 100, 1), .names = "{.col}_pct")) %>%
  ungroup()
# View table
symptom_table

#Table 3---------------------------------------

# -----------------------
# Helper function
# -----------------------
summarize_tabyl <- function(df, var1, var2) {
  tabyl(df, !!sym(var1), !!sym(var2)) %>%
    adorn_percentages("row") %>%
    adorn_pct_formatting(digits = 1) %>%
    adorn_ns()
}

combined_vl_data<- combined_vl_data %>%
  mutate(
         year_range = cut(year_of_notification, 
                          breaks = c(2006, 2010, 2014, 2018, 2023), 
                          labels = c("2007-2009", "2010-2014", "2015-2018", "2019-2023")))
# -----------------------
# Summaries
# -----------------------
# Age sex distribution---------------------------------------- 

summarize_tabyl(combined_vl_data, "Age-Group", "Type of treatments")

# type of treatment × year range
summarize_tabyl(combined_vl_data, "year_range", "Type of treatments")

# type of treatment ×  type of cases at baseline
summarize_tabyl(combined_vl_data,  "Type of cases at baseline", "Type of treatments")

# type of treatment ×  HIV status
summarize_tabyl(combined_vl_data,  "HIV Status", "Type of treatments")

#confirmed HIV cases and type of treatment---------

hiv_positive <- combined_vl_data %>%
  filter(`HIV Status` == "Positive")

summarize_tabyl(hiv_positive,  "HIV Status", "Type of treatments")
summarize_tabyl(hiv_positive,  "HIV Status", "Type of cases at baseline")

#==============================================
#Relpase cases------------------------------------
#==============================================
summarize_tabyl <- function(df, var1, var2) {
  tabyl(df, !!sym(var1), !!sym(var2)) %>%
    adorn_percentages("col") %>%
    adorn_pct_formatting(digits = 1) %>%
    adorn_ns()
}

# filter relapse
relapse_cases <- combined_vl_data%>%
  filter(`Type of cases at baseline` == "Relapse")

summarize_tabyl(relapse_cases,  "federative_unit", "Type of cases at baseline")
summarize_tabyl(relapse_cases,  "year_range", "Type of cases at baseline")

# Multi variable logestic reggression relpase cases --------
combined_vl_data <- combined_vl_data %>%
  mutate(
    relapse_binary = ifelse(`Type of cases at baseline` == "Relapse", 1, 0),
    hiv_status_factor = factor(`HIV Status`, levels = c("Negative", "Positive", "Unknown")),
    year_range_factor = factor(year_range, levels = c("2007-2009", "2010-2014", "2015-2018", "2019-2023"))
  )

# Run logistic regression
relapse_model <- glm(
  relapse_binary ~ year_range_factor + hiv_status_factor,
  data = combined_vl_data,
  family = binomial
)

# Summarize results with odds ratios
relapse_model_summary <- tidy(relapse_model, exponentiate = TRUE, conf.int = TRUE)
relapse_model_summary

# relapse among HIV positive cases--------------

# filter relapse
relapse_HIV_cases <- hiv_positive %>%
  filter(`Type of cases at baseline` == "Relapse")

summarize_tabyl(relapse_HIV_cases,  "year_range", "Type of cases at baseline")

table(relapse_HIV_cases$`Type of cases at baseline`)

#===========================================================
#Treatment failure to PA------------------------------------
#===========================================================
combined_vl_data <- combined_vl_data %>%
  mutate(
    is_key_drug = ifelse(
      `Treatment failure` %in% c("Amphotericin B deoxycholate", "Amphotericin B liposomal", "Other"),
      "Yes", "No"
    )
  )


# Flag treatment failures among Pentavalent antimony cases
pa_failure <- combined_vl_data %>%
  filter(`Type of treatments` == "Pentavalent antimony") %>%
  group_by(year_range) %>%
  summarise(
    pa_cases = n(),
    pa_failure_cases = sum(is_key_drug == "Yes"),
    .groups = "drop"
  ) %>%
  mutate(
    failure_percentage = (pa_failure_cases / pa_cases) * 100,
    label = paste0(round(failure_percentage, 1), "% (", pa_failure_cases, "/", pa_cases, ")")
  )

pa_failure


#===========================================================
#Treatment outcomes and mortality------------------------------
#===========================================================
summarize_tabyl <- function(df, var1, var2) {
  tabyl(df, !!sym(var1), !!sym(var2)) %>%
    adorn_percentages("row") %>%
    adorn_pct_formatting(digits = 1) %>%
    adorn_ns()
}

summarize_tabyl(combined_vl_data,  "Evolution of cases", "Type of cases at baseline")


#filter those missing-------------- 
vl_outcomes <- combined_vl_data %>%
  filter(`Evolution of cases` != "Missing")

# Tbale 4-------------------------------------------
table(vl_outcomes $`Evolution of cases`, useNA="ifany")
table(vl_outcomes $`Age-Group`, useNA="ifany")
table(vl_outcomes $`Type of treatments`, useNA="ifany")
table(vl_outcomes $`HIV Status`, useNA="ifany")
table(vl_outcomes $`Type of cases at baseline`, useNA="ifany")
table(vl_outcomes $year_range, useNA="ifany")

summarize_tabyl(vl_outcomes,  "Age-Group", "Evolution of cases")
summarize_tabyl(vl_outcomes,  "Type of treatments", "Evolution of cases")
summarize_tabyl(vl_outcomes, "year_range", "Evolution of cases")
summarize_tabyl(vl_outcomes,  "Type of cases at baseline", "Evolution of cases")
summarize_tabyl(vl_outcomes,  "HIV Status", "Evolution of cases")


#Table 5-----------------------------------------------
death_cases <- combined_vl_data %>%
  filter(`Evolution of cases` %in% c("Death from VL", "Death from other causes"))

table(death_cases $year_range, useNA="ifany")
table(death_cases $`Age-Group`, useNA="ifany")
table(death_cases $`Type of treatments`, useNA="ifany")
table(death_cases $`HIV Status`, useNA="ifany")
table(death_cases $`Evolution of cases`, useNA="ifany")
table(death_cases $`Type of cases at baseline`, useNA="ifany")

# 

#==========================================================
#Figure 2 Age_gender pyramid------------------------------------------
#==========================================================
combined_vl_data<- combined_vl_data %>%
  filter(final_calssification == 1)

# Create age groups and period bins
combined_vl_data <- combined_vl_data %>%
    filter(!is.na(age_in_years), !is.na(year_of_notification)) %>%
    mutate(
        age_group = cut(age_in_years,
                        breaks = seq(0, 105, 5),
                        right = FALSE,
                        labels = paste(seq(0, 100, 5), seq(4, 104, 5), sep = "-")),
        period = case_when(
            year_of_notification %in% 2007:2009 ~ "2007-2009",
            year_of_notification %in% 2010:2014 ~ "2010-2014",
            year_of_notification %in% 2015:2018 ~ "2015-2018",
            year_of_notification %in% 2019:2023 ~ "2019-2023"
        )
    )

# Age-sex pyramid plot
pyramid_data <- combined_vl_data %>%
    filter(!is.na(Gender), !is.na(age_in_years), !is.na(age_group)) %>%
    group_by(period, age_group, Gender) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(count = ifelse(Gender == "Male", -count, count))

pyramid_plot <- ggplot(pyramid_data, aes(x = age_group, y = count, fill = Gender)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(labels = abs) +
  facet_wrap(~period, ncol = 2) +
  labs(title = "A. Age-sex distribution of cases", x = "", y = "") +
  scale_fill_manual(values = c("Male" = "#56B4E9", "Female" = "darkgreen")) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 14)
  )

pyramid_plot

# Export the final plot
tiff("~/Desktop/pyramid_plot.tiff", 
     width = 6000, height = 5000, res = 500,
     compression = "lzw", type = "cairo")
print(pyramid_plot)
dev.off()
#==========================================================
# Figure 2 Median Age by Year (Boxplot)--------------------------------
#==========================================================

filtered_age <- combined_vl_data %>%
    filter(!is.na(age_in_years), !is.na(year_of_notification))

median_age <- ggplot(filtered_age, aes(x = factor(year_of_notification), y = age_in_years)) +
    geom_boxplot(fill = "steelblue", color = "black", outlier.color = "red", outlier.shape = 1) +
    labs(
        title = "B. Age distribution over time",
        x = "",
        y = "Age (years)"
    ) +
    theme_minimal(base_size = 16) +
    theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 14)
    )
median_age 


#==========================================================
#Figure 2 Relapse Proportion by Year and Federative Unit-------------------------
#==========================================================

relapse_cases_per_year <- combined_vl_data %>%
    filter(`Type of cases at baseline` == "Relapse") %>%
    group_by(year_of_notification, federative_unit) %>%
    summarise(relapse_cases = n(), .groups = "drop")

total_cases_per_year <- combined_vl_data %>%
    filter(!is.na(year_of_notification)) %>%
    group_by(year_of_notification, federative_unit) %>%
    summarise(total_cases = n(), .groups = "drop")

relapse_proportion <- left_join(relapse_cases_per_year, total_cases_per_year,
                                by = c("year_of_notification", "federative_unit")) %>%
    mutate(
        percentage = (relapse_cases / total_cases) * 100,
        year_label = paste0(year_of_notification, " (r: ", relapse_cases, " / t: ", total_cases, ")")
    ) %>%
    filter(!is.na(year_of_notification))

overall_trend <- combined_vl_data %>%
    filter(!is.na(year_of_notification)) %>%
    group_by(year_of_notification) %>%
    summarise(
        total_relapse = sum(`Type of cases at baseline` == "Relapse", na.rm = TRUE),
        total_cases = n(),
        percentage = (total_relapse / total_cases) * 100,
        .groups = "drop"
    )

relapse_plot_spaghetti <- ggplot() +
    geom_line(data = relapse_proportion, aes(x = year_of_notification, y = percentage, group = federative_unit),
              color = "skyblue", alpha = 0.5) +
    geom_line(data = overall_trend, aes(x = year_of_notification, y = percentage),
              color = "red", size = 1.5) +
    scale_x_continuous(breaks = seq(min(relapse_proportion$year_of_notification), max(relapse_proportion$year_of_notification), by = 1)) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 75)) +
    labs(
        title = "C. Relapse at baseline, by federative unit",
        x = "",
        y = "Relapse at baseline (%)"
    ) +
    theme_minimal(base_size = 16) +
    theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1)
    )

relapse_plot_spaghetti

 # Export the final plot
tiff("~/Desktop/relapse_plot_spaghetti.tiff", 
     width = 10000, height = 6000, res = 600,
     compression = "lzw", type = "cairo")
print(relapse_plot_spaghetti)
dev.off()
#==========================================================
#Figure 2 Bar graph for VL-HIV co-infection over the year----------------------------
#==========================================================

# Step 1: Prepare data
hiv_status_bar_data <- combined_vl_data %>%
    filter(!is.na(`HIV Status`), !is.na(year_of_notification)) %>%
    group_by(year_of_notification, `HIV Status`) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(year_of_notification) %>%
    mutate(
        total = sum(count),
        Proportion = count / total
    )


# Step 2: Define color-blind-friendly palette
hiv_colors <- c("Negative"= "yellow","Unknown" = "#0072B2", "Positive" = "#D55E00")  # Blue and vermilion

# Step 3: Plot stacked bar chart
hiv_status_bar_plot <- ggplot(hiv_status_bar_data, aes(
    x = factor(year_of_notification),
    y = Proportion,
    fill = `HIV Status`
)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = hiv_colors) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
        title = "D. VL-HIV co-infection",
        x = "",
        y = "VL-HIV co-infection cases",
        fill = "HIV Status"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 12)
    )

# Display the plot
hiv_status_bar_plot

# Export the final plot
tiff("~/Desktop/rhiv_status_bar_plot", 
     width = 10000, height = 6000, res = 600,
     compression = "lzw", type = "cairo")
print(hiv_status_bar_plot)
dev.off()

# Combine plots-------------------------------
# Combine plots: Pyramid on the first column, the rest stacked in the second
final_combined_plot <- plot_grid(
  # First column
  pyramid_plot + theme(plot.margin = margin(10, 10, 10, 10)),       
  # Second column: stack three plots vertically
  plot_grid(
    median_age + theme(plot.margin = margin(10, 10, 10, 10)),     
    relapse_plot_spaghetti + theme(plot.margin = margin(10, 10, 10, 10)),   
    hiv_status_bar_plot + theme(plot.margin = margin(10, 10, 10, 10)),
    ncol = 1,
    align = "v"
  ),
  ncol = 2,
  rel_widths = c(1, 1),  # Adjust relative width if you want the pyramid bigger
  labels = c("A", "B"),
  label_size = 16
)

# Export the final plot
tiff("~/Desktop/FIGURE_2.tiff", 
     width = 14000, height = 10000, res = 600,
     compression = "lzw", type = "cairo")
print(final_combined_plot)
dev.off()
#==========================================================
# Figure 4 Drug regimens used over time (2007-2023), by age-group-------------------------
#==========================================================

# Filter unknown age
filtered_age <- combined_vl_data %>%
  filter(!is.na(age_in_years))

levels = c(
  "Age <1y",
  "Age 1–<5y",
  "Age 5–<15y",
  "Age 15–50y",
  "Age 50+y",
)

# Recode Age_Group into DISCRETE categories (no overlap)
treatment_by_age <- filtered_age %>%
  mutate(
    Age_Group = case_when(
      age_in_years < 1                      ~ "Age <1y",
      age_in_years >= 1  & age_in_years < 5  ~ "Age 1–<5y",
      age_in_years >= 5  & age_in_years < 15 ~ "Age 5–<15y",
      age_in_years >= 15 & age_in_years < 50 ~ "Age 15–50y",
      age_in_years >= 50                     ~ "Age 50+y",
      TRUE ~ "Unknown"
    ),
    Age_Group = factor(
      Age_Group,
      levels = c(
        "Age <1y",
        "Age 1–<5y",
        "Age 5–<15y",
        "Age 15–50y",
        "Age 50+y"
      )
    )
  )

rx_over_year_age <- treatment_by_age %>%
  filter(!is.na(year_of_notification), !is.na(Age_Group)) %>%
  group_by(year_of_notification, Age_Group, `Type of treatments`) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(year_of_notification, Age_Group) %>%
  mutate(
    total = sum(count),
    Proportion = count / total
  )

custom_colors <- c(
  "Amphotericin B deoxycholate" = "#F1E64B",  # yellow
  "Amphotericin B liposomal"   = "#F4A300",  # orange
  "Missing"                    = "#7EC7F2",  # light blue
  "Pentavalent antimony"       = "#1F78B4",  # dark blue
  "Pentamidine"                = "#1B9E77",  # green
  "Treated but not documented" = "#D95F02"   # red-orange
)

# Define full year range (prevents dropped years)
all_years <- sort(unique(rx_over_year_age$year_of_notification))

rx_year_age <- ggplot(
  rx_over_year_age,
  aes(
    x = factor(year_of_notification, levels = all_years),
    y = Proportion,
    fill = `Type of treatments`
  )
) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = "",
    y = "Proportion of cases",
    fill = "Treatment type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  facet_wrap(~ Age_Group)

rx_year_age





# Export with high resolution over all rx change overtime-----
tiff("~/Desktop/change_RX_over_years.tiff", 
     width = 12000, height = 7000, res = 600,
     compression = "lzw", type = "cairo")
print(rx_year_age)
dev.off()

# supplementary figure trend of relpase among HIV positive-------------------------

relapse_HIV_per_year <- hiv_positive %>%
  filter(`Type of cases at baseline` == "Relapse") %>%
  group_by(year_of_notification, federative_unit) %>%
  summarise(relapse_cases = n(), .groups = "drop")

total_HIV_per_year <- hiv_positive %>%
  filter(!is.na(year_of_notification)) %>%
  group_by(year_of_notification, federative_unit) %>%
  summarise(total_cases = n(), .groups = "drop")

relapse_proportion <- left_join(relapse_HIV_per_year, total_HIV_per_year,
                                by = c("year_of_notification", "federative_unit")) %>%
  mutate(
    percentage = (relapse_cases / total_cases) * 100,
    year_label = paste0(year_of_notification, " (r: ", relapse_cases, " / t: ", total_cases, ")")
  ) %>%
  filter(!is.na(year_of_notification))

overall_trend <- hiv_positive %>%
  filter(!is.na(year_of_notification)) %>%
  group_by(year_of_notification) %>%
  summarise(
    total_relapse = sum(`Type of cases at baseline` == "Relapse", na.rm = TRUE),
    total_cases = n(),
    percentage = (total_relapse / total_cases) * 100,
    .groups = "drop"
  )

relapse_HIV_plot_spaghetti <- ggplot() +
  geom_line(data = relapse_proportion, aes(x = year_of_notification, y = percentage, group = federative_unit),
            color = "skyblue", alpha = 0.5) +
  geom_line(data = overall_trend, aes(x = year_of_notification, y = percentage),
            color = "red", size = 1.5) +
  scale_x_continuous(breaks = seq(min(relapse_proportion$year_of_notification), max(relapse_proportion$year_of_notification), by = 1)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 75)) +
  labs(
    title = "A. Proportion of relapse cases among HIV positive cases over time",
    x = "",
    y = "Relapse proportion (%)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

relapse_HIV_plot_spaghetti

# Export with high resolution over all rx change overtime-----
tiff("~/Desktop/relapse_HIV_plot_spaghetti.tiff", 
     width = 8000, height = 5000, res = 600,
     compression = "lzw", type = "cairo")
print(relapse_HIV_plot_spaghetti)
dev.off()

#==========================================================
# Treatment failure to PA-------------
#===========================================================
# Flag key drug failures
combined_vl_data <- combined_vl_data %>%
  mutate(
    is_key_drug = ifelse(
      `Treatment failure` %in% c("Amphotericin B deoxycholate", "Amphotericin B liposomal", "Other"),
      "Yes", "No"
    )
  )

# Filter Pentavalent Antimony cases
PA <- combined_vl_data %>%
  filter(`Type of treatments` == "Pentavalent antimony") %>%
  filter(!is.na(year_of_notification), !is.na(federative_unit))

# Aggregate by year and federative unit
drug_trend_unit <- PA %>%
  group_by(year_of_notification, federative_unit) %>%
  summarise(
    total_cases = n(),
    key_drug_cases = sum(is_key_drug == "Yes"),
    proportion_key_drug = key_drug_cases / total_cases,
    .groups = "drop"
  )

# Overall trend
overall_trend <- PA %>%
  group_by(year_of_notification) %>%
  summarise(
    total_cases = n(),
    key_drug_cases = sum(is_key_drug == "Yes"),
    proportion_key_drug = key_drug_cases / total_cases,
    .groups = "drop"
  )

# Spaghetti plot
treatment_failure_spaghetti <- ggplot() +
  geom_line(data = drug_trend_unit,
            aes(x = year_of_notification, y = proportion_key_drug, group = federative_unit),
            color = "skyblue", alpha = 0.5) +
  geom_line(data = overall_trend,
            aes(x = year_of_notification, y = proportion_key_drug),
            color = "red", size = 1.5) +
  scale_x_continuous(breaks = sort(unique(PA$year_of_notification))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "B. Proportion of treatment failure (PA) over time",
    x = "Year",
    y = "Proportion of treatment failure (PA)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14)
  )

treatment_failure_spaghetti

# Combine plots
final_combined_plot <- plot_grid(
  relapse_HIV_plot_spaghetti + theme(plot.margin = margin(10, 10, 10, 10)),       # A. Median age
  treatment_failure_spaghetti + theme(plot.margin = margin(10, 10, 10, 10)),     # B. Relapse
  labels = c("A", "B"),
  ncol = 1,
  align = "v",
  label_size = 16
)

# Export the final plot
tiff("~/Desktop/suppl_fig.tiff", 
     width = 10000, height = 7000, res = 600,
     compression = "lzw", type = "cairo")
print(final_combined_plot)
dev.off()

#===========================================================
# ------------------ Mortality proportion -------------------
#===========================================================

# Ensure all cases (not just deaths) are included for denominator
total_cases_per_year_unit <- combined_vl_data %>%
    group_by(year_of_notification, federative_unit) %>%
    summarise(total_cases = n(), .groups = "drop")

# Get death counts
mortality <- combined_vl_data %>%
    mutate(death_outcome = ifelse(!is.na(date_of_death), 1, 0)) %>%
    filter(death_outcome == 1)

mortality_cases_per_year_unit <- mortality %>%
    group_by(year_of_notification, federative_unit) %>%
    summarise(total_death = n(), .groups = "drop")

# Join mortality and total cases
mortality_proportion_unit <- left_join(mortality_cases_per_year_unit,
                                       total_cases_per_year_unit,
                                       by = c("year_of_notification", "federative_unit")) %>%
    filter(!is.na(total_cases)) %>%
    mutate(percentage = (total_death / total_cases) * 100)

# Overall trend
overall_mortality <- mortality %>%
    group_by(year_of_notification) %>%
    summarise(total_death = n(), .groups = "drop")

overall_total <- combined_vl_data %>%
    group_by(year_of_notification) %>%
    summarise(total_cases = n(), .groups = "drop")

overall_trend <- left_join(overall_mortality, overall_total, by = "year_of_notification") %>%
    mutate(percentage = (total_death / total_cases) * 100)

# Plot
scale_x_continuous(breaks = seq(min(mortality_proportion_unit$year_of_notification, na.rm = TRUE),
                                max(mortality_proportion_unit$year_of_notification, na.rm = TRUE), by = 1))

mortality_plot_spaghetti <- ggplot() +
    geom_line(data = mortality_proportion_unit,
              aes(x = year_of_notification, y = percentage, group = federative_unit),
              color = "skyblue", alpha = 0.5) +
    geom_line(data = overall_trend,
              aes(x = year_of_notification, y = percentage),
              color = "red", size = 1.5) +
    scale_x_continuous(
        breaks = seq(min(mortality_proportion_unit$year_of_notification, na.rm = TRUE),
                     max(mortality_proportion_unit$year_of_notification, na.rm = TRUE), by = 1)
    ) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 75)) +
    labs(
        title = "Proportion of deaths (case fatality rate) by federative unit over time",
        x = "Year",
        y = "Case fatality rate (%)"
    ) +
    theme_minimal(base_size = 16) +
    theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)  # Optional for readability
    )
mortality_plot_spaghetti

# Export the final plot with high resolution
tiff("~/Desktop/PLOTS/mortality_plot_spaghetti.tiff", 
     width = 5000, height = 3000, res = 500,
     compression = "lzw", type = "cairo")
print(mortality_plot_spaghetti)

dev.off()

#===========================================================
# Mortality prediction-----------------------------------------
#===========================================================

# Remove rows with NA in relevant variables before fitting the model
mortality_complete <- mortality %>%
    filter(!is.na(death_outcome) & !is.na(age_in_years) & !is.na(federative_unit))

# Fit the GLMM
glmm_model <- glmer(death_outcome ~ bs(age_in_years, df = 4) + (1 | federative_unit), 
                    data = mortality_complete, 
                    family = binomial(link = "logit"))

# Generate new data for prediction
age_seq <- seq(min(mortality_complete$age_in_years, na.rm = TRUE), 
               max(mortality_complete$age_in_years, na.rm = TRUE), 
               length.out = 100)

pred_data <- data.frame(
    age_in_years = age_seq,
    federative_unit = mortality_complete$federative_unit[1]  # placeholder value
)

# Get predicted values on the link (logit) scale with standard errors
pred_link <- predict(glmm_model, newdata = pred_data, type = "link", se.fit = TRUE, re.form = NA)

# Calculate 95% confidence intervals on the link scale
pred_data <- pred_data %>%
    mutate(
        fit_link = pred_link$fit,
        se_link = pred_link$se.fit,
        lower_link = fit_link - 1.96 * se_link,
        upper_link = fit_link + 1.96 * se_link,
        predicted_prob = plogis(fit_link),
        lower_prob = plogis(lower_link),
        upper_prob = plogis(upper_link)
    )

# Plot with confidence interval ribbon
prediction<-ggplot(pred_data, aes(x = age_in_years, y = predicted_prob)) +
    geom_line(color = "blue", size = 1.2) +
    geom_ribbon(aes(ymin = lower_prob, ymax = upper_prob), alpha = 0.2, fill = "blue") +
    labs(title = "Predicted probability of death by age with 95% CI",
         x = "Age in years",
         y = "Predicted probability of death") +
    theme_minimal()

prediction

# save high resolution map-------------------------------------------
tiff("~/Desktop/plots/pridiction.tiff", 
     width = 5000, height = 3000, res = 500,
     compression = "lzw", type = "cairo")
print(prediction)
dev.off()

#=======================================================================
# ==== Median age spaghetti plot ====
#median age by year and federative unit--------------------------------

median_age_by_fu <- combined_vl_data %>%
  filter(!is.na(age_in_years), !is.na(year_of_notification), !is.na(federative_unit)) %>%
  group_by(year_of_notification, federative_unit) %>%
  summarise(median_age = median(age_in_years, na.rm = TRUE), .groups = "drop")

# Compute overall median age per year
overall_median_age <- combined_vl_data %>%
  filter(!is.na(age_in_years), !is.na(year_of_notification)) %>%
  group_by(year_of_notification) %>%
  summarise(median_age = median(age_in_years, na.rm = TRUE), .groups = "drop")

# Spaghetti plot
median_age_spaghetti <- ggplot() +
  geom_line(data = median_age_by_fu,
            aes(x = year_of_notification, y = median_age, group = federative_unit),
            color = "skyblue", alpha = 0.5) +
  geom_line(data = overall_median_age,
            aes(x = year_of_notification, y = median_age),
            color = "red", size = 1.5) +
  scale_x_continuous(breaks = seq(min(median_age_by_fu$year_of_notification),
                                  max(median_age_by_fu$year_of_notification), by = 1)) +
  labs(
    title = "",
    x = "Year of Notification",
    y = "Median Age (years)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

median_age_spaghetti

# Export the final plot
tiff("~/Desktop/median_age_spaghetti .tiff", 
     width = 10000, height = 6000, res = 600,
     compression = "lzw", type = "cairo")
print(median_age_spaghetti )
dev.off()


#===========================================================
#Evolution of treatment outcomes by drug over time (2007-2023)--------------------------------------
#===========================================================

# Step 1: Clean and prepare the dataset
combined_vl_data <- combined_vl_data %>%
    filter(!is.na(year_of_notification)) %>%
    mutate(
        `Evolution of cases` = ifelse(is.na(`Evolution of cases`), "Missing", `Evolution of cases`),
        `Type of treatments` = ifelse(is.na(`Type of treatments`), "Missing", `Type of treatments`)
    )

# Step 2: Count and calculate proportions per year and treatment type
age_year_data <- combined_vl_data %>%
    group_by(year_of_notification, `Type of treatments`, `Evolution of cases`) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(year_of_notification, `Type of treatments`) %>%
    mutate(
        total = sum(count),
        Proportion = count / total
    )

# Step 3: Define custom colors for specific outcomes
# Assign specific colors by outcome
outcome_colors <- c(
    "Cured" = "#009E73",                 # Green
    "Death from VL" = "#D55E00",         # Red
    "Death from other causes" = "#E69F00", # Orange
    "Abandonment" = "#F0E442",           # Yellow
    "Transfer" = "#0072B2",              # Blue
    "Missing" = "#999999"                # Grey
)

# Ensure all levels are included (even if not listed above)
unique_outcomes <- unique(age_year_data$`Evolution of cases`)
missing_levels <- setdiff(unique_outcomes, names(outcome_colors))

# Assign remaining levels to additional colors
if (length(missing_levels) > 0) {
    extra_colors <- c("#56B4E9", "#CC79A7", "#000000")  # Add more if needed
    outcome_colors <- c(outcome_colors, setNames(extra_colors[1:length(missing_levels)], missing_levels))
}

# Step 5: Create faceted stacked bar plot with vertical x-axis text and custom fill
rx_by_outcome_overall <- ggplot(age_year_data, aes(
    x = factor(year_of_notification, levels = all_years),
    y = Proportion,
    fill = `Evolution of cases`
)) +
    geom_bar(stat = 'identity', position = 'stack') +
    facet_wrap(~`Type of treatments`, scales = 'free_y') +
    labs(
        title = "Treatment Outcomes by Drug Regimen",
        x = "",
        y = "Proportion of Cases",
        fill = "Outcome"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold"),  # Fully vertical
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 13, face = "bold"),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = outcome_colors)

# Display the updated plot
rx_by_outcome_overall

tiff("~/Desktop/PLOTS/outcomeBy_RX_overYears.tiff", 
     width = 9000, height = 4000, res = 600,
     compression = "lzw", type = "cairo")
print(rx_by_outcome_overall)
dev.off()

#===========================================================
#Evolution of treatment outcomes by drug over time (2007-2023)--------------------
#Under-five--------------------------------------------------
#===========================================================

UNDERFIVE<- combined_vl_data%>%
    filter(age_in_years<5)

# Replace NA values with "Missing" in both Evolution and Treatment type
UNDERFIVE <- UNDERFIVE %>%
    mutate(
        `Evolution of cases` = ifelse(is.na(`Evolution of cases`), "Missing", `Evolution of cases`),
        `Type of treatments` = ifelse(is.na(`Type of treatments`), "Missing", `Type of treatments`)
    )

# Remove rows with NA in the year_of_notification column
UNDERFIVE <- UNDERFIVE %>% filter(!is.na(year_of_notification))


# Assuming combined_vl_data is your data
# Prepare the data
age_year_data <- UNDERFIVE %>%
    group_by(year_of_notification, `Type of treatments`, `Evolution of cases`) %>%
    summarise(count = n(), .groups = 'drop')

# Create the plot
rx_by_outcome_under5 <- ggplot(age_year_data, aes(x = factor(year_of_notification), y = count, fill = `Evolution of cases`)) +
    geom_bar(stat = 'identity', position = 'stack') +
    facet_wrap(~`Type of treatments` , scales = 'free_y') + 
    labs(
        x = "Year of notification",
        y = "Count of cases",
        fill = "Type of outcome") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_brewer(palette = "Set1")  # Colorblind-friendly palette

rx_by_outcome_under5 

#===========================================================
#Evolution of treatment outcomes by drug over time (2007-2023)
#Above-five--------------------------------------------------
#===========================================================

aboveFIVE<- combined_vl_data%>%
    filter(age_in_years>5)

# Replace NA values with "Missing" in both Evolution and Treatment type
aboveFIVE <-aboveFIVE %>%
    mutate(
        `Evolution of cases` = ifelse(is.na(`Evolution of cases`), "Missing", `Evolution of cases`),
        `Type of treatments` = ifelse(is.na(`Type of treatments`), "Missing", `Type of treatments`)
    )

# Remove rows with NA in the year_of_notification column
aboveFIVE <-aboveFIVE %>% filter(!is.na(year_of_notification))


# Assuming combined_vl_data is your data
# Prepare the data
age_year_data <-aboveFIVE %>%
    group_by(year_of_notification, `Type of treatments`, `Evolution of cases`) %>%
    summarise(count = n(), .groups = 'drop')

# Create the plot
rx_by_outcome_above5 <- ggplot(age_year_data, aes(x = factor(year_of_notification), y = count, fill = `Evolution of cases`)) +
    geom_bar(stat = 'identity', position = 'stack') +
    facet_wrap(~`Type of treatments` , scales = 'free_y') + 
    labs(
        x = "Year of notification",
        y = "Count of cases",
        fill = "Type of outcome") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_brewer(palette = "Set1")  # Colorblind-friendly palette

rx_by_outcome_above5

#===========================================================
# Regression--------------------------------------
# Load dataset------------------------------------ 
#===========================================================
combined_vl_data <- read.csv("~/Desktop/from PC/IDDO Work/VL in brazil/scripts/combined_dataset_after_drop.csv") %>%
    clean_names()

# filter confirmed cases only ---------------------------------------
combined_vl_data <- combined_vl_data %>%
    filter( final_calssification == 1)


# Addressing empty or NA values ------------------------------------------

# Scale age for convergence and interpret per Std. Dev change in age
model.fever <- glmer(
    fever == 1 ~ scale(age_in_years, center = FALSE, scale = TRUE) +
        (1 | federative_unit),
    data = combined_vl_data  %>% filter(fever %in% c(1, 2)),
    family = binomial
)

# Model summary
summary(model.fever)

# Confidence intervals using the profile method
confint(model.fever, method = "profile")
exp(fixef(model.fever))  # Odds ratios
exp(confint(model.fever, method = "profile"))  # OR with 95% CI



# Define your list of symptoms
symptoms <- c(
    "fever", "weakness", "edema", "weight_loss", "cough_and_dyspnea ", 
    "pallor", "spleen_enlargement", "increase_infection", 
    "hemorrhagic_phenomena", "liver_enlargement", "jaundice", "other"
)


# Loop through symptoms
# Loop over each symptom and run the model
for (symptom in symptoms) {
    cat("\n\n==============================\n")
    cat("Symptom:", symptom, "\n")
    cat("==============================\n")
    
    # Create a clean formula dynamically
    formula_text <- paste0(symptom, " == 1 ~ scale(age_in_years, center = FALSE, scale = TRUE) + (1 | federative_unit)")
    
    # Filter data to keep only valid values (1 = Yes, 2 = No)
    filtered_data <- combined_vl_data %>% filter(.data[[symptom]] %in% c(1, 2))
    
    # Fit the model
    model <- glmer(
        formula = as.formula(formula_text),
        data = filtered_data,
        family = binomial
    )
    
    # Output model summary
    print(summary(model))
    
    # Confidence intervals
    cat("\n95% Confidence Intervals:\n")
    print(confint(model, method = "profile"))
    
    # Odds Ratios
    cat("\nOdds Ratios:\n")
    print(exp(fixef(model)))
    
    # OR with CI
    cat("\nOR with 95% CI:\n")
    print(exp(confint(model, method = "profile")))
}


#===========================================================
#possion for incidence of cases over time-----------------------------------
#===========================================================
case_summary <- combined_vl_data %>%
    filter(!is.na(federative_unit) & !is.na(year_of_notification)) %>%
    group_by(federative_unit, year_of_notification) %>%
    dplyr::summarize(total_cases = n(), .groups = "drop")

data <- case_summary %>%
    left_join(Pop2007_2023, by = c("federative_unit", "year_of_notification"))

# poision----------
model.2 <- glm(total_cases ~ year_of_notification + federative_unit, offset = log(total_pop_sum), family = poisson(link = "log"), data = data)
summary(model.2)


## Check dispersion parameter with quasi-Poisson regression
model.2q <- glm(total_cases ~ year_of_notification + federative_unit, offset = log(total_pop_sum), family = quasipoisson(link = "log"), data = data)
summary(model.2q)

#why negative binomial?------------------------------
#Explicitly models overdispersion with a variance function:
    Var(Y) = mu+ mu^2/theta
#theta = 0.59, showing substantial overdispersion, but now modeled directly.
#Residual deviance = 487.8 on 399 df → much closer to degrees of freedom compared to Poisson (59,652).
#AIC = 4656.9, dramatically lower than Poisson (61,945) → much better fit.

#negative binomial due to over dispersion
model.3nb <- glm.nb(total_cases ~ year_of_notification + federative_unit + offset(log(total_pop_sum)), data = data)
summary(model.3nb)
confint(model.3nb)
exp(coef(model.3nb))

exp(confint.default(model.3nb))

#===========================================================
# Figure 3 spatio-temporal trends--------------------------------------------
#===========================================================

# Read the dataset and clean column names
combined_vl_data <- read.csv("~/Desktop/from PC/IDDO Work/VL in brazil/scripts/combined_dataset_after_drop.csv") %>%
  clean_names()

combined_vl_data$age_in_years <- with(combined_vl_data, {
  # Extract first digit (unit)
  unit <- as.numeric(substr(age_as_pt_repored, 1, 1))
  
  # Extract remaining digits (value)
  value <- as.numeric(substr(age_as_pt_repored, 2, 4))
  
  # Convert all to age in years based on unit
  ifelse(unit == 1, value / (24 * 365),  # hours to years
         ifelse(unit == 2, value / 365,         # days to years
                ifelse(unit == 3, value / 12,          # months to years
                       ifelse(unit == 4, value, NA))))        # years
})


# Filter confirmed underfive children--------------------------
all_cases<- combined_vl_data %>%
  filter(final_calssification == 1)

# Define the years and required federative unit IDs
years <- 2007:2023
required_ids <- c(12, 13, 15, 16, 17, 21, 22, 23, 24, 25, 26, 27, 28, 29, 
                  31, 32, 33, 35, 41, 42, 43, 50, 51, 52, 53, 11, 14)

# Calculate new cases per federative unit and year
# Ensure both datasets have the same type for 'federative_unit'

# Calculate new cases per federative unit and year
all_cases <- all_cases %>%
  group_by(federative_unit, year_of_notification) %>%
  summarise(total_cases = sum(final_calssification, na.rm = TRUE), .groups = 'drop') %>%
  mutate(federative_unit = as.character(federative_unit)) %>%
  ungroup() %>%
  complete(
    year_of_notification = years, 
    federative_unit = as.character(required_ids),  # Ensure federative_unit is character
    fill = list(total_cases = 0)
  )

# Population data merging------------------------------------------------

Pop2007_2023 <-  read_csv("~/Desktop/from PC/IDDO Work/VL in brazil/scripts/new script/population data 2007_2023.csv") %>%
  clean_names()

# Ensure 'federative_unit' is a character in Pop2007_2022
Pop2007_2023 <- Pop2007_2023 %>%
  mutate(federative_unit = as.character(federative_unit))


# Join the data
all_cases <-all_cases %>%
  left_join(Pop2007_2023, by = c("federative_unit", "year_of_notification"))


# Calculate the incidence of VL per 10,000 population

# Ensure year_range is properly created
all_cases<- all_cases %>%
  mutate(incidence = (total_cases / total_pop) * 100000,
         year_range = cut(year_of_notification, 
                          breaks = c(2006, 2010, 2014, 2018, 2023), 
                          labels = c("2007-2009", "2010-2014", "2015-2018", "2019-2023")))


# Group by federative_unit and year_range to calculate average annual incidence
all_cases$incidence <- as.numeric(all_cases$incidence)

# average annual incidence 
average_annual_incidence <-all_cases%>%
  dplyr::group_by(federative_unit, year_range) %>%
  dplyr::summarize(average_annual_incidence = mean(incidence, na.rm = TRUE), 
                   .groups = "drop")

# chacking for NA
table(average_annual_incidence$year_range, useNA = "ifany")

# Filter out the NA
average_annual_incidence <- average_annual_incidence %>%
  filter(!is.na(year_range))

# Load state data and prepare for joining
states <- read_state(year = 2020, showProgress = FALSE)

states <-states%>%
    rename(federative_unit = code_state) %>%
  mutate(federative_unit = as.character(federative_unit))

# states<- sf::st_layers("~/Desktop/from PC/IDDO Work/VL in brazil/scripts/new script/geoBoundaries-BRA-ADM1-all/geoBoundaries-BRA-ADM1_simplified.shp")
# 
# states<- states%>%
#     rename(federative_unit = shape) %>%
    # mutate(federative_unit = as.character(federative_unit))

# Join state data with new cases
all_case_abovefive <- average_annual_incidence %>%
  left_join(states, by = "federative_unit")


# Save the data as a shapefile for Geoda--------------
st_write(all_case_abovefive, path.expand("~/Downloads/data_shape_underfive.shp"))


# Filter and join data by year ranges
d1 <- states %>%
  left_join(average_annual_incidence %>% filter(year_range == "2007-2009"), by = "federative_unit") %>%
  st_as_sf()

d2 <- states %>%
  left_join(average_annual_incidence %>% filter(year_range == "2010-2014"), by = "federative_unit") %>%
  st_as_sf()

d3 <- states %>%
  left_join(average_annual_incidence %>% filter(year_range == "2015-2018"), by = "federative_unit") %>%
  st_as_sf()

d4 <- states %>%
  left_join(average_annual_incidence %>% filter(year_range == "2019-2023"), by = "federative_unit") %>%
  st_as_sf()


# Define breaks and class labels
breaks <- c(-Inf, 0, 2.4, 4.4, Inf)  # Four breaks for four categories
class_labels <- c("None", ">0-2.4 Sporadic", "2.4-4.4 Moderate", ">4.4 Intense")

# Apply the 'cut' function to classify the incidence variable
d1$class <- cut(
  d1$average_annual_incidence,
  breaks = breaks,
  include.lowest = TRUE,
  labels = class_labels
)

d2$class <- cut(
  d2$average_annual_incidence,
  breaks = breaks,
  include.lowest = TRUE,
  labels = class_labels
)

d3$class <- cut(
  d3$average_annual_incidence,
  breaks = breaks,
  include.lowest = TRUE,
  labels = class_labels
)

d4$class <- cut(
  d4$average_annual_incidence,
  breaks = breaks,
  include.lowest = TRUE,
  labels = class_labels
)

# Ensure the 'class' variable is a factor for proper ordering in plots
d1$class <- factor(d1$class, levels = class_labels)
d2$class <- factor(d2$class, levels = class_labels)
d3$class <- factor(d3$class, levels = class_labels)
d4$class <- factor(d4$class, levels = class_labels)

# Create map function with state names


create_map_with_state_names <- function(data_shape, title) {
  ggplot(data_shape) + 
    geom_sf(aes(fill = class), color = "white", size = 0.2) +
    geom_sf_text(aes(label = name_state), size = 3.5,check_overlap = TRUE,
                 color = "black", fontface = "bold") +  # ← change here
    labs(title = title) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray90", linetype = "dotted"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 14),
      plot.caption = element_text(size = 14, face = "italic")
    ) +
    annotation_north_arrow(
      location = "tr",
      which_north = "true",
      style = north_arrow_fancy_orienteering,
      height = unit(2, "cm"),
      width = unit(2, "cm"),
      pad_x = unit(0.4, "cm"),
      pad_y = unit(0.4, "cm")
    ) +
    annotation_scale(
      location = "bl",
      width_hint = 0.45,
      bar_cols = c("gray60", "white"),
      line_width = 0.8,
      text_cex = 0.8
    )
}

# Create and display maps for each year range
# Replace d1–d4 with your real data (sf objects with `class` and `NM_UF` columns)
plot_2007_2009_with_names <- create_map_with_state_names(d1, "2007–2009")
plot_2010_2014_with_names <- create_map_with_state_names(d2, "2010–2014")
plot_2015_2018_with_names <- create_map_with_state_names(d3, "2015–2018")
plot_2019_2023_with_names <- create_map_with_state_names(d4, "2019–2023")


# Step 1: Define the shared color scale
custom_fill <- scale_fill_manual(
  name = "Incidence Classification",
  values = c(
    "None" = "#009E73",
    ">0-2.4 Sporadic" = "#56B4E9",
    "2.4-4.4 Moderate" = "#E69F00",
    ">4.4 Intense" = "#D55E00"
  )
)

# Step 2: Add the scale to each plot and remove legends
plot_2007_2009_with_names <- plot_2007_2009_with_names + custom_fill + theme(legend.position = "none")
plot_2010_2014_with_names <- plot_2010_2014_with_names + custom_fill + theme(legend.position = "none")
plot_2015_2018_with_names <- plot_2015_2018_with_names + custom_fill + theme(legend.position = "none")
plot_2019_2023_with_names <- plot_2019_2023_with_names + custom_fill + theme(legend.position = "none")

# Step 3: Extract legend from one plot
legend_plot <- plot_2007_2009_with_names + theme(legend.position = "right")
legend <- get_legend(legend_plot)

# Step 4: Combine plots without legends
combined_plots <- plot_grid(
  plot_2007_2009_with_names,
  plot_2010_2014_with_names,
  plot_2015_2018_with_names,
  plot_2019_2023_with_names,
  labels = c("A", "B", "C", "D"),
  ncol = 2
)

# Step 5: Add the legend to the right
combined_plot_with_legend <- plot_grid(
  combined_plots,
  legend,
  rel_widths = c(1, 0.2)
)

# save high resolution map-------------------------------------------
tiff("~/Desktop/MAP_all_age.tiff", 
     width = 14000, height = 12000, res = 700,
     compression = "lzw", type = "cairo")
print(combined_plot_with_legend)
dev.off()

#=========================================================
# Regression for all cause mortality----------------------
#=========================================================

# Convert date_of_death to binary (1 = death, 0 = no death)
death_data <- combined_vl_data %>%
  mutate(death_outcome = ifelse(!is.na(date_of_death), 1, 0))

# Factor coding with reference groups
death_data$`HIV Status` <- relevel(as.factor(death_data$`HIV Status`), ref = "Negative")
death_data$`Age-Group` <- relevel(as.factor(death_data$`Age-Group`), ref = "Age Group 5-15 Years")
death_data$`Type of cases at baseline` <- relevel(as.factor(death_data$`Type of cases at baseline`), ref = "New cases")
death_data$federative_unit <- as.factor(death_data$federative_unit)
death_data$year_range <- relevel(as.factor(death_data$year_range), ref = "2007-2009")
death_data$`Type of treatments` <- relevel(as.factor(death_data$`Type of treatments`), ref = "Pentavalent antimony")

# Variables for univariable analysis
predictors <- c("HIV Status", "Age-Group", "Type of cases at baseline", 
                "year_range", "Type of treatments")

# Function to fit univariable mixed logistic regression and extract OR + CI
univariable_results <- lapply(predictors, function(var) {
  # Build formula safely with backticks if needed
  formula <- as.formula(paste("death_outcome ~", sprintf("`%s`", var), "+ (1 | federative_unit)"))
  
  model <- glmer(formula, data = death_data, family = binomial)
  
  # Extract ORs and CIs
  or <- exp(fixef(model))
  ci <- exp(confint(model, parm = "beta_", method = "Wald"))
  
  data.frame(
    Variable = names(or),
    OR = or,
    Lower_CI = ci[,1],
    Upper_CI = ci[,2],
    row.names = NULL
  )
})

# Combine into one table
univariable_table <- do.call(rbind, univariable_results)

# Display results
print(univariable_table)

#=========================================================
# multivariable model--------------------------------------
#=========================================================

# Multivariable mixed logistic regression
multi_model <- glmer(
  death_outcome ~ `HIV Status` + `Age-Group` + `Type of cases at baseline` +
    year_range + `Type of treatments` + (1 | federative_unit),
  data = death_data,
  family = binomial
)

# Model summary
summary(multi_model)

# Odds ratios (fixed effects only)
or <- exp(fixef(multi_model))

# 95% Wald confidence intervals for fixed effects
ci <- exp(confint(multi_model, parm = "beta_", method = "Wald"))

# Combine into a results table
multi_results <- data.frame(
  Variable = names(or),
  OR = or,
  Lower_CI = ci[,1],
  Upper_CI = ci[,2],
  row.names = NULL
)

print(multi_results)
