### LOAD packages
library(tidyverse)
library(ggsurvfit)
library(survminer)
library(survival)
library(broom)
library(gtsummary)

### LOAD data set

# from ZfKd # Skript von Mike Klora -> Statistikgruppe
GG <- left_join(Tumor, Patient) %>% 
  filter(substr(Diagnose_ICD10_Code,1,3) == "C34" & Diagnosejahr %in% c(2020, 2021)) %>% 
  mutate(Diagnosedatum = as.Date.character(Diagnosedatum),
         DatumVitalStatus = as.Date.character(DatumVitalStatus),
         sex = as.integer(case_when(Geschlecht == "M" ~ 1, Geschlecht == "W" ~ 2, TRUE ~ 0))) %>% 
  mutate(Jahre_Last_Follow_Up = time_length(interval(Diagnosedatum , DatumVitalStatus), "years"),
         Vitalstatus = as.integer(case_when(Verstorben == "J" ~ 1, TRUE ~ 0 ))) %>% 
  filter(!is.na(Vitalstatus) & sex %in% c(1,2)) %>% 
  mutate(T_pc = case_when(!is.na(T_p_2) ~ T_p_2, TRUE ~ T_c_2)) %>% 
  filter(T_pc %in% c("0", "1", "2", "3" ,"4", "x"))

# Example cox-model
res.cox <- coxph(Surv(Jahre_Last_Follow_Up, Vitalstatus) ~ Diagnosealter + Geschlecht + T_pc + Inzidenzort_BL, data =  GG)
summary(res.cox)


##### PACKAGE "forestmodel"
library('forestmodel')

forestmodel::forest_model(res.cox)

##### PACKAGE "forestploter"
# umständlich

library(forestploter)
library(grid)
library(finalfit)

explanatory = c("Diagnosealter", "Geschlecht", "T_pc", "Inzidenzort_BL")
dependent = "Surv(Jahre_Last_Follow_Up, Vitalstatus)"

table <- 
  GG %>%
  finalfit::finalfit( # random effects possible
    dependent = dependent,
    explanatory = explanatory,
    condense = FALSE) %>%
  dplyr::rename(
    "Faktor" = base::paste0("Dependent: ", dependent), 
    "N (%)"="all") %>% 
  dplyr::rename(Kategorie = base::colnames(.)[2]) %>%
  dplyr::mutate(
    HR..multivariable. = base::as.numeric(HR..multivariable.),
    se = (log(U95.y) - log(HR..multivariable.))/1.96,
    " " = paste(rep(" ", 13), collapse = " "),
    `HR (95% CI)` = ifelse(is.na(se), "",
                           sprintf("%.2f (%.2f - %.2f)",
                                   HR..multivariable., L95.y, U95.y)),
    p.y = round(as.numeric(p.y), digits = 3),
    p = ifelse(is.na(p.y), "", p.y)

  )

glimpse(table)

# Define theme
tm <- forest_theme(base_size = 10,
                   refline_col = "red",
                   arrow_type = "closed",
                   footnote_gp = gpar(col = "blue", cex = 0.6))

p <- 
  forest(
  data = table[,c(1:3, 13:15)],
  est = table$HR..multivariable.,
            lower = table$L95.y,
            upper = table$U95.y,
            sizes = 0.4,
            ci_column = 4,
            ref_line = 1,
            arrow_lab = c("Besseres Überleben", "Schlechteres Überleben"),
            xlim = c(0, 4),
            ticks_at = c(0.5, 1, 2, 3),

            theme = tm)

# Print plot
plot(p)  

##### forplo
library(forplo)

forplo(res.cox)
