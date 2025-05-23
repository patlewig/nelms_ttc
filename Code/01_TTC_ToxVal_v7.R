########################################################################################
#
# Author: Mark D. Nelms, Ph.D., nelms.mark@epa.gov
#
# Version: 1.0 12th September 2018
#
# Description: Extracts SMILES from DSSTox that are also present in ToxVal. Once run through
#             Toxtree chems are then split into their corresponding TTC category. 
#             2nd part - IDs study info for each Cramer class for ToxVal. 
#             3rd part - Generates boxplots & removes extreme outliers (1.5x IQR)
#             4th part - Compare distributions using Kolmorogov-Smirnov test
#             5th part - Calculate 5th %ile values for each Cramer class (ToxVal)
#             6th part - Recreate Munro & use K-S to compare ToxVal & Munro distributions
#             7th part - Calculate if 5th %ile values are different for each Cramer class for ToxVal & Munro
#             8th part - Calculate confidence intervals around 5th %ile values using bootstrapping
#             9th part - Begins investigating discrepancy in 5th %ile values for Class III chemcials
#                         - this is done using Ryan's chemotype enrichment workflow
#
# Notes:
#
#
# Potential Issues: None known
#########################################################################################
library(DBI)
library(here)
library(scales)
library(caret)
library(janitor)
library(fitdistrplus)
library(ggthemes)
library(WRS2)
library(tidyverse)

Input <- "InputData"
Inter <- "IntermediateData"
Output <- "OutputFiles"


###   Extracts SMILES from MySQL database & associates it w/DTXSID from ToxVal    ###
###   Writes SMILES & DTXSID out to be run through Toxtree    ###


## This contains info to log in to MySQL (including username and password)
mycnf <- here(Input, "my.cnf")

## ID which mySQL database to connect to
dscon <- dbConnect(drv = RMySQL::MySQL(), default.file = mycnf, group = "toxval")

## Gather the GSID, DTXSID, CAS, and SMILES from the
ToxVal <- dbGetQuery(dscon, "select * from dev_toxval_v7.toxval
                     join dev_toxval_v7.chemical using (chemical_id)")  # Joins QSAR ready SMILES to ToxVal w/out duplicate chemical_id col

## Don't need all columns as there are duplicate
col.to.drop <- c("chemical_id", "toxval_id", "species_id", "study_id", "priority_id","toxval_uuid",
                "parent_toxval_uuid", "manual_qc_flag", "source_url", "subsource_url", "record_url", "source_source_id", 
                "toxval_numeric_converted", "toxval_units_converted", "chemical_id_external")

## Read in ToxVal csv file - rename dsstox_substance_id column - then convert to data table
## Can't use "fread" - get error when using "rename"
ToxVal <- ToxVal %>%
  rename(., DTXSID = dsstox_substance_id) %>%
  rename(species_common = species_original) %>%
  select(-col.to.drop) %>%
  select(-ends_with("_original")) %>%
  select(DTXSID, name, casrn, qsar_ready_smiles, everything())
## Write out ToxVal
write_csv(ToxVal, here(Input, "ToxVal_v7.csv"))

## Keep only the smiles and DTXSID columns
SMILES <- ToxVal %>%
  distinct(DTXSID, .keep_all = TRUE) %>%
  select(qsar_ready_smiles, DTXSID, name) %>% 
  ## Remove "|^1:0|" string present after some of the SMILES
  mutate(qsar_ready_smiles = str_replace(qsar_ready_smiles, "\\|\\^1:0\\|", "")) %>% 
  ## Remove rows that contain NA in qsar_ready_smiles col
  ## Remove rows w/ "|" in smiles col - they are Markush struc not recognised when converted to SDF
  filter(!is.na(qsar_ready_smiles) | qsar_ready_smiles == "NULL") %>%
  filter(!str_detect(qsar_ready_smiles, "\\|"))

## Write the SMILES variable as a tab separated file
write.table(SMILES, here(Inter, "ToxVal_v7_QSAR_SMILES.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)


####                                                 ####                                                     ####
##    Pass tsv through KNIME to generate SDF - Think ChemmineR generates V3000 SDF - not compatible w/Toxtree   ##
##    Pass SDF through Toxtree to run through Cramer, Kroes, Carbamate, OP, and Steroid decision trees          ##
####                                                 ####                                                     ####

## Read in results from Cramer workflow
Cramer <- read_csv(here(Input, "Cramer_results.csv"), col_types = cols()) %>% 
  clean_names() %>% 
  rename(DTXSID = dtxsid,
         qsar_ready_smiles = smiles) %>%
  filter(!grepl("DTXSID", cramer_rules)) %>% ## Remove chemicals that didn't go through workflow properly
  select(DTXSID, qsar_ready_smiles, name, cramer_rules) ## Keep only these 3 cols in this order

## Read in results from Kroes workflow
Kroes <- read_csv(here(Input, "Kroes_results.csv"), col_types = cols()) %>%
  clean_names() %>% 
  rename(DTXSID = dtxsid,
         kroes_rules = kroes_ttc_decision_tree, 
         kroes_explanation = contains("explanation")) %>%
  select(., c(DTXSID, kroes_rules, kroes_explanation)) 

## Read in results from Carbamate workflow
Carbamate <- read_csv(here(Input, "Carbamate_results.csv"), col_types = cols()) %>%
  clean_names() %>% 
  rename(DTXSID = dtxsid,
         carbamate_rules = alkyl_carbamate_and_thiocarbamate) %>%
  select(., c(DTXSID, carbamate_rules)) 

## Read in results from OP workflow
OP <- read_csv(here(Input, "OP_results.csv"), col_types = cols()) %>%
  clean_names() %>% 
  rename(DTXSID = dtxsid, 
         OP_rules = o_ps, 
         OP_explanation = contains("explanation")) %>%
  select(., c(DTXSID, OP_rules, OP_explanation)) 

## Read in results from Steroid workflow
Steroid <- read_csv(here(Input, "Steroid_results.csv"), col_types = cols()) %>%
  clean_names() %>% 
  rename(DTXSID = dtxsid,
         steroid_rules = steroids) %>%
  select(., c(DTXSID, steroid_rules)) 

## Merge Cramer w/other data by DTXSID
TTC <- Cramer %>%
  inner_join(., Kroes, by = "DTXSID") %>%
  inner_join(., Carbamate, by = "DTXSID") %>%
  inner_join(., OP, by = "DTXSID") %>%
  inner_join(., Steroid, by = "DTXSID")

## Filtering those excluded from TTC
Excluded <- TTC %>%
  filter(., kroes_rules == "Risk assessment requires compound-specific toxicity data") %>%
  select(., DTXSID, name, qsar_ready_smiles, kroes_rules, kroes_explanation)
write_csv(Excluded, here(Inter, "Excluded_from_TTC_211218.csv"))

## Filtering for High potency carcinogens
HighPotCarc <- TTC %>%
  filter(., kroes_rules == "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that contain "Risk ..."
  filter(., str_detect(kroes_explanation, pattern = "Q2Y") & str_detect(kroes_explanation, pattern = "Q3Y")) %>% ## Filter for rows answering Yes to Q2 & Q3
  select(., DTXSID, name, qsar_ready_smiles, kroes_rules, kroes_explanation)
write_csv(HighPotCarc, here(Inter, "HighPotCarc_211218.csv"))

## Filtering for Steroids
Ster <- TTC %>% 
  filter(., str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that contain "Negligible"
  filter(., str_detect(steroid_rules, "Default Class 1")) %>% ## Filter for rows designated as Carbamates
  select(., DTXSID, name, qsar_ready_smiles, kroes_rules, steroid_rules)
write_csv(Ster, here(Inter, "TTC_Steroids_211218.csv"))

## Filtering for Organophosphates
Organo <- TTC %>%
  filter(., kroes_rules == "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that contain "Risk ..."
  filter(., str_detect(kroes_explanation, pattern = "Q1Y")) %>% ## Filter for rows answering Yes to Q1
  filter(., OP_rules == "Default Class 1") %>% ## Filter for rows designated as OPs
  select(., DTXSID, name, qsar_ready_smiles, kroes_rules, kroes_explanation, OP_rules)
write_csv(Organo, here(Inter, "TTC_OPs_211218.csv"))

## Filtering for Carbamates
Carb <- TTC %>%
  filter(., str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that contain "Negligible"
  filter(., str_detect(carbamate_rules, "Default Class 1")) %>% ## Filter for rows designated as Carbamates
  select(., DTXSID, name, qsar_ready_smiles, kroes_rules, carbamate_rules)
write_csv(Carb, here(Inter, "TTC_Carbamates_211218.csv"))

## Combining OPs and Carbamates - i.e. anticholinesterase inhibitors
Anti.CHE <- full_join(Carb, Organo)
write_csv(Anti.CHE, here(Inter, "Anti_CHE.csv"))

## Filtering for chems that:
## 1) Require compound-specific tox data & answered Yes to Q1,
## 2) are NOT OPs,
## 3) ARE steroids
NA.for.TTC <- TTC %>%
  filter(., kroes_rules == "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that contain "Risk ..."
  filter(., str_detect(kroes_explanation, pattern = "Q1Y")) %>% ## Filter for rows answering Yes to Q1
  filter(., str_detect(OP_rules, "Default Class 2")) %>% ## Filter for rows that are NOT designated as OPs
  full_join(., Ster, by = "DTXSID") %>% ## Join with chemicals ID'd as Steroids
  full_join(., HighPotCarc, by = "DTXSID") ## Join with chemicals pulled out as High Potency Carcinogens
write_csv(NA.for.TTC, here(Inter, "NA_for_TTC_211218.csv"))

## Filtering for Genetox that requires a TTC of 0.15ug/day
GenTox.TTC <- TTC %>%
  filter(., str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that contain "Negligible"
  filter(., str_detect(steroid_rules, "Default Class 2")) %>% ## Filter for rows that are NOT steroids
  filter(., str_detect(carbamate_rules, "Default Class 2")) ## Filter for rows that are NOT carbamates
write_csv(GenTox.TTC, here(Inter, "GenTox_TTC_211218.csv"))

## Filtering for Cramer Class I chemicals
Class.I <- TTC %>%
  filter(., kroes_rules != "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that do NOT contain "Risk ..."
  filter(., !str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that do NOT contain "Negligible" 
  filter(., str_detect(cramer_rules, "Low")) ## Filter for rows that are in Cramer Class I
write_csv(Class.I, here(Inter, "Class_I_211218.csv"))

## Filtering for Cramer Class II chemicals
Class.II <- TTC %>%
  filter(., kroes_rules != "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that do NOT contain "Risk ..."
  filter(., !str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that do NOT contain "Negligible"
  filter(., str_detect(cramer_rules, "Intermediate")) ## Filter for rows that are in Cramer Class II
write_csv(Class.II, here(Inter, "Class_II_211218.csv"))

## Filtering for Cramer Class III chemicals
Class.III <- TTC %>%
  filter(., kroes_rules != "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that do NOT contain "Risk ..."
  filter(., !str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that do NOT contain "Negligible"
  filter(., str_detect(cramer_rules, "High")) ## Filter for rows that are in Cramer Class III
write_csv(Class.III, here(Inter, "Class_III_211218.csv"))


######################                  ######################                  ######################                  ######################
###   ID chemicals in ToxVal present in each Cramer class & extract study data    ###
###   These can then be used to calculate 5th %ile NO(A)EL & TTC thresholds    ###


## Read in ToxVal data and TTC group info
ToxVal <- read_csv(here(Input, "ToxVal_v7.csv"),
                   col_types = cols(
                     study_duration_value = col_double())) %>% # Need to specify double - will default to INT
  mutate(., species_common = as.factor(tolower(species_common)),
         study_duration_class = as.factor(tolower(study_duration_class)),
         study_type = as.factor(tolower(study_type)),
         toxval_type = as.factor(tolower(toxval_type))) # change these columns to lowercase
  
Class.I <- read_csv(here(Inter, "Class_I_211218.csv"), col_types = cols()) 
Class.II <- read_csv(here(Inter, "Class_II_211218.csv"), col_types = cols()) 
Class.III <- read_csv(here(Inter, "Class_III_211218.csv"), col_types = cols())
#Gentox <- read_csv(here(Inter, "GenTox_TTC_211218.csv"))
#ACH <- read_csv(here(Inter, "Anti_CHE_211218.csv")) 
#NotApp <- read_csv(here(Inter, "NA_for_TTC_211218.csv")) 

## Set the study parameters of interest here - only have to change info here if want to narrow or broaden scope
## Data source
#study.source <- "ToxRefDB"
## Type of study (sub)chronic, developmental, reproductive
study.duration <- c("subchronic", "chronic", "reproductive", "developmental", "multigeneration")
## Type of data wanting
data.type <- c("noael", "noel")
## Exposure route
route <- "oral"
## Animal - we want rodent species & human. Any other species we do/don't want?
## Exact species names - if use partial matching these bring in many unwanted species
species <- c("albino rat", "albino laboratory rat", "mus booduga", "mus domesticus",
             "mus musculus", "mus platythrix", "mus sp.", "rat", "rat, dog", "rat & dog", "rat, young")
## Partial species names - these won't bring in unwanted species
species.par <- c("cavia", "gerbil", "guinea pig", "hamster", "mesocricetus", "microtus", 
                 "mouse", "oryctolagus", "peromyscus", "rabbit", "rats", "rattus", "rodent",
                 "sylvilagus", "vole")
## Measurement units
ttc.unit <- "mg/kg-day"
## Set up so can divide subchronic NO(A)ELs by 3
subchr <- c("subchronic")

## Function IDs chemicals in a specific Cramer class
## Then extracts ToxVal data that follow study parameters outlined above
id_cramer_in_toxval <- function(cramer, df = ToxVal) {
  ## This code filters ToxVal by chems in class of interest (ID'ed above)
  ## Then filters by specific study info (source, duration, data type, species, and route of exposure)
  ## Where study duration == subchronic divides toxval number by 3 (safety factor from subchronic to chronic)
  #
  ## Args:
  ## cramer: Cramer class csv to use to filter ToxVal
  ## df: dataset to use - default is ToxVal
  c.class <- gsub("\\.", " ", deparse(substitute(cramer))) ## Takes the name of cramer & swaps the "." for a space
  
  x <- df %>%
    filter(., DTXSID %in% cramer$DTXSID) %>%
    filter(., #source %in% study.source &
           toxval_type %in% data.type &
             exposure_route %in% route &
             toxval_units %in% ttc.unit &
             toxval_numeric > 0) %>% # Can't have negative toxicity values!
    # Filter species using either exact or partial names
    # str_detect used so don't have to list out every single variant of animal name
    filter(., species_common %in% species | str_detect(species_common, 
                                                       pattern = paste(species.par, collapse = "|"))) %>%
    # Filter study type for those study types of interest & w/ length containing "chronic" or a "-"
    filter(., str_detect(study_type, pattern = paste(study.duration, collapse = "|")) &
             (str_detect(study_duration_class, pattern = "chronic") | study_duration_class == "-")) %>%
    mutate(toxval_numeric = ifelse(study_type %in% subchr, toxval_numeric / 3, toxval_numeric)) %>%
    mutate(log_noel = log10(toxval_numeric)) %>% 
    mutate(structural_class = c.class) %>% ## Use the variable generated above to add Cramer class column
    select(structural_class, DTXSID, name, casrn, qsar_ready_smiles, toxval_numeric, log_noel, everything())
}


### CLASS I ###
Class.I.data <- id_cramer_in_toxval(Class.I)

### CLASS II ###
Class.II.data <- id_cramer_in_toxval(Class.II)

### CLASS III ###
Class.III.data <- id_cramer_in_toxval(Class.III)


######################                  ######################                  ######################                  ######################
###   Generate boxplots to show distribution of NO(A)ELs across chemicals   ###
###   w/ multiple NO(A)EL values for each Cramer Classification   ###


## Colours for each Cramer class - for use in plots
## Defining here enables easy changing of colours if desired
class.colours = c("#009988", "#EE7733", "#CC3311")  # Teal, Orange, Red

## Cramer Class I
Class.I.data %>%
  group_by(DTXSID) %>% 
  arrange(DTXSID) %>% 
  filter(n() > 1) %>% 
  ggplot(., aes(x = DTXSID, y = toxval_numeric)) +
  geom_boxplot(fill = class.colours[1]) +
  coord_trans(y = "log10") +
  scale_y_continuous(
    breaks = scales::log_breaks(n = 6, base = 10), # Change n to change # labels shown
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  xlab("Chemical DSSTox ID") +
  ylab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here(Output, "Class_I_NOAEL_distribution.png"), device = "png", 
       width = 9, height = 6, units = "in", dpi = 600, bg = "transparent")

## Cramer Class II
Class.II.data %>% 
  group_by(DTXSID) %>% 
  arrange(DTXSID) %>% 
  filter(n() > 1) %>% 
  ggplot(., aes(x = DTXSID, y = toxval_numeric)) +
  geom_boxplot(fill = class.colours[2]) +
  coord_trans(y = "log10") +
  scale_y_continuous(
    breaks = scales::log_breaks(n = 6, base = 10), # Change n to change # labels shown
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  xlab("Chemical DSSTox ID") +
  ylab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here(Output, "Class_II_NOAEL_distribution.png"), device = "png", 
       width = 9, height = 6, units = "in", dpi = 600, bg = "transparent")


## Cramer Class III
Class.III.data %>% 
  group_by(DTXSID) %>% 
  arrange(DTXSID) %>% 
  filter(n() > 1) %>% 
  ggplot(., aes(x = DTXSID, y = toxval_numeric)) +
  geom_boxplot(fill = class.colours[3]) +
  coord_trans(y = "log10") +
  scale_y_continuous(
    breaks = scales::log_breaks(n = 7, base = 10), # Change n to change # labels shown
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  xlab("Chemical DSSTox ID") +
  ylab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here(Output, "Class_III_NOAEL_distribution.png"), device = "png", 
       width = 9, height = 6, units = "in", dpi = 600, bg = "transparent")


## Removal of data points outside of Tukey fences from boxplots
isnt_out_tukey <- function(x, k = 1.5, na.rm = TRUE) {
  ## ID which values are outside Tukey fence (i.e. 1.5x InterQuartile Range)
  #
  ## Args:
  ## x: data & column to use to calculate IQR for
  ## k: constant to ID Tukey fence - default is 1.5
  ## Returns: Boolean of if value is w/in Tukey fence 
  ## (TRUE = inside, FALSE = outlier)
  quar <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- diff(quar)
  
  tuk <- (quar[1] - k * iqr <= x) & (x <= quar[2] + k * iqr)
  tukey <- tuk
}

## ID values within Tukey fence (i.e. non extreme values)
## Remove outliers & take min NO(A)EL for each chem
## CLASS I
Class.I.min <- Class.I.data %>%
  group_by(DTXSID) %>%
  mutate(tukey = isnt_out_tukey(toxval_numeric)) %>% 
  filter(tukey == TRUE) %>% 
  slice(which.min(toxval_numeric)) # Change to summarise(toxval_numeric = mean|median(toxval_numeric))
#write_csv(Class.I.min, here(Output, "Class_I_min_NOAEL.csv"))

## CLASS II
Class.II.min <- Class.II.data %>%
  group_by(DTXSID) %>%
  mutate(tukey = isnt_out_tukey(toxval_numeric)) %>% 
  filter(tukey == TRUE) %>% 
  slice(which.min(toxval_numeric))
#write_csv(Class.II.min, here(Output, "Class_II_min_NOAEL.csv"))


## CLASS III
Class.III.min <- Class.III.data %>%
  group_by(DTXSID) %>%
  mutate(tukey = isnt_out_tukey(toxval_numeric)) %>% 
  filter(tukey == TRUE) %>% 
  slice(which.min(toxval_numeric))
#write_csv(Class.III.min, here(Output, "Class_III_min_NOAEL.csv"))


## Join all Cramer classes together ready for plotting CDFs below
Class.com <- Class.I.min %>%
  full_join(., Class.II.min) %>%
  full_join(., Class.III.min) %>%
  select(structural_class, DTXSID, name, toxval_numeric, log_noel, casrn)
#write_csv(Class.com, here(Output, "Complete_Cramer_class_min_NOAEL.csv"))


######################                  ######################                  ######################                  ######################
###   Using Kolmorogov-Smirnov test to compare distributions between    ###
###   Class I, II, & III chemicals from ToxVal   ###


###   Kolmorogov-Smirnov tests between ToxVal Cramer classes    ###
## Class I vs Class II
ks.test(Class.I.min$toxval_numeric, Class.II.min$toxval_numeric, alternative = "two.sided")
# D = 0.12149, p-value = 0.6545

## Class I vs Class III
ks.test(Class.I.min$toxval_numeric, Class.III.min$toxval_numeric, alternative = "two.sided")
# D = 0.23628, p-value = 1.332e-15

## Class II vs Class III
ks.test(Class.II.min$toxval_numeric, Class.III.min$toxval_numeric, alternative = "two.sided")
# D = 0.2248, p-value = 0.06615


######################                  ######################                  ######################                  ######################
### Generation of CDF plots for each of the Cramer class  ###

### CLASS I ###

## Calculate log normal distribution for Class I chemicals
CI.lnrm <- fitdistr(Class.I.min$toxval_numeric, "lognormal")

## Generate scatterplot with log normal distribution
ggplot(Class.I.min, aes(x = toxval_numeric)) +
  stat_ecdf(geom = "point", pad = FALSE, colour = class.colours[1]) +  
  # Add log normal line to cumulative distribution
  stat_function(fun = plnorm, args = list(CI.lnrm$estimate[1], CI.lnrm$estimate[2]),
                colour = class.colours[1]) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 6, base = 10), # Change n to change # labels shown
    labels = scales::trans_format("log10", scales::math_format(.x))
    ) +
  annotation_logticks(sides = "b") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = percent_format()) + # This puts into % scale to 0 decimal places
  xlab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  ylab("Cumulative Frequency") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA)
  )


### CLASS II ###

## Calculate log normal distribution for Class II chemicals
CII.lnrm <- fitdistr(Class.II.min$toxval_numeric, "lognormal")

## Generate scatterplot with log normal distribution
ggplot(Class.II.min, aes(x = toxval_numeric)) +
  stat_ecdf(geom = "point", pad = FALSE, colour = class.colours[2]) + 
  # Add log normal line to cumulative distribution
  stat_function(fun = plnorm, args = list(CII.lnrm$estimate[1], CII.lnrm$estimate[2]),
                colour = class.colours[2]) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 6, base = 10),
    labels = scales::trans_format("log10", scales::math_format(.x))
  ) +
  annotation_logticks(sides = "b") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = percent_format()) + # This puts into % scale to 0 decimal places
  xlab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  ylab("Cumulative Frequency") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA)
  )


### CLASS III ###

## Calculate log normal distribution for Class III chemicals
CIII.lnrm <- fitdistr(Class.III.min$toxval_numeric, "lognormal")

## Generate scatterplot with log normal distribution
ggplot(Class.III.min, aes(x = toxval_numeric)) +
  stat_ecdf(geom = "point", pad = FALSE, colour = class.colours[3]) + 
  # Add log normal line to cumulative distribution
  stat_function(fun = plnorm, args = list(CIII.lnrm$estimate[1], CIII.lnrm$estimate[2]),
                colour = class.colours[3]) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 6, base = 10),
    labels = scales::trans_format("log10", scales::math_format(.x))
  ) +
  annotation_logticks(sides = "b") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = percent_format()) + # This puts into % scale to 0 decimal places
  xlab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  ylab("Cumulative Frequency") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA)
  )


###   Cumulative distributions of all Cramer classifications    ###
## Plot all CDFs on same graph
ggplot(Class.com, aes(x = toxval_numeric, colour = structural_class)) +
  stat_ecdf(geom = "point", pad = FALSE) +
  stat_function(fun = plnorm, args = list(CI.lnrm$estimate[1], CI.lnrm$estimate[2]),
                aes(colour = "Class I")) +
  stat_function(fun = plnorm, args = list(CII.lnrm$estimate[1], CII.lnrm$estimate[2]),
                aes(colour = "Class II")) +
  stat_function(fun = plnorm, args = list(CIII.lnrm$estimate[1], CIII.lnrm$estimate[2]),
                aes(colour = "Class III")) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 6, base = 10),
    labels = scales::trans_format("log10", scales::math_format(.x))) +
  annotation_logticks(sides = "b") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = percent_format()) + # This puts into % scale to 0 decimal places
  ggtitle("CDFs of Cramer Classes \nfrom ToxVal data") +
  xlab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  ylab("Cumulative Frequency") +
  scale_colour_manual(values = class.colours, name = "Cramer Class", breaks = c("Class I", "Class II", "Class III"), 
                      labels = c("Class I", "Class II", "Class III")) +
  theme_fivethirtyeight() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_rect(fill = "transparent"),
        #plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(1,1), 
        legend.position = c(1,0.2), 
        legend.title = element_blank(), 
        legend.direction = "vertical")
ggsave(here(Output, "TTC_v7.png"), device = "png", width = 6.5, height = 7, units = "in", dpi = 600, bg = "transparent")


######################                  ######################                  ######################                  ######################
###   Use Cramer Class designations to calculate 5th %ile NO(A)EL values    ###
###   These can then be used to calculate TTC thresholds    ###


## Class I
quantile(fitdist(Class.I.min$toxval_numeric, "lnorm"), probs = 0.05)
## Results @ 5%ile NOAEL min (565 obs) = 3.73

## Class II
quantile(fitdist(Class.II.min$toxval_numeric, "lnorm"), probs = 0.05)
## Results @ 5%ile NOAEL min (39 obs) = 3.46

## Class III
quantile(fitdist(Class.III.min$toxval_numeric, "lnorm"), probs = 0.05)
## Results @ 5%ile NOAEL min (700 obs) = 0.39


######################                  ######################                  ######################                  ######################
###   Munro dataset   ###
###   Calculate 5%ile NO(A)EL values for each Cramer class    ###
###   Use Kolmorogov-Smirnov test to compare distributions between Cramer classes between ToxVal and Munro datasets   ###


###   Comparing ToxVal and Munro datasets for each Cramer class   ###
## List of chemical names that also require having their NO(A)EL divided by 3
exceptions <- c("Bis(2-ethylhexyl)phthalate", "Di(2-ethylhexyl)adipate", "Dimethyl terephthalate", "Dimethyldicarbonate", "Isopropyl alcohol",
                "Tertiary butyl hydroquinone", "Xylitol", "Acrilic acid", "Butylated hydroxytoluene", "Ethyl maltol", "Aldicarb", 
                "Chlorofructose, 6-", "Ipronidazole", "Methyl-4-chlorophenoxyacetic acid, 2-", "Metolachlor", "Pravadoline")

## Read in Munro data
Munro <- read_csv(here(Input, "Munro_dataset.csv"), col_types = cols()) %>%
  clean_names() %>% 
  rename(., name = name_munro_1996,
         casrn = cas_original) %>% 
  mutate(study_type = trimws(study_type, which = "right")) %>% # NEED TO REMOVE TRAILING WHITESPACE - helps with ifelse() below
  mutate(structural_class = case_when(structural_class == "1" ~ "Class I",
                                      structural_class == "2" ~ "Class II",
                                      structural_class == "3" ~ "Class III")) %>% 
  # If study_type is subchronic or part of named exceptions above divide munro by 3
  mutate(calculated_noel = case_when(study_type == "sub" ~ (noel_calculated_munro_mg_kg_day / 3),
                                     name %in% exceptions ~ (noel_calculated_munro_mg_kg_day / 3),
                                     TRUE ~ noel_calculated_munro_mg_kg_day),
         log_noel = log10(calculated_noel)) %>% 
  # Read in csv with InChI keys - need to use this to compare overlap with ToxVal
  inner_join(read_csv(here(Input, "Munro_InChi.csv"), col_types = cols()) %>% 
               select(NAME_Munro_1996, InChI_Key), 
             by = c("name" = "NAME_Munro_1996"))

## Separate data into respective Cramer classes
Munro.Class.I <- Munro %>%
  filter(structural_class == "Class I")

Munro.Class.II <- Munro %>% 
  filter(structural_class == "Class II")

Munro.Class.III <- Munro %>% 
  filter(structural_class == "Class III")

## Calculation of 5th %ile NO(A)EL
## Class I
quantile(fitdist(Munro.Class.I$calculated_noel, "lnorm"), probs = 0.05)
## Results @ 5%ile NOAEL min (137 obs) = 2.94

## Class II
quantile(fitdist(Munro.Class.II$calculated_noel, "lnorm"), probs = 0.05)
## Results @ 5%ile NOAEL min (28 obs) = 0.95

## Class III
quantile(fitdist(Munro.Class.III$calculated_noel, "lnorm"), probs = 0.05)
## Results @ 5%ile NOAEL min (448 obs) = 0.15


## Compare ToxVal Class I with Munro Class I distribuions
## Class I
ks.test(Class.I.min$toxval_numeric, Munro.Class.I$calculated_noel, alternative = "two.sided")
# D = 0.12121, p-value = 0.07834

## Class II
ks.test(Class.II.min$toxval_numeric, Munro.Class.II$calculated_noel, alternative = "two.sided")
# D = 0.27289, p-value = 0.1764

## Class III
ks.test(Class.III.min$toxval_numeric, Munro.Class.III$calculated_noel, alternative = "two.sided")
# D = 0.20455, p-value = 2.361e-10


######################                  ######################                  ######################                  ######################
### Identify overlap of chemicals in ToxVal and Munro ###


## Join ToxVal with InChI Keys to compare overlap with Munro
Class.com <- Class.com %>%  
  # Read in csv with InChI keys - only select DTXSID and InChI_Key cols
  inner_join(read_csv(here(Input, "ToxVal_InChI.csv"), col_types = cols()) %>% 
               select(DTXSID, InChI_Key), by = "DTXSID")

## Use InChI Keys to ID overlap between ToxVal and Munro
overlap <- Class.com %>% 
  select(-casrn) %>% 
  inner_join(Munro %>% 
               # Only keep these cols from Munro
               select(structural_class, calculated_noel, log_noel, InChI_Key), 
             by = "InChI_Key") %>% 
  # Rename cols ending w/ ".x" to "_toxval"
  rename_at(vars(ends_with(".x")), list( ~ sub(".x", "_toxval", .))) %>% 
  # Rename cols ending w/ ".y" to "_munro"
  rename_at(vars(ends_with(".y")), list( ~ sub(".y", "_munro", .))) %>% 
  mutate(noel_diff = log_noel_munro - log_noel_toxval)
#write_csv(overlap, here(Output, "Chemical_Overlap_ToxVal_Munro.csv"))

## Scatter plot of ToxVal vs Munro NOAEL values
ggplot(overlap, aes(x = log10(toxval_numeric), y = log10(calculated_noel))) +
  geom_jitter() +
  # Add zero variance line
  geom_abline(aes(slope = 1, intercept = 0, linetype = "solid"), alpha = 0.3, show.legend = TRUE) +
  # Add +0.5 log unit variance line
  geom_abline(aes(slope = 1, intercept = 0.5, linetype = "dashed"), alpha = 0.3, show.legend = TRUE) +
  # Add -0.5 log unit variance line
  geom_abline(aes(slope = 1, intercept = -0.5, linetype = "dashed"), alpha = 0.3, show.legend = TRUE) +
  geom_vline(xintercept = 0, alpha = 0.3, show.legend = FALSE) +
  geom_hline(yintercept = 0, alpha = 0.3, show.legend = FALSE) +
  xlab(expression(" ToxValDB NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  ylab(expression(" Munro NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  scale_linetype_manual(labels = c("±0.5 log \nvariance", "Zero \nvariance"), values = c("dashed", "solid")) +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(1,1), 
        legend.position = c(0.96, 0.3), 
        legend.title = element_blank(), 
        legend.direction = "vertical",
        legend.key.size = unit(0.85, "cm"))
ggsave(here(Output, "Overlap_Toxval_Munro_scatter.png"), device = "png", width = 7, height = 7, units = "in", dpi = 600)

## Use K-S test to compare ToxVal and Munro NOAELs for overlapping chemicals
ks.test(overlap$toxval_numeric, overlap$calculated_noel, alternative = "two.sided")
# D = 0.10502, p-value = 0.1785

## Calculate lognormal fitted distribution for Munro and ToxVal
M.lnrm <- fitdistr(overlap$calculated_noel, "lognormal")
T.lnrm <- fitdistr(overlap$toxval_numeric, "lognormal")

## Plot ECDF for overlapping chemicals
overlap %>% 
  select(DTXSID, toxval_numeric, calculated_noel) %>% 
  gather("source", "value", -DTXSID) %>% 
  mutate(source = case_when(source == "toxval_numeric" ~ "ToxVal",
                            source == "calculated_noel" ~ "Munro")) %>% 
  ggplot(., aes(x = value, shape = source)) +
  geom_point(stat = "ecdf") +
  stat_function(fun = plnorm, linetype = "dashed", args = list(M.lnrm$estimate[1], M.lnrm$estimate[2])) +
  stat_function(fun = plnorm, args = list(T.lnrm$estimate[1], T.lnrm$estimate[2])) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 6, base = 10),
    labels = scales::trans_format("log10", scales::math_format(.x))) +
  annotation_logticks(sides = "b") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = percent_format()) + # This puts into % scale to 0 decimal places
  ggtitle("CDFs of overlapping chemicals \nfrom ToxValDB and Munro data") +
  xlab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  ylab("Cumulative Frequency") +
  scale_shape_manual(values = c("ToxVal" = 16, "Munro" = 4), 
                     name = "Source", 
                     labels = c("ToxVal", "Munro"), 
                     breaks = c("ToxVal", "Munro")) +
  theme_fivethirtyeight() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_rect(fill = "transparent"),
        #plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(1,1), 
        legend.position = c(1,0.2), 
        legend.title = element_blank(), 
        legend.direction = "vertical")
ggsave(here(Output, "Overlap_Toxval_Munro_CDF.png"), device = "png", width = 6.5, height = 7, units = "in", dpi = 600)


######################                  ######################                  ######################                  ######################
###   Generate Cumulative Distribution Frequency graph comparing   ###
###   each Cramer Classification for ToxVal and Munro datasets   ###


## Rename toxval_numeric col so easier to combine data sets
Class.com <- Class.com %>% 
  rename(calculated_noel = toxval_numeric)

## Combine ToxVal and Munro data
ToxVal.Munro <- Munro %>%
  select(structural_class, calculated_noel, casrn, name) %>% 
  # Bind dfs together & add source col - this will add names rather than number
  bind_rows("Munro" = ., "ToxVal" = Class.com, .id = "source") %>% 
  mutate(log_noel = log10(calculated_noel)) 
  

munro1.lnrm <- fitdistr(Munro.Class.I$calculated_noel, "lognormal")
munro2.lnrm <- fitdistr(Munro.Class.II$calculated_noel, "lognormal")
munro3.lnrm <- fitdistr(Munro.Class.III$calculated_noel, "lognormal")

## Plot each cumulative distribution function on one plot
## By putting "shape" in aes() generates ecdf for ToxVal and Munro separately
ggplot(ToxVal.Munro, aes(x = calculated_noel, colour = structural_class, shape = source)) +
  geom_point(stat = "ecdf") +
  stat_function(fun = plnorm, args = list(CI.lnrm$estimate[1], CI.lnrm$estimate[2]),
                aes(colour = "Class I")) +
  stat_function(fun = plnorm, args = list(CII.lnrm$estimate[1], CII.lnrm$estimate[2]),
                aes(colour = "Class II")) +
  stat_function(fun = plnorm, args = list(CIII.lnrm$estimate[1], CIII.lnrm$estimate[2]),
                aes(colour = "Class III")) +
  stat_function(fun = plnorm, linetype = "dashed", args = list(munro1.lnrm$estimate[1], munro1.lnrm$estimate[2]),
                aes(colour = "Class I")) +
  stat_function(fun = plnorm, linetype = "dashed", args = list(munro2.lnrm$estimate[1], munro2.lnrm$estimate[2]),
                aes(colour = "Class II")) +
  stat_function(fun = plnorm, linetype = "dashed", args = list(munro3.lnrm$estimate[1], munro3.lnrm$estimate[2]),
                aes(colour = "Class III")) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 6, base = 10),
    labels = scales::trans_format("log10", scales::math_format(.x))) +
  annotation_logticks(sides = "b") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = percent_format()) + # This puts into % scale to 0 decimal places
  xlab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  ylab("Cumulative Frequency") +
  scale_shape_manual(values = c("ToxVal" = 16, "Munro" = 4), 
                     name = "Source", 
                     labels = c("ToxVal", "Munro"), 
                     breaks = c("ToxVal", "Munro")) +
  scale_colour_manual(values = class.colours, name = "Cramer Class", breaks = c("Class I", "Class II", "Class III"), 
                      labels = c("Class I", "Class II", "Class III")) +
  theme_fivethirtyeight() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_rect(fill = "transparent"),
        #plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(1,1), 
        legend.position = c(1,0.5), 
        legend.title = element_text(), 
        legend.direction = "vertical")
ggsave(here(Output, "Toxval_Munro_CDF_v7.png"), device = "png", width = 6.5, height = 7, units = "in", dpi = 600)


## Generate CDF for each Cramer Class & use grid.arrange to plot them together
## Class I
c1 <- ToxVal.Munro %>% 
  filter(structural_class == "Class I") %>% # Filter for only Class I chemicals
  ggplot(., aes(x = calculated_noel, colour = structural_class, shape = source)) +
  geom_point(stat = "ecdf") +
  stat_function(fun = plnorm, args = list(CI.lnrm$estimate[1], CI.lnrm$estimate[2]),
                aes(colour = "Class I")) +
  stat_function(fun = plnorm, linetype = "dashed", args = list(munro1.lnrm$estimate[1], munro1.lnrm$estimate[2]),
                aes(colour = "Class I")) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 6, base = 10),
    labels = scales::trans_format("log10", scales::math_format(.x))) +
  annotation_logticks(sides = "b") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = percent_format()) + # This puts into % scale to 0 decimal places
  ggtitle("Cramer Class I") +
  xlab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  ylab("Cumulative Frequency") +
  scale_shape_manual(values = c("ToxVal" = 16, "Munro" = 4), 
                     name = "Source", 
                     labels = c("ToxVal", "Munro"), 
                     breaks = c("ToxVal", "Munro")) +
  scale_colour_manual(values = class.colours[1], name = "Cramer Class") +
  theme_fivethirtyeight() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_rect(fill = "transparent"),
        #plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(1,1), 
        legend.position = c(1,0.25), 
        legend.title = element_text(), 
        legend.direction = "vertical") +
  guides(colour = FALSE) # Removes colour from legend
#ggsave(here(Output, "Toxval_Munro_Class_I_CDF.png"), device = "png", width = 9, height = 6, units = "in", dpi = 600)

## Class II
c2 <- ToxVal.Munro %>% 
  filter(structural_class == "Class II") %>% # Filter for only Class II chemicals
  ggplot(., aes(x = calculated_noel, colour = structural_class, shape = source)) +
  geom_point(stat = "ecdf") +
  stat_function(fun = plnorm, args = list(CII.lnrm$estimate[1], CII.lnrm$estimate[2]),
                aes(colour = "Class II")) +
  stat_function(fun = plnorm, linetype = "dashed", args = list(munro2.lnrm$estimate[1], munro2.lnrm$estimate[2]),
                aes(colour = "Class II")) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 6, base = 10),
    labels = scales::trans_format("log10", scales::math_format(.x)),
    limits = c(1e-2, 1e4)) +
  annotation_logticks(sides = "b") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = percent_format()) + # This puts into % scale to 0 decimal places
  ggtitle("Cramer Class II") +
  xlab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  ylab("Cumulative Frequency") +
  scale_shape_manual(values = c("ToxVal" = 16, "Munro" = 4), 
                     name = "Source", 
                     labels = c("ToxVal", "Munro"), 
                     breaks = c("ToxVal", "Munro")) +
  scale_colour_manual(values = class.colours[2], name = "Cramer Class") +
  theme_fivethirtyeight() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_rect(fill = "transparent"),
        #plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(1,1), 
        legend.position = c(1,0.25), 
        legend.title = element_text(), 
        legend.direction = "vertical") +
  guides(colour = FALSE)
#ggsave(here(Output, "Toxval_Munro_Class_II_CDF.png"), device = "png", width = 9, height = 6, units = "in", dpi = 600)

## Class III
c3 <- ToxVal.Munro %>% 
  filter(structural_class == "Class III") %>% # Filter for only Class III chemicals
  ggplot(., aes(x = calculated_noel, colour = structural_class, shape = source)) +
  geom_point(stat = "ecdf") +
  stat_function(fun = plnorm, args = list(CIII.lnrm$estimate[1], CIII.lnrm$estimate[2]),
                aes(colour = "Class III")) +
  stat_function(fun = plnorm, linetype = "dashed", args = list(munro3.lnrm$estimate[1], munro3.lnrm$estimate[2]),
                aes(colour = "Class III")) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 6, base = 10),
    labels = scales::trans_format("log10", scales::math_format(.x))) +
  annotation_logticks(sides = "b") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = percent_format()) + # This puts into % scale to 0 decimal places
  ggtitle("Cramer Class III") +
  xlab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  ylab("Cumulative Frequency") +
  scale_shape_manual(values = c("ToxVal" = 16, "Munro" = 4), 
                     name = "Source", 
                     labels = c("ToxVal", "Munro"), 
                     breaks = c("ToxVal", "Munro")) +
  scale_colour_manual(values = class.colours[3], name = "Cramer Class") +
  theme_fivethirtyeight() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_rect(fill = "transparent"),
        #plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(1,1), 
        legend.position = c(1,0.25), 
        legend.title = element_text(), 
        legend.direction = "vertical") +
  guides(colour = FALSE)
#ggsave(here(Output, "Toxval_Munro_Class_III_CDF.png"), device = "png", width = 9, height = 6, units = "in", dpi = 600)

## Arrange the 3 plots together
all_cdf <- plot_grid(c1, c2, c3, ncol = 3)
save_plot(here(Output, "ToxVal_Munro_All_CDF.png"), all_cdf,  ncol = 3, base_height = 6)


######################                  ######################                  ######################                  ######################
###   Calculate whether the log of the 5th %ile NO(A)EL values are different   ###
###   for each Cramer Classification for ToxVal and Munro datasets   ###


CI_fifth_comp <- ToxVal.Munro %>% 
  filter(structural_class == "Class I") #%>% 
  #mutate(log_noel = log10(calculated_noel))
qcomhd(log_noel ~ source, data = CI_fifth_comp, q = .05, nboot = 5000)
#   q  n1  n2   est1   est2 est1-est.2 ci.low ci.up p.crit p.value
#0.05 137 565 0.3923 0.3926     -3e-04 -0.701 0.361   0.05    0.95


CII_fifth_comp <- ToxVal.Munro %>% 
  filter(structural_class == "Class II") #%>% 
  #mutate(log_noel = log10(calculated_noel))
qcomhd(log_noel ~ source, data = CII_fifth_comp, q = .05, nboot = 5000)
#   q n1 n2    est1   est2 est1-est.2  ci.low  ci.up p.crit p.value
#0.05 28 39 -0.0825 0.3687    -0.4512 -1.1389 0.2185   0.05   0.214

CIII_fifth_comp <- ToxVal.Munro %>% 
  filter(structural_class == "Class III") #%>% 
  #mutate(log_noel = log10(calculated_noel))
qcomhd(log_noel ~ source, data = CIII_fifth_comp, q = .05, nboot = 5000)
#   q  n1  n2    est1    est2 est1-est.2  ci.low  ci.up p.crit p.value
#0.05 448 700 -0.9001 -0.5678    -0.3322 -0.6725 -0.014   0.05  0.0412

## Generate density plots of log NO(A)ELs for ToxVal & Munro
ggplot(CI_fifth_comp, aes(x = log_noel, colour = source)) +
  geom_density() +
  ggtitle("Density function log10 NO(A)EL \nClass I chemicals") +
  theme_fivethirtyeight() +
  scale_color_fivethirtyeight() +
  theme(axis.title = element_text())

ggplot(CII_fifth_comp, aes(x = log_noel, colour = source)) +
  geom_density() +
  ggtitle("Density function log10 NO(A)EL \nClass II chemicals") +
  theme_fivethirtyeight() +
  scale_color_fivethirtyeight() +
  theme(axis.title = element_text())

ggplot(CIII_fifth_comp, aes(x = log_noel, colour = source)) +
  geom_density() +
  ggtitle("Density function log10 NO(A)EL \nClass III chemicals") +
  theme_fivethirtyeight() +
  scale_color_fivethirtyeight() +
  theme(axis.title = element_text())


######################                  ######################                  ######################                  ######################
###   Calculate the confidence intervals for each of the 5th %ile   ###
###   NO(A)EL values in both ToxVal and Munro & plot them as   ###
###   a point graph with associated error bars    ###

fifth_CI_boot <- function(x, prob = 0.05, iter = 5000) {
  bln <- bootdist(fitdist(x, "norm"), bootmethod = "param", niter = iter)
  quantile(bln, probs = prob)
}

## ToxVal Class I chemicals
ToxVal_CI_boot <- fifth_CI_boot(Class.I.min$log_noel)

## ToxVal Class II chemicals
ToxVal_CII_boot <- fifth_CI_boot(Class.II.min$log_noel)

## ToxVal Class III chemicals
ToxVal_CIII_boot <- fifth_CI_boot(Class.III.min$log_noel)

## Munro Class I chemicals
Munro_CI_boot <- fifth_CI_boot(Munro.Class.I$log_noel)

## Munro Class II chemicals
Munro_CII_boot <- fifth_CI_boot(Munro.Class.II$log_noel)

## Munro Class III chemicals
Munro_CIII_boot <- fifth_CI_boot(Munro.Class.III$log_noel)

h <- data.frame("source" = c("ToxVal", "Munro"), 
                "struc_class" = c("Class I", "Class I", "Class II", "Class II", "Class III", "Class III"),
                "fifth_percentile" = c(0.572, 0.468, 0.539, -0.021, -0.405, -0.815), 
                "lower" = c(0.473, 0.233, 0.176, -0.490, -0.529, -0.973), 
                "upper" = c(0.680, 0.725, 0.936, 0.480, -0.276, -0.657))
## This will take the information from the bootstrapping and generate a data.frame that can
## be used below to generate a dot & whisker plot
fifth.CI <- data.frame("source" = c("ToxVal", "ToxVal", "ToxVal", "Munro", "Munro", "Munro"), 
                "struc_class" = c("Class I", "Class II", "Class III", "Class I", "Class II", "Class III"),
                "fifth_percentile" = c(ToxVal_CI_boot$quantmedian, ToxVal_CII_boot$quantmedian, ToxVal_CIII_boot$quantmedian,
                                       Munro_CI_boot$quantmedian, Munro_CII_boot$quantmedian, Munro_CIII_boot$quantmedian), 
                "lower" = c(ToxVal_CI_boot$quantCI[1,], ToxVal_CII_boot$quantCI[1,], ToxVal_CIII_boot$quantCI[1,], 
                            Munro_CI_boot$quantCI[1,], Munro_CII_boot$quantCI[1,], Munro_CIII_boot$quantCI[1,]), 
                "upper" = c(ToxVal_CI_boot$quantCI[2,], ToxVal_CII_boot$quantCI[2,], ToxVal_CIII_boot$quantCI[2,], 
                            Munro_CI_boot$quantCI[2,], Munro_CII_boot$quantCI[2,], Munro_CIII_boot$quantCI[2,]))

ggplot(fifth.CI, aes(x = struc_class, colour = source)) +
  geom_point(aes(y = fifth_percentile), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(x = struc_class, ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
  xlab("Cramer Class") +
  ylab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  theme_fivethirtyeight() +
  theme(axis.title = element_text()) +
  scale_colour_fivethirtyeight()
ggsave(here(Output, "Toxval_Munro_5th_Percentile_CI.png"), device = "png", width = 9, height = 6, units = "in", dpi = 600)


######################                  ######################                  ######################                  ######################
###   Compare ToxPrints between ToxVal &  ###
###   Munro Cramer class III chemicals    ###


## Extract SMILES from ToxVal
ToxVal %>% 
  select(DTXSID, name, qsar_ready_smiles) %>% 
  write_tsv(., here(Inter, "ToxVal_v7_SMILES.tsv"), quote_escape = FALSE)

## Extract SMILES from Munro
Munro %>% 
  select(id_munro, name, smiles) %>% 
  write_tsv(., here(Inter, "Munro_SMILES.tsv"), quote_escape = FALSE)


####                                          ####                                      ####
##    Pass Munro_SMILES tsv through KNIME to dearomatise & desalt SMILES & generate SDF   ##
##    - ToxVal doesn't need this as they are already dearomatised & desalted              ##
##    Run through the ChemoTyper to get toxprints                                         ##
####                                          ####                                      ####


## Read in ToxVal Toxprints & filter for Class III chemicals
ToxVal.tp <- read_tsv(here(Inter, "ToxVal_v7_toxprints.tsv"), col_types = cols()) %>%
  rename(., DTXSID = M_NAME) %>% 
  filter(DTXSID %in% Class.III.min$DTXSID)

## Read in Munro Toxprints & filter for Class III chemicals
Munro.tp <- read_tsv(here(Inter, "Munro_toxprints.tsv"), col_types = cols()) %>%
  rename(., name = M_name) %>% 
  filter(name %in% Munro.Class.III$name)

## Read in ToxPrint info to pull out Level 2 ToxPrint names
tp.info <- read_csv(here(Input, "toxprint_V2.0_r711_5Levels.csv"), col_types = cols()) %>%
  mutate(tmp = as.character(`Level 2 full`)) %>%
  distinct(tmp) %>%
  pull(.)

## Function to calculate number of chemicals in data frame that contain particular set of ToxPrints
## This uses Level 2 ToxPrint names to gather columns
partialStringMax <- function(dataFrame, partialString) {
  tmp <- dataFrame[grepl(partialString, names(dataFrame), fixed = TRUE)]  # Find string in colnames & generate temp df
  if (ncol(tmp) == 0) { return("NA") # If no col is found w/string return NA
    } 
  else {
    return(tmp[,colSums(tmp) == max(colSums(tmp)), drop = FALSE]) # Else return
    } 
}

## Generate data frame comparing chemicals in ToxVal and Munro Cramer class III 
## Set up blank data frame
comparison <- data.frame()
## For each Level 2 chemotype name run partialStringMax function - calculate frequency of that chemotype in dataset
##
for (i in 1:length(tp.info)) {
  m <- colSums(partialStringMax(Munro.tp, tp.info[i]), na.rm = TRUE) / nrow(Munro.tp)
  new.df <- data.frame(chemotype = tp.info[i], freq = m, data = "Munro")
  comparison <- rbind(comparison, new.df)
}

for (i in 1:length(tp.info)) {
  m <- colSums(partialStringMax(ToxVal.tp, tp.info[i]), na.rm = TRUE) / nrow(ToxVal.tp)
  new.df <- data.frame(chemotype = tp.info[i], freq = m, data = "ToxVal")
  comparison <- rbind(comparison, new.df)
}
## Keep unique rows (using all rows)
comparison <- distinct(comparison, .keep_all = TRUE)

## Generate bar plot that compares frequencies of Cramer Class III chemicals
ggplot(comparison) +
  geom_col(aes(x = chemotype, y = freq, fill = data),
           position = "dodge") + 
  coord_flip() +
  xlim(rev(levels(comparison$chemotype))) +
  xlab("ToxPrint") +
  ylab("Frequency") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = percent_format(accuracy = 1)) +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(),
        legend.justification = c(1,1), 
        legend.position = c(1,0.55), 
        legend.title = element_blank(), 
        legend.direction = "vertical") +
  scale_fill_fivethirtyeight()
ggsave(here(Output, "ToxVal_Munro_Class_III_Comparison.png"), device = "png", 
       width = 7.25, height = 9, units = "in", dpi = 600, bg = "transparent")


######################                  ######################                  ######################                  ######################
### Extract ToxVal & Munro Class III chems ###
### Run through Ryan's chemotype enrichment workflow ###
### Remove chemicals with enriched chemotypes & redo analysis ###


## Read in ToxVal Toxprints & keep only Class III chemicals
ToxVal.tp <- read.delim(here(Inter, "ToxVal_v7_toxprints.tsv"), check.names = FALSE) %>%
  rename(., DTXSID = M_NAME) %>% 
  filter(DTXSID %in% Class.III.min$DTXSID) %>% # Keep on Class III chemicals
  inner_join(.,Class.III.min[,c("DTXSID", "name")], by = "DTXSID") %>% # Join DTXSID to name
  mutate(from_munro = 0) %>% # Add column indicating chems not from Munro data set
  select(name, from_munro, everything(), -DTXSID) # Remove DTXSIDs

## Read in Munro Toxprints & keep only Class III chemicals
Munro.tp <- read.delim(here(Inter, "Munro_toxprints.tsv"), check.names = FALSE) %>%
  rename(., name = M_name) %>% 
  filter(name %in% Munro.Class.III$name) %>% 
  mutate(from_munro = 1) %>% 
  select(name, from_munro, everything())

## Join all rows from ToxVal.tp and Munro.tp for chemotype enrichment analysis - CLASS III CHEMICALS ONLY
all.tp <- Munro.tp %>% 
  full_join(., ToxVal.tp) %>% 
  rowid_to_column() %>% 
  select(name, rowid, from_munro, everything()) %>% 
  write_table(all.tp, here(Inter, "ToxVal_Munro_Toxprints.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

####                        ####                        ####
##    Run through Ryan's chemotype enrichment workflow    ##
####                        ####                        ####

## Read in results from chemotype enrichment
enrich <- read.delim(here(Inter, "enrichment_table_result.tsv"), check.names = FALSE)

## ID chemotypes only present in Munro & not in ToxVal
chem_enrich <- enrich %>% 
  filter(`Odds Ratio` == "Inf") %>% 
  select(Fingerprint_ID) %>% 
  deframe()


tmp <- data.frame()
## Using chemotypes from chem_enrich
## ID chems from Munro that contain any of those chemotypes
for (chemotype in chem_enrich) {
  df1 <- Munro.tp %>% 
    filter(get(chemotype) == 1)
  
  tmp <- rbind(tmp, df1)
}
## Chems may have multiple chemotypes from list - keep only unique chems
tmp <- tmp %>% 
  distinct(name, .keep_all = TRUE)

## Filter Munro Class III list for chems NOT ID'd above
## i.e. chemicals that DON'T contain ToxPrints not in ToxVal
Munro.no.Toxprint <- Munro.Class.III %>% 
  filter(!name %in% tmp$name)

## Calculate 5th %ile
quantile(fitdist(Munro.no.Toxprint$calculated_noel, "lnorm"), probs = 0.05)
# Results @ 5%ile NOAEL (306 obs) =  0.22 mg/kg/day

## Are the distributions statistically different?
ks.test(Class.III.min$toxval_numeric, Munro.no.Toxprint$calculated_noel, alternative = "two.sided")
# D = 0.22388, p-value = 1.074e-09

ToxVal.Munro.C3 <- Munro.no.Toxprint %>%
  select(structural_class, calculated_noel, casrn, name) %>% 
  # Bind dfs together & add source col - this will add names rather than number
  bind_rows("Munro" = ., "ToxVal" = Class.com, .id = "source") %>% 
  mutate(log_noel = log10(calculated_noel)) %>% 
  filter(structural_class == "Class III") 

qcomhd(log_noel ~ source, data = ToxVal.Munro.C3, q = .05, nboot = 5000)
#   q  n1  n2    est1    est2 est1-est.2  ci.low  ci.up p.crit p.value
#0.05 306 700 -0.6777 -0.5678    -0.1098 -0.4402 0.2004   0.05  0.4988

fifth_CI_boot(ToxVal.Munro.C3$log_noel)


######################                  ######################                  ######################                  ######################
## Compare CDFs of Chemicals with ToxPrints only present in Munro and those
## chemicals that do NOT contain ToxPrints only present in Munro
## How do they compare, is there any overlap?


## ID chemicals that contain ToxPrints only present in Munro
Munro.ToxPrint <- Munro.Class.III %>% 
  filter(name %in% tmp$name)

## Use K-S test to compare distributions between chems w/ and w/out ToxPrints ID above
ks.test(Munro.ToxPrint$calculated_noel, Munro.no.Toxprint$calculated_noel, alternative = "two.sided")
# D = 0.10292, p-value = 0.2557

## Calculate lognormal distributions
No.toxprint.lnrm <- fitdistr(Munro.no.Toxprint$calculated_noel, "lognormal")
Toxprint.lnrm <- fitdistr(Munro.ToxPrint$calculated_noel, "lognormal")

## Plot ECDF of chemicals that contain ToxPrints only found in Munro
## and those chemicals that don't contain those ToxPrints
Munro.Class.III %>% 
  mutate(source = if_else(name %in% tmp$name, "Contains ToxPrint", "No ToxPrint")) %>% 
  select(source, everything()) %>% 
  ggplot(aes(x = calculated_noel, shape = source)) +
  geom_point(stat = "ecdf") +
  stat_function(fun = plnorm, linetype = "dashed", args = list(No.toxprint.lnrm$estimate[1], No.toxprint.lnrm$estimate[2])) +
  stat_function(fun = plnorm, args = list(Toxprint.lnrm$estimate[1], Toxprint.lnrm$estimate[2])) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 6, base = 10),
    labels = scales::trans_format("log10", scales::math_format(.x))) +
  annotation_logticks(sides = "b") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = percent_format()) + # This puts into % scale to 0 decimal places
  ggtitle("CDF of Munro Cramer class III split \nby specific ToxPrints") +
  xlab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  ylab("Cumulative Frequency") +
  scale_shape_manual(values = c("Contains ToxPrint" = 16, "No ToxPrint" = 4), 
                     name = "Source", 
                     labels = c("Contains ToxPrint", "No ToxPrint"), 
                     breaks = c("Contains ToxPrint", "No ToxPrint")) +
  theme_fivethirtyeight() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_rect(fill = "transparent"),
        #plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(1,1), 
        legend.position = c(1,0.3), 
        legend.title = element_text(), 
        legend.direction = "vertical")


