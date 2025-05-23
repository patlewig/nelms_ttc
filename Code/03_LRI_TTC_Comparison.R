########################################################################################
#
# Author: Mark D. Nelms, Ph.D., nelms.mark@epa.gov
#
# Version: 1.0 28th June 2019
#
# Description: Profiling the data from the LRI/Scitovation app through Toxtree v3.1.0 
#           to compare the assigned TTC category when using the old OP and carbamate SMARTS 
#           patterns to the results for the same set of chemicals using the our updated 
#           OP/carbamate SMARTS patterns
#           
#
# Notes:
#
#
# Potential Issues: None known
#########################################################################################
library(janitor)
library(here)
library(ggthemes)
library(tidyverse)

Input <- "InputData"
Inter <- "IntermediateData"
Output <- "OutputFiles"


## Read in LRI csv file that Grace downloaded from the web app and save as tsv to run through KNIME to generate SDF
## Can't select "name" for use in tsv file - name for DTXSID90849581 gets split over multiple lines when put into Excel 
##          this also breaks txt2sdf KNIME workflow
read_csv(here(Input, "TTC_data_LRI.csv"), col_types = cols()) %>% 
  clean_names() %>% 
  rename(DTXSID = dss_tox_substance_id) %>% 
  select(smiles, DTXSID, casrn) %>%
  write_tsv(here(Inter, "LRI_SMILES.tsv"), quote_escape = FALSE)


####                                                 ####                                                     ####
##    Pass tsv through KNIME to generate SDF - Think ChemmineR generates V3000 SDF - not compatible w/Toxtree   ##
##    Pass SDF through Toxtree to run through Cramer, Kroes, Carbamate, OP, and Steroid decision trees          ##
####                                                 ####                                                     ####


######################                  ######################                  ######################                  ######################
###   Assign chemicals to TTC categories based upon new OP/Carbamate Toxtree modules   ###


LRI <- read_csv(here(Input, "TTC_data_LRI.csv"), col_types = cols()) %>% 
  clean_names() %>% 
  rename(DTXSID = dss_tox_substance_id)

## Read in results from Cramer workflow
Cramer <- read_csv(here(Input, "LRI_Cramer_results.csv"), col_types = cols()) %>%
  clean_names() %>% 
  rename(DTXSID = dtxsid) %>%
  filter(!grepl("DTXSID", cramer_rules)) %>% ## Remove chemicals that didn't go through workflow properly
  inner_join(select(LRI, DTXSID, name), by = "DTXSID") %>% 
  select(DTXSID, smiles, name, cramer_rules) ## Keep only these 3 cols in this order

## Read in results from Kroes workflow
Kroes <- read_csv(here(Input, "LRI_Kroes_results.csv"), col_types = cols()) %>%
  clean_names() %>% 
  rename(DTXSID = dtxsid,
         kroes_rules = kroes_ttc_decision_tree, 
         kroes_explanation = contains("explanation")) %>%
  select(DTXSID, kroes_rules, kroes_explanation) 

## Read in results from Carbamate workflow
Carbamate <- read_csv(here(Input, "LRI_Carbamate_results.csv"), col_types = cols()) %>%
  clean_names() %>% 
  rename(DTXSID = dtxsid,
         carbamate_rules = mark_carbamates) %>%
  select(DTXSID, carbamate_rules) 

## Read in results from OP workflow
OP <- read_csv(here(Input, "LRI_OP_results.csv"), col_types = cols()) %>%
  rename(OP_rules = Mark_OPs, 
         OP_explanation = contains("explanation")) %>%
  select(DTXSID, OP_rules, OP_explanation) 

## Read in results from Steroid workflow
Steroid <- read_csv(here(Input, "LRI_Steroid_results.csv"), col_types = cols()) %>%
  rename(steroid_rules = Steroids) %>%
  select(DTXSID, steroid_rules) 

## Read in results from old OP workflow
Old_OP <- read_csv(here(Input, "LRI_Old_OP_results.csv"), col_types = cols()) %>%
  rename(Old_OP_rules = OPs, 
         OP_explanation = contains("explanation")) %>%
  select(DTXSID, Old_OP_rules, OP_explanation) 

## Read in results from old Carbamate workflow
Old_Carbamate <- read_csv(here(Input, "LRI_Old_Carbamate_results.csv"), col_types = cols()) %>%
  rename(Old_carbamate_rules = `Alkyl carbamate and thiocarbamate`) %>%
  select(DTXSID, Old_carbamate_rules) 

## Merge Cramer w/other data by DTXSID
TTC <- Cramer %>%
  inner_join(., Kroes, by = "DTXSID") %>%
  inner_join(., Carbamate, by = "DTXSID") %>%
  inner_join(., OP, by = "DTXSID") %>%
  inner_join(., Steroid, by = "DTXSID") %>% 
  inner_join(., Old_Carbamate, by = "DTXSID") %>%
  inner_join(., Old_OP, by = "DTXSID")

## Filtering those excluded from TTC
Excluded <- TTC %>%
  filter(., kroes_rules == "Risk assessment requires compound-specific toxicity data") %>%
  select(., DTXSID, name, smiles, kroes_rules, kroes_explanation)
#write_csv(Excluded, here(Inter, "LRI_Excluded_from_TTC.csv"))

## Filtering for High potency carcinogens
HighPotCarc <- TTC %>%
  filter(., kroes_rules == "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that contain "Risk ..."
  filter(., str_detect(kroes_explanation, pattern = "Q2Y") & str_detect(kroes_explanation, pattern = "Q3Y")) %>% ## Filter for rows answering Yes to Q2 & Q3
  select(., DTXSID, name, smiles, kroes_rules, kroes_explanation)
#write_csv(HighPotCarc, here(Inter, "LRI_HighPotCarc.csv"))

## Filtering for Steroids
Ster <- TTC %>% 
  filter(., str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that contain "Negligible"
  filter(., str_detect(steroid_rules, "Default Class 1")) %>% ## Filter for rows designated as Carbamates
  select(., DTXSID, name, smiles, kroes_rules, steroid_rules)
#write_csv(Ster, here(Inter, "LRI_TTC_Steroids.csv"))

## Filtering for Organophosphates
Organo <- TTC %>%
  filter(., kroes_rules == "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that contain "Risk ..."
  filter(., str_detect(kroes_explanation, pattern = "Q1Y")) %>% ## Filter for rows answering Yes to Q1
  filter(., OP_rules == "Default Class 1") %>% ## Filter for rows designated as OPs
  select(., DTXSID, name, smiles, kroes_rules, kroes_explanation, OP_rules)
#write_csv(Organo, here(Inter, "LRI_TTC_OPs.csv"))

## Filtering for Carbamates
Carb <- TTC %>%
  #filter(., str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that contain "Negligible"
  filter(., carbamate_rules == "Default Class 1") %>% ## Filter for rows designated as Carbamates
  select(., DTXSID, name, smiles, kroes_rules, carbamate_rules)
#write_csv(Carb, here(Inter, "LRI_TTC_Carbamates.csv"))

## Combining OPs and Carbamates - i.e. anticholinesterase inhibitors
Anti.CHE <- full_join(Carb, Organo)
#write_csv(Anti.CHE, here(Inter, "LRI_Anti_CHE.csv"))

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
#write_csv(NA.for.TTC, here(Inter, "LRI_NA_for_TTC.csv"))

## Filtering for Genetox that requires a TTC of 0.15ug/day
GenTox.TTC <- TTC %>%
  filter(., str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that contain "Negligible"
  filter(., str_detect(steroid_rules, "Default Class 2")) %>% ## Filter for rows that are NOT steroids
  filter(., str_detect(carbamate_rules, "Default Class 2")) ## Filter for rows that are NOT carbamates
#write_csv(GenTox.TTC, here(Inter, "LRI_GenTox_TTC.csv"))

## Filtering for Cramer Class I chemicals
Class.I <- TTC %>%
  filter(., kroes_rules != "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that do NOT contain "Risk ..."
  filter(., !str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that do NOT contain "Negligible"
  filter(str_detect(carbamate_rules, "Default Class 2")) %>% ## Filter for rows that are NOT carbamates
  filter(., str_detect(cramer_rules, "Low")) ## Filter for rows that are in Cramer Class I
#write_csv(Class.I, here(Inter, "LRI_Class_I.csv"))

## Filtering for Cramer Class II chemicals
Class.II <- TTC %>%
  filter(., kroes_rules != "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that do NOT contain "Risk ..."
  filter(., !str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that do NOT contain "Negligible"
  filter(str_detect(carbamate_rules, "Default Class 2")) %>% ## Filter for rows that are NOT carbamates
  filter(., str_detect(cramer_rules, "Intermediate")) ## Filter for rows that are in Cramer Class II
#write_csv(Class.II, here(Inter, "LRI_Class_II.csv"))

## Filtering for Cramer Class III chemicals
Class.III <- TTC %>%
  filter(., kroes_rules != "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that do NOT contain "Risk ..."
  filter(., !str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that do NOT contain "Negligible"
  filter(str_detect(carbamate_rules, "Default Class 2")) %>% ## Filter for rows that are NOT carbamates
  filter(., str_detect(cramer_rules, "High")) ## Filter for rows that are in Cramer Class III
#write_csv(Class.III, here(Inter, "LRI_Class_III.csv"))


######################                  ######################                  ######################                  ######################
###   Assign chemicals to TTC categories based upon OLD OP/Carbamate Toxtree modules   ###


## Add in the results from the old OPs and Carbamates module from Toxtree
## Filtering for Organophosphates
Old_Organo <- TTC %>%
  filter(., kroes_rules == "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that contain "Risk ..."
  filter(., str_detect(kroes_explanation, pattern = "Q1Y")) %>% ## Filter for rows answering Yes to Q1
  filter(., Old_OP_rules == "Default Class 1") %>% ## Filter for rows designated as OPs
  select(., DTXSID, name, smiles, kroes_rules, kroes_explanation, Old_OP_rules)
#write_csv(Organo, here(Inter, "LRI_TTC_OPs.csv"))

## Filtering for Carbamates
Old_Carb <- TTC %>%
  filter(., str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that contain "Negligible"
  filter(., Old_carbamate_rules == "Default Class 1") %>% ## Filter for rows designated as Carbamates
  select(., DTXSID, name, smiles, kroes_rules, Old_carbamate_rules)
#write_csv(Carb, here(Inter, "LRI_TTC_Carbamates.csv"))

## Combining OPs and Carbamates - i.e. anticholinesterase inhibitors
Old_Anti.CHE <- full_join(Old_Carb, Old_Organo)
#write_csv(Anti.CHE, here(Inter, "LRI_Anti_CHE.csv"))

## Filtering for chems that:
## 1) Require compound-specific tox data & answered Yes to Q1,
## 2) are NOT OPs,
## 3) ARE steroids
Old_NA.for.TTC <- TTC %>%
  filter(., kroes_rules == "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that contain "Risk ..."
  filter(., str_detect(kroes_explanation, pattern = "Q1Y")) %>% ## Filter for rows answering Yes to Q1
  filter(., str_detect(Old_OP_rules, "Default Class 2")) %>% ## Filter for rows that are NOT designated as OPs
  full_join(., Ster, by = "DTXSID") %>% ## Join with chemicals ID'd as Steroids
  full_join(., HighPotCarc, by = "DTXSID") ## Join with chemicals pulled out as High Potency Carcinogens
#write_csv(NA.for.TTC, here(Inter, "LRI_NA_for_TTC.csv"))

## Filtering for Genetox that requires a TTC of 0.15ug/day
Old_GenTox.TTC <- TTC %>%
  filter(., str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that contain "Negligible"
  filter(., str_detect(steroid_rules, "Default Class 2")) %>% ## Filter for rows that are NOT steroids
  filter(., str_detect(Old_carbamate_rules, "Default Class 2")) ## Filter for rows that are NOT carbamates
#write_csv(GenTox.TTC, here(Inter, "LRI_GenTox_TTC.csv"))

## Filtering for Cramer Class I chemicals
Old_Class.I <- TTC %>%
  filter(., kroes_rules != "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that do NOT contain "Risk ..."
  filter(., !str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that do NOT contain "Negligible"
  filter(str_detect(Old_carbamate_rules, "Default Class 2")) %>% ## Filter for rows that are NOT carbamates
  filter(., str_detect(cramer_rules, "Low")) ## Filter for rows that are in Cramer Class I
#write_csv(Class.I, here(Inter, "LRI_Class_I.csv"))

## Filtering for Cramer Class II chemicals
Old_Class.II <- TTC %>%
  filter(., kroes_rules != "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that do NOT contain "Risk ..."
  filter(., !str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that do NOT contain "Negligible"
  filter(str_detect(Old_carbamate_rules, "Default Class 2")) %>% ## Filter for rows that are NOT carbamates
  filter(., str_detect(cramer_rules, "Intermediate")) ## Filter for rows that are in Cramer Class II
#write_csv(Class.II, here(Inter, "LRI_Class_II.csv"))

## Filtering for Cramer Class III chemicals
Old_Class.III <- TTC %>%
  filter(., kroes_rules != "Risk assessment requires compound-specific toxicity data") %>% ## Filter for rows that do NOT contain "Risk ..."
  filter(., !str_detect(kroes_rules, "Negligible")) %>% ## Filter for rows that do NOT contain "Negligible"
  filter(str_detect(Old_carbamate_rules, "Default Class 2")) %>% ## Filter for rows that are NOT carbamates
  filter(., str_detect(cramer_rules, "High")) ## Filter for rows that are in Cramer Class III
#write_csv(Class.III, here(Inter, "LRI_Class_III.csv"))


## Add TTC assignments after using old and new OP/Carbamates module
LRI_up <- LRI %>% 
  mutate(mark_ttc_class = case_when(DTXSID %in% NA.for.TTC$DTXSID ~ "NA for TTC", 
                                    DTXSID %in% GenTox.TTC$DTXSID ~"Genotoxicity Alert",
                                    DTXSID %in% Anti.CHE$DTXSID ~ "OP/Carbamate",
                                    DTXSID %in% Class.I$DTXSID ~ "Class I",
                                    DTXSID %in% Class.II$DTXSID ~ "Class II",
                                    DTXSID %in% Class.III$DTXSID ~ "Class III",
                                    TRUE ~ NA_character_)) %>%
  ## Add reason why assigned to NA.for.TTC class - Steroid, High Potency Carcinogen, or answering Yes to Q1 in Kroes workflow
  mutate(reason_for_NA = case_when(DTXSID %in% Ster$DTXSID ~ "Steroid", 
                                   DTXSID %in% HighPotCarc$DTXSID ~ "High Potency Carcinogen",
                                   DTXSID %in% NA.for.TTC$DTXSID & !DTXSID %in% Organo$DTXSID ~ "Yes to Q1",
                                   TRUE ~ NA_character_)) %>% 
  mutate(old_ttc_class = case_when(DTXSID %in% Old_NA.for.TTC$DTXSID ~ "NA for TTC", 
                                   DTXSID %in% Old_GenTox.TTC$DTXSID ~"Genotoxicity Alert",
                                   DTXSID %in% Old_Anti.CHE$DTXSID ~ "OP/Carbamate",
                                   DTXSID %in% Old_Class.I$DTXSID ~ "Class I",
                                   DTXSID %in% Old_Class.II$DTXSID ~ "Class II",
                                   DTXSID %in% Old_Class.III$DTXSID ~ "Class III",
                                   TRUE ~ NA_character_))
write_csv(LRI_up, here(Output, "LRI_TTC_updated.csv"))

## Generate tile plot comparing frequency of assignment after using old and new OP/Carbamate modules
LRI_up %>% 
  select(mark_ttc_class, old_ttc_class) %>% 
  mutate(mark_ttc_class = fct_explicit_na(mark_ttc_class, na_level = "No Assignment"),
         old_ttc_class = fct_explicit_na(old_ttc_class, na_level = "No Assignment")) %>% 
  mutate_all(as_factor) %>%
  group_by(mark_ttc_class, old_ttc_class, .drop = FALSE) %>% 
  summarise(Freq = n()) %>% 
  # Where Freq is 0 replace with NA - will remove 0s from tile plot
  mutate(Freq = case_when(Freq == 0 ~ NA_integer_,
                          TRUE ~ Freq)) %>%
  ggplot(., aes(x = mark_ttc_class, y = old_ttc_class, fill = log10(Freq))) +
  geom_tile() +
  geom_text(aes(label = Freq)) + # Add text with frequency to tiles
  ylab('Original OP/Carbamate Module') +
  xlab('\nUpdated OP/Carbamate Module') +
  scale_fill_gradientn(name = expression("Log"[10]*"Frequency"), colours = c("#ffffff", "#93ADAC", "#314654"), na.value = "grey90") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text())
ggsave(here(Output, "Original_Updated_OP_Carb_Comparison.png"), device = "png", width = 9, height = 9, units = "in", dpi = 600)
  
  
  
  
