########################################################################################
#
# Author: Mark D. Nelms, Ph.D., nelms.mark@epa.gov
#
# Version: 1.0 12th September 2018
#
# Description: Continues Mark_TTC_ToxVal_v7 file investigation of discrepancy of 5th %ile
#             NOAEL values between ToxVal and Munro Cramer Class III chemicals
#             1st part - Set up Munro SMILES for running chemicals through Toxtree
#             2nd part - Calculation of 5th %ile value after removing OPs/carbamates
#                       ID'ed using original SMARTS patterns for Munro Class III chemicals
#             3rd part - Calculation of 5th %ile value after removing OP/carbamates ID'ed
#                       using updated version SMARTS patterns based upon manually ID'd 
#                       OPs/carbamates by Munro et al (1999) and EFSA
#
# Notes:
#
#
# Potential Issues: None known
#########################################################################################
library(scales)
library(here)
library(janitor)
library(fitdistrplus)
library(ggthemes)
library(WRS2)
library(tidyverse)

Input <- "InputData"
Inter <- "IntermediateData"
Output <- "OutputFiles"


###   Read in Munro data from csv file   ###
## List of chemical names that also require having their NO(A)EL divided by 3
exceptions <- c("Bis(2-ethylhexyl)phthalate", "Di(2-ethylhexyl)adipate", "Dimethyl terephthalate", "Dimethyldicarbonate", "Isopropyl alcohol",
                "Tertiary butyl hydroquinone", "Xylitol", "Acrilic acid", "Butylated hydroxytoluene", "Ethyl maltol", "Aldicarb", 
                "Chlorofructose, 6-", "Ipronidazole", "Methyl-4-chlorophenoxyacetic acid, 2-", "Metolachlor", "Pravadoline")

## Read in Munro data
Munro <- read_csv(here(Input, "Munro_dataset.csv"), col_types = cols()) %>%
  clean_names() %>% 
  rename(name = name_munro_1996,
         casrn = cas_original) %>% 
  mutate(study_type = trimws(study_type, which = "right")) %>% # NEED TO REMOVE TRAILING WHITESPACE - helps with ifelse() below
  mutate(structural_class = case_when(structural_class == "1" ~ "Class I",
                                      structural_class == "2" ~ "Class II",
                                      structural_class == "3" ~ "Class III")) %>% 
  # If study_type is subchronic or part of named exceptions above divide munro by 3
  mutate(calculated_noel = ifelse(study_type == "sub", noel_calculated_munro_mg_kg_day / 3, 
                                  ifelse(name %in% exceptions, noel_calculated_munro_mg_kg_day / 3, noel_calculated_munro_mg_kg_day))) %>% 
  mutate(log_noel = log10(calculated_noel))

## Separate data into respective Cramer classes
Munro.Class.III <- Munro %>% 
  filter(structural_class == "Class III")

## Save Class III chemical info & run through KNIME workflow to dearomatise & desalt SMILES
## Then run through Toxtree
write.table(Munro.Class.III[,c(smiles, id_munro, name)], here(Inter, "Munro_Class_III_smiles.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

####                                                 ####                                                     ####
##    Pass tsv through KNIME to dearomatise & desalt SMILES & generate SDF                                      ##
##    - Think ChemmineR generates V3000 SDF - not compatible w/Toxtree                                          ##
##    Pass SDF through Toxtree to run through Cramer, Kroes, Carbamate, OP, and Steroid decision trees          ##
####                                                 ####                                                     ####

## Read in results from Cramer workflow
Cramer <- read_csv(here(Input, "Munro_CIII_Cramer_results.csv"), col_types = cols()) %>% 
  clean_names() %>% 
  rename(SMILES = smiles) %>% 
  filter(!grepl("id_munro", cramer_rules)) %>% ## Remove chemicals that didn't go through workflow properly
  select(id_munro, name, SMILES, cramer_rules) 

## Read in results from Carbamate workflow
Carbamate <- read_csv(here(Input, "Munro_CIII_Carbamate_results.csv"), col_types = cols()) %>%
  clean_names() %>% 
  rename(carbamate_rules = alkyl_carbamate_and_thiocarbamate) %>%
  select(id_munro, carbamate_rules) 

## Read in results from OP workflow
OP <- read_csv(here(Input, "Munro_CIII_OP_results.csv"), col_types = cols()) %>%
  clean_names() %>% 
  rename(OP_rules = o_ps, 
         OP_explanation = contains("explanation")) %>%
  select(id_munro, OP_rules, OP_explanation) 

## Merge Cramer w/other data by DTXSID
TTC <- Cramer %>%
  inner_join(., Carbamate, by = "id_munro") %>%
  inner_join(., OP, by = "id_munro")

## Filtering for Organophosphates
Organo <- TTC %>%
  filter(., OP_rules == "Default Class 1") %>% ## Filter for rows designated as OPs
  select(., id_munro, name, SMILES, OP_rules)

## Filtering for Carbamates
Carb <- TTC %>%
  filter(., str_detect(carbamate_rules, "Default Class 1")) %>% ## Filter for rows designated as Carbamates
  select(., id_munro, name, SMILES, carbamate_rules)

## Combining OPs and Carbamates - i.e. anticholinesterase inhibitors
Anti.CHE <- full_join(Carb, Organo)


######################                  ######################                  ######################                  ######################
###   After removal of OPs/Carbamates what happens to:    ###
###   1) 5th %ile NOAEL
###   2) Comparison of distribution with ToxVal Class III chemicals - statistically different?
###   3) Comparison of 5th %ile NOAEL compared to ToxVal Class III chems - statistically different?
###   4) Plot CDF of Munro class III OPs and non-OPs for comparison


## Read in Class.III.min data created by Mark_TTC_ToxVal_v7.R
Class.III.min <- read_csv(here(Output, "Class_III_min_NOAEL.csv"), col_types = cols()) %>% 
  rename(calculated_noel = toxval_numeric) %>% 
  select(structural_class, calculated_noel, casrn, name)

## Remove OPs/Carbamates & retrieve NOAEL data
Munro.III.no.OP <- TTC %>%
  filter(!id_munro %in% Anti.CHE$id_munro) %>% ## Filter for chemicals that are NOT OPs/carbamates
  #filter(., str_detect(cramer_rules, "High")) %>% ##
  inner_join(Munro.Class.III, by = c("id_munro", "name")) %>% 
  select(id_munro, name, SMILES, everything(), 
         -cramer_rules, -carbamate_rules, -OP_rules, -OP_explanation)

## Calculate 5th %ile NOAEL
quantile(fitdist(Munro.III.no.OP$calculated_noel, "lnorm"), probs = 0.05)
# Results @ 5%ile NOAEL (386 obs) = 0.20 mg/kg/day

## Compare distribution to ToxVal Class III chemicals
ks.test(Class.III.min$calculated_noel, Munro.III.no.OP$calculated_noel, alternative = "two.sided")
# D = 0.17815, p-value = 2.769e-07

## Combine ToxVal and Munro Class III data
CIII_fifth_comp <- Munro.III.no.OP %>%
  select(structural_class, calculated_noel, casrn, name) %>% 
  # Bind dfs together & add source col - this will add names rather than number
  bind_rows("Munro" = ., "ToxVal" = Class.III.min, .id = "source") %>% 
  mutate(log_noel = log10(calculated_noel)) 

## Are the 5th percentile NO(A)ELs significantly different?
qcomhd(log_noel ~ source, data = CIII_fifth_comp, q = .05, nboot = 5000)
#   q  n1  n2    est1    est2 est1-est.2  ci.low  ci.up p.crit p.value
#0.05 386 700 -0.7651 -0.5678    -0.1972 -0.5648 0.1328   0.05  0.2444

## Filter for those chemicals that are considered to be OPs from Munro class III
Munro.III.OP <- Munro.Class.III %>% 
  filter(id_munro %in% Anti.CHE$id_munro) %>% 
  mutate(source = "OPs")

## Calculate 5th %ile value of OPs/Carbamates removed from Munro Class III  
quantile(fitdist(Munro.III.OP$calculated_noel, "lnorm"), probs = 0.05)
# Results @ 5%ile NOAEL (62 obs) = 0.056

## Run K-S test to compare Munro class III w/ and w/out OPs and carbamates
ks.test(Munro.III.OP$calculated_noel, Munro.III.no.OP$calculated_noel, alternative = "two.sided")
# D = 0.31681, p-value = 4.404e-05

## Combine OP and non-OP chemicals ready to be
## plotted as a CDF of OP and non-OP chemicals from Munro class III
Munro.OP.and.no.OP <- Munro.III.no.OP %>% 
  mutate(source = "No OPs") %>% 
  full_join(Munro.III.OP)

## Calculate lognormal distributions for OPs and non-OPs
OP.lnrm <- fitdistr(Munro.III.OP$calculated_noel, "lognormal")
no.OP.lnrm <- fitdistr(Munro.III.no.OP$calculated_noel, "lognormal")

## CDF of OP and non-OP chemicals from Munro class III
ggplot(Munro.OP.and.no.OP, aes(x = calculated_noel, shape = source)) +
  geom_point(stat = "ecdf") +
  stat_function(fun = plnorm, linetype = "dashed", args = list(no.OP.lnrm$estimate[1], no.OP.lnrm$estimate[2])) +
  stat_function(fun = plnorm, args = list(OP.lnrm$estimate[1], OP.lnrm$estimate[2])) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 6, base = 10),
    labels = scales::trans_format("log10", scales::math_format(.x))) +
  annotation_logticks(sides = "b") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = percent_format()) + # This puts into % scale to 0 decimal places
  ggtitle("CDF of Munro Cramer class III split \ninto OPs and non-OPs") +
  xlab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  ylab("Cumulative Frequency") +
  scale_shape_manual(values = c("OPs" = 16, "No OPs" = 4), 
                     name = "Source", 
                     labels = c("OPs", "No OPs"), 
                     breaks = c("OPs", "No OPs")) +
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


######################                  ######################                  ######################                  ######################
###   Using my version of SMARTS remove OPs/Carbamates, what happens to:    ###
###   1) 5th %ile NOAEL
###   2) Comparison of distribution with ToxVal Class III chemicals - statistically different?
###   3) Comparison of 5th %ile NOAEL compared to ToxVal Class III chems - statistically different?
###   4) Plot CDF of Munro class III split into OPs and non-OPs 


## Run Class III chemicals through OP/Carbamate decision tree using my SMARTS

## Read in results from Carbamate workflow
Mark.Carbamate <- read_csv(here(Input, "Munro_CIII_Mark_Carbamate_results.csv"), col_types = cols()) %>%
  rename(carbamate_rules = Mark_Carbamates) %>%
  select(id_munro, carbamate_rules) 

## Read in results from OP workflow
Mark.OP <- read_csv(here(Input, "Munro_CIII_Mark_OP_results.csv"), col_types = cols()) %>%
  rename(OP_rules = Mark_OPs, 
         OP_explanation = contains("explanation")) %>%
  select(id_munro, OP_rules, OP_explanation) 

## Merge Cramer w/other data by DTXSID
TTC <- Cramer %>%
  inner_join(., Mark.Carbamate, by = "id_munro") %>%
  inner_join(., Mark.OP, by = "id_munro") 

## Filtering for Organophosphates - Mark's SMARTS
Organo <- TTC %>%
  filter(OP_rules == "Default Class 1") %>% ## Filter for rows designated as OPs
  select(id_munro, name, SMILES, OP_rules)

## Filtering for Carbamates - Mark's SMARTS
Carb <- TTC %>%
  filter(str_detect(carbamate_rules, "Default Class 1")) %>% ## Filter for rows designated as Carbamates
  select(id_munro, name, SMILES, carbamate_rules)

## Combining OPs and Carbamates - i.e. anticholinesterase inhibitors
Anti.CHE.updated <- full_join(Carb, Organo)

## Remove OPs/Carbamates & retrieve NOAEL data
Munro.III.no.OP.updated <- TTC %>%
  filter(!id_munro %in% Anti.CHE$id_munro) %>% ## Filter for chemicals that are NOT OPs/carbamates
  #filter(., str_detect(cramer_rules, "High")) %>% ##
  inner_join(Munro.Class.III, by = c("id_munro", "name")) %>% 
  select(id_munro, name, SMILES, everything(), 
         -carbamate_rules, -OP_rules, -OP_explanation) 

## Calculate 5th %ile NOAEL
quantile(fitdist(Munro.III.no.OP.updated$calculated_noel, "lnorm"), probs = 0.05)
# Results @ 5%ile NOAEL (397 obs) =  0.23 mg/kg/day

## Compare distribution to ToxVal Class III chemicals
ks.test(Class.III.min$toxval_numeric, Munro.III.no.OP.updated$calculated_noel, alternative = "two.sided")
# D = 0.17924, p-value = 2.905e-07

## Combine ToxVal and Munro Class III data
ToxVal.Munro.updated <- Munro.III.no.OP.updated %>%
  select(structural_class, calculated_noel, casrn, name) %>% 
  # Bind dfs together & add source col - this will add names rather than number
  bind_rows("Munro" = ., "ToxVal" = Class.III.min, .id = "source") %>% 
  mutate(log_noel = log10(calculated_noel)) 

CIII_fifth_comp <- ToxVal.Munro.updated %>% 
  filter(structural_class == "Class III") #%>% 
  #mutate(log_noel = log10(calculated_noel))

## Are the 5th percentile NO(A)ELs significantly different?
qcomhd(log_noel ~ source, data = CIII_fifth_comp, q = .05, nboot = 5000)
#   q  n1  n2  est1    est2 est1-est.2 ci.low  ci.up p.crit p.value
#0.05 377 700 -0.6736 -0.5678    -0.1058 -0.4871 0.1977   0.05  0.4864

## Calculate 5th %ile value of OPs/Carbamates removed from Munro Class III
Munro.III.OP.updated <- Munro.Class.III %>% 
  filter(id_munro %in% Anti.CHE.updated$id_munro) %>% 
  mutate(source = "OPs")

quantile(fitdist(Munro.III.OP.updated$calculated_noel, "lnorm"), probs = 0.05)
# Results @ 5%ile NOAEL (51 obs) = 0.037

## Combine OP and non-OP chemicals ready to be
## plotted as a CDF of OP and non-OP chemicals from Munro class III
Munro.OP.and.no.OP.updated <- Munro.III.no.OP.updated %>% 
  mutate(source = "No OPs") %>% 
  full_join(Munro.III.OP.updated, by = c("id_munro", "name"))

## Calculate lognormal distributions for OPs and non-OPs
OP.lnrm <- fitdistr(Munro.III.OP.updated$calculated_noel, "lognormal")
no.OP.lnrm <- fitdistr(Munro.III.no.OP.updated$calculated_noel, "lognormal")

## CDF of OP and non-OP chemicals from Munro class III
ggplot(Munro.OP.and.no.OP.updated, aes(x = calculated_noel, shape = source)) +
  geom_point(stat = "ecdf") +
  stat_function(fun = plnorm, linetype = "dashed", args = list(no.OP.lnrm$estimate[1], no.OP.lnrm$estimate[2])) +
  stat_function(fun = plnorm, args = list(OP.lnrm$estimate[1], OP.lnrm$estimate[2])) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 6, base = 10),
    labels = scales::trans_format("log10", scales::math_format(.x))) +
  annotation_logticks(sides = "b") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = percent_format()) + # This puts into % scale to 0 decimal places
  ggtitle("CDF of Munro Cramer class III split \ninto OPs and non-OPs") +
  xlab(expression("NO(A)EL Log"[10]*"(mg/kg bw/day)")) +
  ylab("Cumulative Frequency") +
  scale_shape_manual(values = c("OPs" = 16, "No OPs" = 4), 
                     name = "Source", 
                     labels = c("OPs", "No OPs"), 
                     breaks = c("OPs", "No OPs")) +
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


