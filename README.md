README
================
This is the code repository associated with the manuscript: Nelms MD, Pradeep P and Patlewicz. **2019**
Evaluating potential refinements to existing Threshold of Toxicological Concern (TTC) values for environmentally-relevant compounds. *Regul Toxicol Pharmacol* 104505  doi: 10.1016/j.yrtph.2019.104505

To re-run the analysis - three data folders need to be downloaded and placed in the same project folder. These are InputData, IntermediateData and OutputFiles. The files are available on Zenodo at 10.5281/zenodo.15498325

The same information is available here: https://gaftp.epa.gov/Comptox/CCTE_Publication_Data/CCED_Publication_Data/PatlewiczGrace/CompTox-TTC_evaluation_refinements-regtoxpharm/ 

The supplementary information provided in the manuscript is also provide as a separate zip file in the same zenodo address.


## Evaluating potential refinements to existing Threshold of Toxicological Concern (TTC) values for environmentally-relevant chemicals

July 31 2019. Contact: <patlewicz.grace@epa.gov>, <nelms.mark@epa.gov>

This repository contains code and input files associated with the
manuscript of the same title that was published by Nelms et al (2019) in
… see … for more details.

In order to generate the files that are required to completely run this
project the files contained in the Code folder should be run in
numerical order (i.e. 00 -\> 03). All intermediate and output files can
be generated using the information present in the code and input data
folders.

Please note that R(v3.5.2), Toxtree (v3.1.0) and KNIME (v3.4.2) were
used throughout this project. Where you
see

``` r
####                                              ####                                            ####
##                                                                                                  ##
##                                                                                                  ##
####                                              ####                                            ####
```

indicates instances where external programs (Toxtree and/or KNIME) are
required to be used before continuing on with the R code.
