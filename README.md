## Habitat preference of an herbivore shapes the habitat distribution of its host plant
Alexandre, N. M., P. T. Humphrey, A. D. Gloss, J. Lee, J. Frazier, H. A. Affeldt III, and N. K. Whiteman. 2018.  Ecosphere. 10.1002/ ecs2.2372 (full citation pending).

### Files & Descriptions

#### Analysis R code

R scripts and associated narrative description of analysis can be found for all experiments presented in the main text in the files below. We provide Rmarkdown (`.Rmd`) files and corresponding compiled `.html` output.


1. Herbivory survey analysis: `01_Herbivory_survey_v1.html`, `01_Herbivory_survey_v1.Rmd`
2. Host Choice I: Sun- vs. shade-derived bittercress: `02_Host_source_choice_v1.html`, `02_Host_source_choice_v1.Rmd`
3. Host Choice II: Effects of light and temperature: `03_Habitat_type_choice_v1.html`, `03_Habitat_type_choice_v1.Rmd`

The file `_sun_shade_header.R` contains code necessary for all analysis code to be run. 

Data used in each of these analyses can be loaded from the `data` directory in the repo root. Users must set the working directory of R to the root of this cloned `github` repo in order to properly load the data.

The directory `model_output` is where output from the above analyses gets written. Most of the write steps have been commented out so users do not over-write the raw output used to generate Table 1 and Fig. 1 of the main text.

Finally, the directory `Appendix_S3` contains the `.tex` source, `.bib` BibTeX file, and associated figures required to fully compile the Statistical Supplement (Appendix S3) as published by the journal.

Please contact the owner of this repo (PT Humphrey; phumphrey [at] g.harvard.edu) with any questions.
