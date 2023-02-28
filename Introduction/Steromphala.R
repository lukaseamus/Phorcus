################################################
#### Ecology of Phorcus turbinatus          ####
#### Calculation of introductory statistics ####
#### Luka Seamus Wright                     ####
################################################

S <- read.csv("~/PATH/Steromphala.csv")

require(psych)
Lstats <- with(S, describeBy(lower, species, mat = T))
Ustats <- with(S, describeBy(upper, species, mat = T))

detach(package:psych)
rm(list = ls())
graphics.off()
cat("\014")
