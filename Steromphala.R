###############################################################
#### Ecology of Phorcus turbinatus                         ####
#### Steromphala umbilicalis and S. cineraria distribution ####
#### Luka Seamus Wright                                    ####
###############################################################

S <- read.csv("~/PATH/Steromphala.csv")

require(psych)
Lstats <- with(S, describeBy(Lower, Species, mat = T))
Ustats <- with(S, describeBy(Upper, Species, mat = T))
