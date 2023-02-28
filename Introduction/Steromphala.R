
S <- read.csv("~/PATH/Steromphala.csv")

require(psych)
Lstats <- with(S, describeBy(Lower, Species, mat = T))
Ustats <- with(S, describeBy(Upper, Species, mat = T))

detach(package:psych)
rm(list = ls())
graphics.off()
cat("\014")
