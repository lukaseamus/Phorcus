S <- read.csv("~/Desktop/Projects/Phorcus/Data/Steromphala.csv")

require(psych)
Lstats <- with(S, describeBy(Lower, Species, mat = T))
Ustats <- with(S, describeBy(Upper, Species, mat = T))

detach(package:psych)