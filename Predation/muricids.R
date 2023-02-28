########################################
#### Ecology of Phorcus turbinatus  ####
#### Analysis of laboratory data    ####
#### Luka Seamus Wright             ####
########################################

#### 1.  Data preparation ####
mu <- read.csv("~/PATH/muricids.csv")
mu <- mu[1:34,] # remove Stramonita haemastoma
mass <- with(mu, prey.mass/predator.mass) # prey:predator mass ratio
level <- mu$level # aquarium water level
p <- mu$predation # predation risk

#### 2.  Data analysis ####
m1 <- glm(p ~ level * mass, family = binomial(link = "cauchit"))
m2 <- glm(p ~ level * mass, family = binomial(link = "cloglog"))
AIC(m1, m2) # m2 (model with cloglog link function) fits better

drop1(m2, test = "Chisq") # no interaction
m3 <- glm(p ~ level + mass, family = binomial(link = "cloglog"))

par(mfrow = c(1, 2))
plot(resid(m3) ~ fitted(m3))
abline(0, 0) # good homogeneity
qqnorm(resid(m3))
qqline(resid(m3)) # ok normality
par(mfrow = c(1, 1))

require(car)
Anova(m3, type = 2)
# Response: p
#       LR Chisq Df Pr(>Chisq)
# level  1.82464  1     0.1768
# mass   0.44881  1     0.5029

summary(m3)
# level, z = -1.302, p = 0.19
# mass, z = 0.742, p = 0.46

coef(m3)
# a = -0.4918968
# b1 = -6.0522056
# b2 = 0.7739422
# The cloglog equation is derived from the cumulative distribution function of 
# the Gumbel distribution and expressed as y = 1-exp(-exp(a + b1*x1 + b2*x2)),
# so in this case y = 1-exp(-exp(-0.49 - 6.05*x1 + 0.77*x2)) describes the 
# relationship between predation probability, aquarium water level and
# and prey:predator mass ratio. Since I am interested in predicting the effect 
# of water leevel on predation while keeping the mass ratio constant, the equation is
# y = 1-exp(-exp(-0.49 - 6.05*x1 + 0.77*mean(mass)))

#### 3.  Data visualisation ####
mean(p) # 35.29 % average predation risk
mean(mass) # average prey:predator mass ratio = 0.5065038

stat <- aggregate(p ~ level, length, data = mu) # calculate sample sizes
colnames(stat)[2] <- "n"
stat$ones <- aggregate(p ~ level, function(x) sum(x == 1), data = mu)$p # count ones

require(binom)
se <- with(stat, 
           binom.confint(ones, n, conf.level = (pnorm(1)-0.5)*2,     # calculate Wilson
                         method = "wilson", type = "central"))[,3:6] # standard errors
se$level <- stat$level

new <- data.frame(level = seq(min(level), max(level), 0.01), # generate new data
                  mass = rep(mean(mass), 21))
inv <- family(m3)$linkinv # calculate inverse of link function
predicted <- predict(m3, type = "link", se.fit = TRUE, newdata = new) # predict

new$fit <- inv(predicted$fit)
new$upper <- inv(predicted$fit + predicted$se.fit * qnorm(0.975))
new$lower <- inv(predicted$fit - predicted$se.fit * qnorm(0.975))

require(ggplot2)
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .3, .2, .2),"cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 15, hjust = 0),
                 axis.text = element_text(size = 12, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.key.size = unit(.3, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.spacing.x = unit(.1, "cm"),
                 legend.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title.align = 0,
                 legend.title = element_text(size = 12, face = "bold"),
                 text = element_text(family = "Helvetica Neue"))

mp <- ggplot(data = new) +
        geom_line(aes(level, fit), colour = "#a1d5cf", lty = 5) +
        geom_ribbon(aes(level, ymin = lower, ymax = upper), fill = "#a1d5cf", alpha = 0.5) +
        geom_linerange(data = se, aes(level, ymin = lower, ymax = upper),
                       size = 0.5, colour = "#a1d5cf") +
        geom_point(data = se, aes(level, mean, size = n), colour = "#a1d5cf") +
        scale_size_continuous(labels = c(expression(italic("n ")*"= 4"),
                                         expression(italic("n ")*"= 5"),
                                         expression(italic("n ")*"= 6"),
                                         expression(italic("n ")*"= 7"),
                                         expression(italic("n ")*"= 8"))) +
        geom_rug(aes(0.25, mean(p)), colour = "#a1d5cf", sides = "r",
                 length = unit(.25, "cm")) +
        annotate("text", x = c(0.155, 0.025), y = c(0.801, 0.03), size = 4.2, hjust = 0, parse = T,
                 label = c("italic('N ')*'= 34'", "bold(Equation)"), family = "Helvetica Neue") +
        geom_text(aes(0.01, 0.033), label = expression(italic("y ")*"= 1 – e"^-e^{"–6.05"*italic("x ")*"– 0.1"}),
                  family = "Helvetica Neue", size = 4.2, hjust = 0, check_overlap = T) +
        labs(y = expression("Predation risk (d"^-1*")"),
             x = "Aquarium water level (m)") +
        coord_flip(xlim = c(0, 0.25), ylim = c(0, 1), expand = FALSE) +
        scale_y_continuous(breaks = seq(0, 1, by = 0.1), position = "right") +
        scale_x_continuous(breaks = seq(0, 0.25, by = 0.05)) +
        mytheme +
        theme(legend.position = c(0.9, 0.82),
              legend.title = element_blank())
        

mp # dimensions: 4 x 4 in

#### 4.  Cleanup ####
detach(package:car)
detach(package:ggplot2)
detach(package:binom)
rm(list = ls())
graphics.off()
cat("\014")
