########################################
#### Ecology of Phorcus turbinatus  ####
#### Analysis of physiological data ####
#### Luka Seamus Wright             ####
########################################

#### 1.  Data preparation ####
#### 1.1 Length:mass data ####
l.m <- read.csv("~/Desktop/Projects/Phorcus/Data/length.mass.csv")

require(ggplot2)
ggplot(l.m, aes(length, mass)) + 
  geom_point() + 
  facet_wrap(~species)
# there is one clear outlier in the H. trunculus data that should be removed
l.m[20,] <- NA # remove outlier
l.m <- l.m[complete.cases(l.m),] # remove NA
rownames(l.m) <- NULL # reset rownames
ggplot(l.m, aes(length, mass)) + 
  geom_point() + 
  facet_wrap(~species)
# the trend for P. turbinatus is clearly nonlinear
# but a power curve may fit all trends best

l.ml <- l.m$length
l.mm <- l.m$mass
sp <- factor(l.m$species)

#### 1.2 Physiology data ####
p <- read.csv("~/Desktop/Projects/Phorcus/Data/physiology.csv")
pl <- p$length
pm <- p$mass
d <- p$day
s <- p$survival

#### 2.  Data analysis ####
#### 2.1 Length:mass data ####
m1 <- lm(l.mm ~ l.ml * sp)
par(mfrow = c(1,2))
plot(resid(m1) ~ fitted(m1)) # heteroskedasticity: there is a clear nonlinear trend
abline(0, 0)
boxplot(resid(m1) ~ sp) # heterogeneity: difference in residual variance between species

hist(resid(m1))
qqnorm(resid(m1))
qqline(resid(m1)) # not normal but balanced
# perhaps a gamma generalised linear model will improve the fit

m2 <- glm(l.mm ~ l.ml * sp, family = Gamma(link = "log"))
plot(resid(m2) ~ fitted(m2)) # not much improvement: nonlinear trend still visible
abline(0, 0)
boxplot(resid(m2) ~ sp) # slightly better homogeneity

hist(resid(m2))
qqnorm(resid(m2))
qqline(resid(m2)) # better normality
# nonlinear (power) modelling is required

m3 <- nls(l.mm ~ k[sp]*l.ml^n[sp], 
          start = list(k = c(0.0002, 0.0011, 0.0002), # starting values were obtained from
                       n = c(2.8426, 2.7571, 2.9864)))# standard curve fitting (e.g. Excel)

plot(resid(m3) ~ fitted(m3))
abline(0, 0) # obvious trend in residuals is gone but strong heteroscedasticity
boxplot(resid(m3) ~ sp) # strong heterogeneity

hist(resid(m3))
qqnorm(resid(m3))
qqline(resid(m3)) # not normal but balanced
# try modelling heteroscedasticity with gnls()

require(nlme)
m4 <- gnls(mass ~ k*length^n, 
           start = list(k = c(0.0002, 0.0011, 0.0002), 
                        n = c(2.8426, 2.7571, 2.9864)),
           params = list(k ~ species, n ~ species),
           weights = varPower(),
           data = l.m)

plot(resid(m4,  type = "normalized") ~ fitted(m4,  type = "normalized")) 
abline(0, 0)
boxplot(resid(m4,  type = "normalized") ~ sp) # almost perfectly homogenous

hist(resid(m4,  type = "normalized"))
qqnorm(resid(m4,  type = "normalized"))
qqline(resid(m4,  type = "normalized")) # perfectly normal
par(mfrow = c(1,1))
# m4 is chosen as the optimal model

coef(m4)
# P. turbinatus equation: y = 0.001074418x^2.752005
# S. haemastoma equation: y = 0.0002028907x^2.978715
# H. trunculus equation: y = 0.0002296665x^2.848017

l.m$species <- factor(l.m$species, levels = c("Phorcus turbinatus",
                                              "Hexaplex trunculus",
                                              "Stramonita haemastoma"))
m4 <- gnls(mass ~ k*length^n, 
           start = list(k = c(0.0011, 0.0002, 0.0002), 
                        n = c(2.7571, 2.8426, 2.9864)),
           params = list(k ~ species, n ~ species),
           weights = varPower(),
           data = l.m)
summary(m4)
# P. turbinatus, k: t = 15.36255, p < 0.001
#                n: t = 116.16905, p < 0.001

l.m$species <- factor(l.m$species, levels = c("Hexaplex trunculus",
                                              "Phorcus turbinatus",
                                              "Stramonita haemastoma"))
m4 <- gnls(mass ~ k*length^n, 
           start = list(k = c(0.0002, 0.0011, 0.0002), 
                        n = c(2.8426, 2.7571, 2.9864)),
           params = list(k ~ species, n ~ species),
           weights = varPower(),
           data = l.m)

#### 2.2 Physiology data ####
pm[is.na(pm)] <- 0.001074418*pl[is.na(pm)]^2.752005 # replace all NAs with predicted values

m5 <- glm(s ~ d * pm, family = binomial(link = "logit"))
m6 <- glm(s ~ d * pm, family = binomial(link = "probit"))
AIC(m5, m6) # m6 (model with probit link function) fits better

drop1(m6, test = "Chisq") # no interaction
m6 <- glm(s ~ d + pm, family = binomial(link = "probit"))

par(mfrow = c(1, 2))
plot(resid(m6) ~ fitted(m6))
abline(0, 0) # good homogeneity
qqnorm(resid(m6))
qqline(resid(m6)) # good normality
par(mfrow = c(1, 1))

require(car)
Anova(m6, type = 2)
# Response: s
#    LR Chisq Df Pr(>Chisq)    
# d   135.214  1     <2e-16 ***
# pm    0.768  1     0.3809 

summary(m6)
# d, z = -8.341, p < 0.001
# pm, z = -0.916, p = 0.36

coef(m6)
# a = 3.19033764
# b1 = -0.75703629
# b2 = -0.09164321
# The equation for the cumulative distribution function of the Gaussian distribution
# is y = pnorm(a + b1*x1 + b2*x2) so in this case y = pnorm(3.19 + -0.76*x1 -0.09*x2)
# describes the relationship between survival probability, time spent immersed (d) 
# and gastropod mass. Since I am interested in predicting the effect of time on 
# survival while keeping mass constant, the equation is
# y = pnorm(3.19 + -0.76*x -0.09*mean(pm))

#### 3.  Data visualisation ####
#### 3.1 Length:mass data ####
aggregate(l.ml ~ species, min, data = l.m)
aggregate(l.ml ~ species, max, data = l.m)

new <- with(l.m, data.frame(length = c(seq(32.45, 55.43, 0.01), # generate new data
                                       seq(6.55, 26.90, 0.01),
                                       seq(34.07, 46.38, 0.01)),
                            species = c(rep("Hexaplex trunculus", 2299),
                                        rep("Phorcus turbinatus", 2036),
                                        rep("Stramonita haemastoma", 1232))))
new$fit <- predict(m4, newdata = new)

# bootstrap confidence interval
bootfun <- function(newdata) {
            start <- coef(m4)
            boot <- l.m[sample(nrow(l.m), size = nrow(l.m), replace = TRUE),]
            bootfit <- try(update(m4,
                        start = start,
                        data = boot),
                        silent = TRUE)
            if (inherits(bootfit, "try-error")) return(rep(NA, nrow(newdata)))
            predict(bootfit, newdata)
}

bmat <- replicate(1000, bootfun(new))
new$lwr <- apply(bmat, 1, quantile, 0.025, na.rm = TRUE)
new$upr <- apply(bmat, 1, quantile, 0.975, na.rm = TRUE)

new$species <- factor(new$species, levels = c("Phorcus turbinatus",
                                              "Stramonita haemastoma",
                                              "Hexaplex trunculus"))
l.m$species <- factor(l.m$species, levels = c("Phorcus turbinatus",
                                              "Stramonita haemastoma",
                                              "Hexaplex trunculus"))

annotation <- aggregate(l.ml ~ species, length, data = l.m)
colnames(annotation)[2] <- "n"
annotation$n <- c(expression(italic("n ")*"= 291"),
                  expression(italic("n ")*"= 40"),
                  expression(italic("n ")*"= 41"))
annotation$n <- as.character(annotation$n)
annotation$equation <- c(expression("y = 0.001x"^2.75),
                         expression("y = 0.0002x"^2.98),
                         expression("y = 0.0002x"^2.85))
annotation$equation <- as.character(annotation$equation)

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

l.mp <- ggplot(data = new) +
          geom_line(aes(length, fit, colour = species)) +
          geom_ribbon(aes(length, ymin = lwr, ymax = upr, fill = species), alpha = 0.5) +
          geom_point(data = l.m, aes(length, mass, colour = species), 
                     size = 1, shape = 16, alpha = 0.3) +
          scale_colour_manual(values = c("#a1d5cf", "#f5a54a", "#b194c0"),
                              guide = "none") +
          scale_fill_manual(values = c("#a1d5cf", "#f5a54a", "#b194c0"),
                            guide = "none") +
          facet_wrap(~species) +
          geom_text(data = annotation, aes(1.6, 24, label = n), 
                    family = "Helvetica Neue", parse = T, size = 4.2, hjust = 0) +
          geom_text(data = annotation, aes(1.6, 22, label = equation), 
                    family = "Helvetica Neue", parse = T, size = 4.2, hjust = 0) +
          labs(y = "Mass (g)",
               x = "Shell length (mm)") +
          scale_x_continuous(breaks = seq(0, 60, by = 20)) +
          coord_cartesian(xlim = c(0, 60), ylim = c(0, 25), expand = FALSE) +
          mytheme +
          theme(strip.background = element_blank(),
                strip.text = element_text(size = 12, hjust = 0,
                                          face = "italic"),
                panel.spacing = unit(.5, "cm"))

l.mp # dimensions: 4 x 8 in


#### 3.2 Physiology data ####
stat <- aggregate(s ~ d, length, data = p) # calculate sample sizes
colnames(stat)[2] <- "n"
stat$ones <- aggregate(s ~ d, function(x) sum(x == 1), data = p)$s # count ones

require(binom)
se <- with(stat, 
           binom.confint(ones, n, conf.level = (pnorm(1)-0.5)*2,     # calculate Wilson
                         method = "wilson", type = "central"))[,3:6] # standard errors
se$d <- stat$d

new <- data.frame(d = seq(0, max(d), 0.01), # generate new data
                  pm = rep(mean(pm), 701))
inv <- family(m6)$linkinv # calculate inverse of link function
predicted <- predict(m6, type = "link", se.fit = TRUE, newdata = new) # predict

new$fit <- inv(predicted$fit)
new$upper <- inv(predicted$fit + predicted$se.fit * qnorm(0.975))
new$lower <- inv(predicted$fit - predicted$se.fit * qnorm(0.975))

pp <- ggplot(data = new) +
        geom_line(aes(d, fit), colour = "#a1d5cf") +
        geom_ribbon(aes(d, ymin = lower, ymax = upper), fill = "#a1d5cf", alpha = 0.5) +
        geom_linerange(data = se, aes(d, ymin = lower, ymax = upper),
                        size = 0.5, colour = "#a1d5cf") +
        geom_point(data = se, aes(d, mean, size = factor(n)), colour = "#a1d5cf") +
        scale_size_manual(values = c(1.5, 3),
                          labels = c(expression(italic("n ")*"= 10"),
                                     expression(italic("n ")*"= 20"))) +
        annotate("text", x = c(0.17, 0.15), y = c(0.21, 0.105), size = 4.2, hjust = 0, parse = T,
                 label = c("italic('N ')*'= 200'", "bold(Equation)"), family = "Helvetica Neue") +
        geom_text(aes(0.172, 0.05), label = expression("y = "*phi*"(–0.76x + 3)"),
                  family = "Helvetica Neue", size = 4.2, hjust = 0, check_overlap = T) +
        labs(y = "Survival probability",
             x = "Time spent immersed (d)") +
        scale_x_continuous(breaks = seq(0, 7, by = 1)) +
        scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
        coord_cartesian(xlim = c(0, 7), ylim = c(0, 1), expand = FALSE) +
        mytheme +
        theme(legend.position = c(0.121, 0.3),
              legend.title = element_blank())
        

pp # dimensions: 4 x 4 in


#### 4.  Cleanup ####
detach(package:car)
detach(package:cowplot)
detach(package:ggplot2)
detach(package:nlme)
detach(package:binom)
rm(list = ls())
graphics.off()
cat("\014")
