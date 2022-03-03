#####################################################
#### Ecology of Phorcus turbinatus               ####
#### Analysis of predation and distribution data ####
#### Luka Seamus Wright                          ####
#####################################################

#### 1.  Data preparation ####
#### 1.1 Predation data ####
p <- read.csv("~/Desktop/Projects/Phorcus/Data/predation.csv")
ppos <- p$lat.position
mass <- p$mass
pred <- p$predation

#### 1.2 Distribution data ####
d <- read.csv("~/Desktop/Projects/Phorcus/Data/distribution.csv")
d$species <- factor(d$species, levels = c("Phorcus turbinatus",
                                          "Stramonita haemastoma",
                                          "Thalassoma pavo",
                                          "Hermodice carunculata",
                                          "Hexaplex trunculus"))
dpos <- d$lat.position
sp <- d$species

#### 2.  Data analysis ####
#### 2.1 Predation data ####
m1 <- glm(pred ~ ppos * mass, family = binomial(link = "logit"))
m2 <- glm(pred ~ ppos * mass, family = binomial(link = "cauchit"))
AIC(m1, m2) # m2 (model with cauchit link function) fits better

drop1(m2, test = "Chisq") # no interaction
m2 <- glm(pred ~ ppos + mass, family = binomial(link = "cauchit"))

par(mfrow = c(1,2))
plot(resid(m2)  ~ fitted(m2))
abline(0, 0)
qqnorm(resid(m2))
qqline(resid(m2))
par(mfrow = c(1,1))
# a few outliers where individuals sporadically survived at depth
# but the model fit generally looks good

require(car)
Anova(m2, type = 2)
# Response: pred
#      LR Chisq Df Pr(>Chisq)    
# ppos  28.3706  1  1.002e-07
# mass   2.1643  1     0.1412  
  
summary(m2)
# position: z ratio = -2.68, p = 0.007
# mass: z ratio = -1.58, p = 0.11

coef(m2)
# a = 0.7105266
# b1 = -1.8982500
# b2 = -0.3371762
# The equation for the cumulative distribution function of the Cauchy distribution
# is y = 1/pi*atan(a + b*x) + 0.5 so in this case y = 1/pi*atan(0.71 - 1.9x1 - 0.34x2) + 0.5
# describes the relationship between position relative to lowest astronomical tide,
# gastropod mass and predation risk. Since I am interested in predicting the effect
# of position on predation risk while keeping mass constant, the equation is
# y = 1/pi*atan(0.71 - 1.9x1 - 0.34*mean(mass)) + 0.5

#### 2.2 Distribution data ####
m3 <- lm(dpos ~ sp)

par(mfrow = c(1,2))
plot(resid(m3)  ~ fitted(m3))
abline(0, 0)
plot(resid(m3)  ~ sp) # clear heterogeneity

hist(resid(m3))
qqnorm(resid(m3))
qqline(resid(m3))  # deviance from normality at lower end
par(mfrow = c(1,1))

# use nlme::gls() to model heterogeneity
require(nlme)
m4 <- gls(dpos ~ sp, weights = varIdent(form = ~1|sp))

par(mfrow = c(1,2))
plot(resid(m4, type = "normalized")  ~ fitted(m3, type = "normalized"))
abline(0, 0)
plot(resid(m4, type = "normalized")  ~ sp) # more homogenous
hist(resid(m4, type = "normalized"))
qqnorm(resid(m4, type = "normalized"))
qqline(resid(m4, type = "normalized"))  # almost perfectly normal
par(mfrow = c(1,1))

Anova(m4, type = 2)
# Response: dpos
#    Df  Chisq Pr(>Chisq)    
# sp  4 1991.8  < 2.2e-16

summary(m4)
# turbinatus vs. haemastoma, t = -13.32, p < 0.001
# turbinatus vs. pavo, t = -40.37, p < 0.001
# turbinatus vs. carunculata, t = -13.49, p < 0.001
# turbinatus vs. trunculus, t = -9.55, p < 0.001

sp <- factor(sp, levels = c("Stramonita haemastoma", "Thalassoma pavo", "Hermodice carunculata", 
                            "Hexaplex trunculus", "Phorcus turbinatus"))
m4 <- gls(dpos ~ sp, weights = varIdent(form = ~1|sp))
summary(m4)
# haemastoma vs. pavo, t = -23.46, p < 0.001
# haemastoma vs. carunculata, t = -11.31, p < 0.001
# haemastoma vs. trunculus, t = -8.41, p < 0.001

sp <- factor(sp, levels = c("Thalassoma pavo", "Hermodice carunculata", "Hexaplex trunculus", 
                            "Phorcus turbinatus", "Stramonita haemastoma"))
m4 <- gls(dpos ~ sp, weights = varIdent(form = ~1|sp))
summary(m4)
# pavo vs. carunculata, t = -5.82, p < 0.001
# pavo vs. trunculus, t = -5.45, p < 0.001

sp <- factor(sp, levels = c("Hermodice carunculata", "Hexaplex trunculus", "Phorcus turbinatus", 
                            "Stramonita haemastoma", "Thalassoma pavo"))
m4 <- gls(dpos ~ sp, weights = varIdent(form = ~1|sp))
summary(m4)
# carunculata vs. trunculus, t = -2.02, p = 0.04

sp <- factor(sp, levels = c("Phorcus turbinatus", "Stramonita haemastoma", "Thalassoma pavo", 
                            "Hermodice carunculata", "Hexaplex trunculus"))

#### 3.  Data visualisation ####
#### 3.1 Predation data ####
new <- data.frame(ppos = seq(min(ppos), 1, 0.01), # generate new data
                  mass = rep(mean(mass), 998))
inv <- family(m2)$linkinv # calculate inverse of link function
predicted <- predict(m2, type = "link", se.fit = TRUE, newdata = new) # predict

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

pp <- ggplot(data = new) +
        geom_rug(data = p, aes(ppos, colour = factor(pred, levels = c("1", "0")))) +
        scale_colour_manual(values = c("#cdd0d1", "#000000"),
                            labels = c("Yes", "No"),
                            guide = guide_legend(title = "Predation")) +
        geom_line(aes(ppos, fit), colour = "#a1d5cf") +
        geom_ribbon(aes(ppos, ymin = lower, ymax = upper), fill = "#a1d5cf", alpha = 0.5) +
        annotate("text", x = c(-5.1, -5.55), y = c(0.066, 0.063), size = 4.2, hjust = 0, parse = T,
                 label = c("italic('N ')*'= 100'", "bold(Equation)"), family = "Helvetica Neue") +
        geom_text(aes(-5.8, 0.066), label = expression("y = tan"^-1*"(–1.9x – 0.14)/"*pi*" + 0.5"),
                  family = "Helvetica Neue", size = 4.2, hjust = 0, check_overlap = T) +
        labs(y = expression("Predation risk (d"^-1*")"),
             x = "Position in relation to LAT (m)") +
        scale_x_continuous(breaks = seq(-6, 1, by = 0.5)) +
        scale_y_continuous(breaks = seq(0, 1, by = 0.1), position = "right") +
        coord_flip(xlim = c(-6, 1), ylim = c(0, 1), expand = FALSE) +
        theme(legend.position = c(0.18, 0.2)) +
        mytheme

pp

#### 3.2 Distribution data ####
n <- aggregate(dpos ~ sp, length, data = d)
colnames(n) <- c("species", "n")
n # sample size per species

require(hdrcde)
dhdr <- aggregate(dpos ~ sp, function(x){hdr(x, prob = 85, h = 0.15)}, data = d)
dhdr <- with(dhdr, data.frame(sp = sp,
                              l1 = c(dpos[[1]][,1], dpos[[2]][,1], dpos[[3]][,1], 
                                     dpos[[4]][,1], dpos[[5]][,3]),
                              u1 = c(dpos[[1]][,2], dpos[[2]][,2], dpos[[3]][,2], 
                                     dpos[[4]][,2], dpos[[5]][,4]),
                              l2 = c(NA, NA, NA, NA, dpos[[5]][,1]),
                              u2 = c(NA, NA, NA, NA, dpos[[5]][,2])))
dhdr # 85% probability highest density regions per species

require(ggridges)
plot <- ggplot(data = d) +
          geom_density_ridges(aes(x = dpos, y = sp),
                              scale = 1.5, bandwidth = 0.15)

require(dplyr)
density_lines <- d %>% 
  group_by(species) %>% 
  summarise(x_mean = mean(lat.position)) %>% 
  mutate(group = as.integer(species)) %>% 
  left_join(ggplot_build(plot) %>% purrr::pluck("data", 1),
            on = "group") %>% 
  group_by(group) %>%
  summarise(x_mean = first(x_mean), 
            density = approx(x, density, first(x_mean))$y, 
            scale = first(scale), 
            iscale = first(iscale))
# probability density at mean position per species

d$species <- as.character(d$species)
d <- d %>%
  mutate(order = case_when(
    startsWith(species, "P") ~ "a",
    startsWith(species, "S") ~ "b",
    startsWith(species, "T") ~ "c",
    startsWith(species, "Her") ~ "d",
    startsWith(species, "Hex") ~ "e"
  ))
d$species <- factor(d$species)
d$order <- factor(d$order)

require(ggnewscale)
dp <- ggplot() +
        geom_vline(aes(xintercept = 0.23)) +
        geom_density_ridges_gradient(data = d[d$species == "Hexaplex trunculus",],
                                     aes(x = lat.position, y = order, fill = factor(stat(quantile))),
                                     scale = 0.768, bandwidth = 0.15, calc_ecdf = TRUE, colour = "#b194c0",
                                     quantile_lines = TRUE, quantiles = c(0.025, 0.5, 0.975)) +
        scale_fill_manual(values = alpha("#b194c0", c(0.2, 0.5, 0.5, 0.2))) +
        new_scale_fill() +
        geom_density_ridges_gradient(data = d[d$species == "Hermodice carunculata",],
                                     aes(x = lat.position, y = order, fill = factor(stat(quantile))),
                                     scale = 0.768, bandwidth = 0.15, calc_ecdf = TRUE, colour = "#b3061e",
                                     quantile_lines = TRUE, quantiles = c(0.025, 0.5, 0.975)) +
        scale_fill_manual(values = alpha("#b3061e", c(0.2, 0.5, 0.5, 0.2))) +
        new_scale_fill() +
        geom_density_ridges_gradient(data = d[d$species == "Thalassoma pavo",],
                                     aes(x = lat.position, y = order, fill = factor(stat(quantile))),
                                     scale = 0.768, bandwidth = 0.15, calc_ecdf = TRUE, colour = "#44b3df",
                                     quantile_lines = TRUE, quantiles = c(0.025, 0.5, 0.975)) +
        scale_fill_manual(values = alpha("#44b3df", c(0.2, 0.5, 0.5, 0.2))) +
        new_scale_fill() +
        geom_density_ridges_gradient(data = d[d$species == "Stramonita haemastoma",],
                                     aes(x = lat.position, y = order, fill = factor(stat(quantile))),
                                     scale = 0.768, bandwidth = 0.15, calc_ecdf = TRUE, colour = "#f5a54a",
                                     quantile_lines = TRUE, quantiles = c(0.025, 0.5, 0.975)) +
        scale_fill_manual(values = alpha("#f5a54a", c(0.2, 0.5, 0.5, 0.2))) +
        new_scale_fill() +
        geom_density_ridges_gradient(data = d[d$species == "Phorcus turbinatus",],
                                     aes(x = lat.position, y = order, fill = factor(stat(quantile))),
                                     scale = 0.768, bandwidth = 0.15, calc_ecdf = TRUE, colour = "#a1d5cf",
                                     quantile_lines = TRUE, quantiles = c(0.025, 0.5, 0.975)) +
        scale_fill_manual(values = alpha("#a1d5cf", c(0.2, 0.5, 0.5, 0.2))) +
        geom_segment(data = density_lines,
                     aes(x = x_mean, y = group, xend = x_mean,
                         yend = group + density * scale * iscale,
                         colour = factor(group)), lty = 5) +
        scale_colour_manual(values = c("#a1d5cf", "#f5a54a", "#44b3df", "#b3061e", "#b194c0")) +
        annotate("text", x = -5.99, y = 1.2:5.2, size = 4.2, hjust = 1, angle = 270, parse = T,
                 family = "Helvetica Neue", label = c("italic('Phorcus turbinatus')",
                                                      "italic('Stramonita haemastoma')",
                                                      "italic('Thalassoma pavo')",
                                                      "italic('Hermodice carunculata')",
                                                      "italic('Hexaplex trunculus')")) +
        annotate("text", x = 0.92, y = 1.5:4.5, size = 4.2, family = "Helvetica Neue",
                 label = c("***", "***", "***", "*")) +
        labs(y = expression("Predation risk (d"^-1*")"),
             x = "Position in relation to LAT (m)") +
        scale_x_continuous(breaks = seq(-6, 1, by = 0.5), expand = c(0, 0)) +
        scale_y_discrete(position = "right", expand = expansion(add = c(0, 0.6)),
                         labels = c(expression(italic("n ")*"= 720"),
                                    expression(italic("n ")*"= 220"),
                                    expression(italic("n ")*"= 540"),
                                    expression(italic("n ")*"= 70"),
                                    expression(italic("n ")*"= 50"))) +
        coord_flip(xlim = c(-6, 1)) +
        mytheme +
        theme(axis.ticks.x = element_line(colour = FALSE),
              axis.line.x = element_line(colour = FALSE),
              axis.text.x = element_text(hjust = 0),
              axis.title.x = element_text(colour = FALSE),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              legend.position = "none")

dp

#### 3.3 Combined plot ####
require(cowplot)
finalp <- plot_grid(pp, dp, labels = "auto", label_size = 15)
finalp # 6 x 8 in

#### 4.  Cleanup ####
detach(package:car)
detach(package:cowplot)
detach(package:ggnewscale)
detach(package:ggridges)
detach(package:ggplot2)
detach(package:dplyr)
detach(package:hdrcde)
detach(package:nlme)
rm(list = ls())
graphics.off()
cat("\014")
