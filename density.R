########################################
#### Ecology of Phorcus turbinatus  ####
#### Analysis of density data       ####
#### Luka Seamus Wright             ####
########################################

#### 1.  Data preparation ####
#### 1.1 Density data ####
d <- read.csv("~/Desktop/Projects/Phorcus/Data/density.csv")
pred <- d[d$species %in% c("Stramonita haemastoma", "Thalassoma pavo"),]
rownames(pred) <- NULL
d <- d[d$species == "Phorcus turbinatus",]
rownames(d) <- NULL

den <- d$adjusted
den.site <- d$site

#### 1.2 Distribution data ####
dis <- read.csv("~/Desktop/Projects/Phorcus/Data/distribution.csv")
dis <- dis[dis$species == "Phorcus turbinatus",]
rownames(dis) <- NULL

pos <- dis$lat.position
dis.site <- dis$site

#### 2.  Data analysis ####
#### 2.1 Density data ####
m1 <- lm(den ~ den.site)
par(mfrow = c(1,2))
plot(resid(m1) ~ fitted(m1)) # heterogeneity: there is a clear difference in residual spread
abline(0, 0)
boxplot(resid(m1) ~ den.site) # heterogeneity: Ras has higher residual variance

hist(resid(m1))
qqnorm(resid(m1))
qqline(resid(m1)) # slightly right-skewed
# perhaps a gamma generalised linear model will improve the fit (Poisson would be best
# if the count data were left as integers but they were changed to numeric densities)

m2 <- glm(den+1 ~ den.site, family = Gamma(link = "log")) # +1 because model cannot take zeros
plot(resid(m2) ~ fitted(m2)) # clear improvement: homogenous
abline(0, 0)
boxplot(resid(m2) ~ den.site) # all sites now have similar residual variance

hist(resid(m2))
qqnorm(resid(m2))
qqline(resid(m2)) # near perfect normality
# m2 is chosen as optimal

require(car)
Anova(m2, type = 2)
# Response: den + 1
#          LR Chisq Df Pr(>Chisq)    
# den.site   258.38  2  < 2.2e-16 ***

summary(m2)
# Dwejra vs. Ras, t = 7.412, p < 0.001
# Dwejra vs. Xwejni, t = 9.394, p < 0.001

den.site <- factor(den.site, levels = c("Xwejni", "Dwejra", "Ras"))
m2 <- glm(den+1 ~ den.site, family = Gamma(link = "log"))
summary(m2)
# Xwejni vs. Ras, t = 16.806, p < 0.001

#### 2.2 Distribution data ####
m3 <- lm(pos ~ dis.site)
par(mfrow = c(1,2))
plot(resid(m3) ~ fitted(m3)) 
abline(0, 0)
boxplot(resid(m3) ~ dis.site) # good homogeneity

hist(resid(m3))
qqnorm(resid(m3))
qqline(resid(m3)) # good normality
# m3 is chosen as optimal

Anova(m3, type = 2)
# Response: pos
#           Sum Sq  Df F value    Pr(>F)    
# dis.site  5.9062   2  253.62 < 2.2e-16 ***
# Residuals 8.3488 717  

summary(m3)
# Dwejra vs. Ras, t = 20.741, p < 0.001
# Dwejra vs. Xwejni, t = 0.782, p = 0.43

dis.site <- factor(dis.site, levels = c("Xwejni", "Dwejra", "Ras"))
m3 <- lm(pos ~ dis.site)
summary(m3)
# Xwejni vs. Ras, t = 18.649, p < 0.001

#### 3.  Data visualisation ####
#### 3.2 Density data ####
require(psych)
stat <- with(pred, describeBy(adjusted, species, mat = T))
stat

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

d$site <- factor(d$site, levels = c("Xwejni", "Dwejra", "Ras"))

require(ggridges)
plot <- ggplot(data = d) +
  geom_density_ridges(aes(adjusted, site),
                      scale = 1.5, bandwidth = 8)

require(dplyr)
density_lines <- d %>% 
  group_by(site) %>% 
  summarise(x_mean = mean(adjusted)) %>% 
  mutate(group = as.integer(site)) %>% 
  left_join(ggplot_build(plot) %>% purrr::pluck("data", 1),
            on = "group") %>% 
  group_by(group) %>%
  summarise(x_mean = first(x_mean), 
            density = approx(x, density, first(x_mean))$y, 
            scale = first(scale), 
            iscale = first(iscale))
# probability density at mean position per species

d$site <- as.character(d$site)
d <- d %>% 
  mutate(order = case_when(
    startsWith(site, "X") ~ "a",
    startsWith(site, "D") ~ "b",
    startsWith(site, "R") ~ "c",
  ))
d$site <- factor(d$site)
d$order <- factor(d$order)

require(ggnewscale)
denp <- ggplot() +
        geom_density_ridges_gradient(data = d[d$site == "Ras",],
                                     aes(x = adjusted, y = order, fill = factor(stat(quantile))),
                                     scale = 60, bandwidth = 8, calc_ecdf = TRUE, colour = "#a1d5cf",
                                     quantile_lines = TRUE, quantiles = c(0.025, 0.5, 0.975)) +
        scale_fill_manual(values = alpha("#a1d5cf", c(0.2, 0.5, 0.5, 0.2))) +
        new_scale_fill() +
        geom_density_ridges_gradient(data = d[d$site == "Dwejra",],
                                     aes(x = adjusted, y = order, fill = factor(stat(quantile))),
                                     scale = 60, bandwidth = 8, calc_ecdf = TRUE, colour = "#a1d5cf",
                                     quantile_lines = TRUE, quantiles = c(0.025, 0.5, 0.975)) +
        scale_fill_manual(values = alpha("#a1d5cf", c(0.2, 0.5, 0.5, 0.2))) +
        new_scale_fill() +
        geom_density_ridges_gradient(data = d[d$site == "Xwejni",],
                                     aes(x = adjusted, y = order, fill = factor(stat(quantile))),
                                     scale = 60, bandwidth = 8, calc_ecdf = TRUE, colour = "#a1d5cf",
                                     quantile_lines = TRUE, quantiles = c(0.025, 0.5, 0.975)) +
        scale_fill_manual(values = alpha("#a1d5cf", c(0.2, 0.5, 0.5, 0.2))) +
        geom_segment(data = density_lines,
                     aes(x = x_mean, y = group, xend = x_mean,
                         yend = group + density * scale * iscale,
                         colour = factor(group)), lty = 5, colour = "#a1d5cf") +
        annotate("text", x = 500, y = 1.15:3.15, size = 4.2, hjust = 0, angle = 270,
                 family = "Helvetica Neue", label = c("Il-Bajja tax-Xwejni",
                                                      "Il-Qala tad-Dwejra",
                                                      "Id-Dejjaq tar-Ras")) +
        annotate("text", x = 485, y = 1.5:2.5, size = 4.2, family = "Helvetica Neue",
                 label = c("***", "***")) +
        labs(x = expression("Density (m"^-1*")")) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_discrete(position = "right", expand = expansion(add = c(0, 1.1)),
                         labels = c(expression(italic("n ")*"= 100"),
                                    expression(italic("n ")*"= 100"),
                                    expression(italic("n ")*"= 100"))) +
        coord_flip(xlim = c(0, 500)) +
        mytheme +
        theme(axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.text.x = element_text(hjust = 0),
              axis.title.x = element_blank(),
              legend.position = "none")

denp

#### 3.2 Distribution data ####
dis$site <- factor(dis$site, levels = c("Xwejni", "Dwejra", "Ras"))

plot <- ggplot(data = dis) +
  geom_density_ridges(aes(lat.position, site),
                      scale = 1.5, bandwidth = 0.03)

density_lines <- dis %>% 
  group_by(site) %>% 
  summarise(x_mean = mean(lat.position)) %>% 
  mutate(group = as.integer(site)) %>% 
  left_join(ggplot_build(plot) %>% purrr::pluck("data", 1),
            on = "group") %>% 
  group_by(group) %>%
  summarise(x_mean = first(x_mean), 
            density = approx(x, density, first(x_mean))$y, 
            scale = first(scale), 
            iscale = first(iscale))
# probability density at mean position per species

dis$site <- as.character(dis$site)
dis <- dis %>% 
  mutate(order = case_when(
    startsWith(site, "X") ~ "a",
    startsWith(site, "D") ~ "b",
    startsWith(site, "R") ~ "c",
  ))
dis$site <- factor(dis$site)
dis$order <- factor(dis$order)

disp <- ggplot() +
        geom_density_ridges_gradient(data = dis[dis$site == "Ras",],
                                     aes(x = lat.position, y = order, fill = factor(stat(quantile))),
                                     scale = 0.285, bandwidth = 0.03, calc_ecdf = TRUE, colour = "#a1d5cf",
                                     quantile_lines = TRUE, quantiles = c(0.025, 0.5, 0.975)) +
        scale_fill_manual(values = alpha("#a1d5cf", c(0.2, 0.5, 0.5, 0.2))) +
        new_scale_fill() +
        geom_density_ridges_gradient(data = dis[dis$site == "Dwejra",],
                                     aes(x = lat.position, y = order, fill = factor(stat(quantile))),
                                     scale = 0.285, bandwidth = 0.03, calc_ecdf = TRUE, colour = "#a1d5cf",
                                     quantile_lines = TRUE, quantiles = c(0.025, 0.5, 0.975)) +
        scale_fill_manual(values = alpha("#a1d5cf", c(0.2, 0.5, 0.5, 0.2))) +
        new_scale_fill() +
        geom_density_ridges_gradient(data = dis[dis$site == "Xwejni",],
                                     aes(x = lat.position, y = order, fill = factor(stat(quantile))),
                                     scale = 0.285, bandwidth = 0.03, calc_ecdf = TRUE, colour = "#a1d5cf",
                                     quantile_lines = TRUE, quantiles = c(0.025, 0.5, 0.975)) +
        scale_fill_manual(values = alpha("#a1d5cf", c(0.2, 0.5, 0.5, 0.2))) +
        geom_segment(data = density_lines,
                     aes(x = x_mean, y = group, xend = x_mean,
                         yend = group + density * scale * iscale,
                         colour = factor(group)), lty = 5, colour = "#a1d5cf") +
        annotate("text", x = 0.78, y = 1.5:2.5, size = 4.2, family = "Helvetica Neue",
                 label = c("n.s.", "***")) +
        labs(x = expression("Position above LAT (m)")) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_discrete(position = "right", expand = expansion(add = c(0, 1.2)),
                         labels = c(expression(italic("n ")*"= 217"),
                                    expression(italic("n ")*"= 293"),
                                    expression(italic("n ")*"= 210"))) +
        coord_flip(xlim = c(0, 0.8)) +
        mytheme +
        theme(axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.text.x = element_text(hjust = 0),
              axis.title.x = element_blank(),
              legend.position = "none")

disp

#### 3.3 Combined plot ####
require(cowplot)
finalp <- plot_grid(denp, disp, labels = "auto", label_size = 15, nrow = 2)
finalp # 6 x 8 in

#### 4.  Cleanup ####
detach(package:car)
detach(package:psych)
detach(package:cowplot)
detach(package:ggnewscale)
detach(package:ggridges)
detach(package:ggplot2)
detach(package:dplyr)
rm(list = ls())
graphics.off()
cat("\014")
