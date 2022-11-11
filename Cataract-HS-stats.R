library(readxl)
library(tidyverse)
library(DescTools)
library(lme4)
library(broom.mixed)
library(kableExtra)
library(xtable)
library(emmeans)
library(ggsci)
library(sjPlot)
library(coda)
library(rjags)
library(R2jags)
library(superdiag)
library(mcmcplots)
library(ggmcmc)
library(gridExtra)

cats <- read_excel(path = "C:/Users/edmondsonef/Desktop/Cataract/Hess Stats Class/GRSD.cataract.xlsx", sheet = "GRSD.cat")

# remove spaces from column and value names
names(cats) <- str_replace_all(names(cats), " ", "_")
cats <- cats %>%
  mutate(CoatColor = str_replace_all(coat_color, " ", "_"))

# turn categorical vars into factor
cats <- cats %>%
  rename(Age = `days`,
         Weight = Weight,
         Animal = original) %>%
  mutate(Sex = as.factor(sex),
         CoatColor = as.factor(CoatColor),
         Family = as.factor(family), # should this stay a factor? Yes?
         BCS = as.ordered(BCS),
         Treatment = relevel(as.factor(groups), ref = "Unirradiated"),
         #MyeloidLeukemia = as.factor(Myeloid_Leukemia),
         HarderianTumor = as.factor(Harderian_Tumor),
         #PreTLymphoma = as.factor(PreT_Lymphoma),
         Score = as.ordered(Cataract_Score)) # ordinal cat; leave as numeric?

# select vars, add binary conversion for score
cats <- cats %>%
  select(c(Animal, Sex, Weight, CoatColor, Family, BCS, Age, Treatment,Score)) %>%
  mutate(Cataracts = ifelse(Score < 2, 0, 1))

# glance at the dataset
# str(cats)
gsc <- cats %>% group_by(Score) %>%
  count(Treatment) %>%
  pivot_wider(names_from = Score, values_from = n) %>%
  group_by(Treatment) %>%
  mutate(Total = sum(c(`1`, `2`, `3`, `4`)))
xtable(gsc, label = "tab:gtab", caption = "Counts of score by treatment group") %>%
  xtable2kable(booktabs = T,
               include.rownames = FALSE,
               table.placement = NULL) %>%
  add_header_above(c( " " = 1, "Cataract Score" = 4, " " = 1)) %>%
  kable_styling(full_width = F, position = "float_right")
# barplot of proportion of each group with cataracts
sex_trt <- cats %>%
  group_by(Treatment, Sex, Cataracts) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(Treatment, Sex) %>%
  mutate(proportion = round(n/sum(n), digits = 2))
cats_grp <- sex_trt %>%
  filter(Cataracts == 1)
b1 <- ggplot(cats_grp, aes(x = Sex, y = proportion, fill = Treatment)) +
  geom_col(position = "dodge") +
  scale_fill_startrek() +
  theme_light()
b1
# Plot of average score of sex-group-family
gr_score <- cats %>%
  group_by(Treatment, Family) %>%
  summarize(mean_score = mean(Cataracts))
l1 <- ggplot(gr_score, aes(x = Treatment, y = mean_score, color = Family)) +
  geom_line(aes(group = Family)) +
  scale_x_discrete(expand = c(0, .2)) +
  theme_light() +
  theme(legend.position="none") +
  labs(y = "mean score")
l1
mod <- glmer(Cataracts ~ Treatment*Sex + (1|Family), data = cats, family = binomial)
# Specify the model
cat("model{
for(i in 1:N){
CAT[i] ~ dbern(p[i]) # Bernoulli-distributed response
logit(p[i])<- b0+ b1*Gamma[i] + b2*HZE[i] + b3*Male[i] +
b4*Male[i]*Gamma[i] + b5*Male[i]*HZE[i] + a[Family[i]] # likelihood function 
}
for(j in 1:nFam){
a[j] ~ dnorm(0, tau)
}
b0 ~ dnorm(0.0, 1.0E-3) # vaguely informative priors
b1 ~ dnorm(0.0, 1.0E-3)
b2 ~ dnorm(0.0, 1.0E-3)
b3 ~ dnorm(0.0, 1.0E-3)
b4 ~ dnorm(0.0, 1.0E-3)
b5 ~ dnorm(0.0, 1.0E-3)
tau ~ dgamma(1.0E-3,1.0E-3)
sigma2 <- 1/tau # convert precision 'tau' to variance 'sigma2'
}", file = "cat.jag")
# Prepare the data for JAGS
# break Treatment into dummy variables for each group
treatment <- model.matrix(~ Treatment - 1, cats)
sex <- model.matrix(~Sex -1, cats)
colnames(treatment) <- c("Unirradiated", "Gamma", "HZE")
cats <- data.frame(cats, treatment, sex)
# format relevant data as a list  
data <- list(CAT = cats$Cataracts, Gamma = cats$Gamma,
             HZE = cats$HZE, Male = cats$SexM, Family = cats$Family,
             nFam = length(unique(cats$Family)), N = nrow(cats))
nIter <- 60000
nChains <- 3
nThin <- 1
BurnIn <- 10000 
nAdapt <- 1000
ests <- summary(mod)$coef[,1] # pull starting values from frequentist model
var <- as.numeric(as.data.frame(VarCorr(mod))$vcov)
inits <- list(list("tau" = var+0.2, "b0" = ests[1]+0.5, "b1" = ests[2]+0.5, "b2" = ests[3]+0.5,
                   "b3" = ests[4]+0.2, "b4" = ests[5]+0.2, "b5" = ests[6]+0.2),
              list("tau" = var-0.2, "b0" = ests[1]-0.5, "b1" = ests[2]-0.5, "b2" = ests[3]-0.5,
                   "b3" = ests[4]-0.2, "b4" = ests[5]-0.2, "b5" = ests[6]-0.2),
              list("tau" = var, "b0" = ests[1], "b1" = ests[2], "b2" = ests[3],
                   "b3" = ests[4], "b4" = ests[5], "b5" = ests[6]))
# -- Compile and run the model
params <- c("b0", "b1", "b2", "b3", "b4", "b5", "sigma2")
set.seed(556)
model.fit <- jags(data = data,
                  inits = inits,
                  parameters.to.save = params,
                  model.file = "cat.jag",
                  n.chains = nChains,
                  n.iter = nIter,
                  n.burnin = BurnIn,
                  n.thin = nThin)
mcmc.model <- as.mcmc(model.fit)
#summary(mcmc.model)
  
posts <- mcmc.model[[3]][,-7]
#posts <- exp(posts)/(1 + exp(posts))
phpds <- HPDinterval(posts)
posts <- data.frame(posts)
# Create table of posterior estimates
means <- apply(posts, 2, mean)
medians <- apply(posts, 2, median)
mode_fun <- function(x) {
  ux <- round(unique(x), digits = 3)
  return(ux[which.max(tabulate(match(x, ux)))])
}
modes <- apply(round(posts, 4), 2, mode_fun)
sds <- apply(posts, 2, sd)
glmm_probs <- exp(ests) / (1 + exp(ests))
p_var <- exp(var) / (1 + exp(var))
est_tab <- round(data.frame(c(ests, var), means, medians, modes, sds, phpds), digits = 3)
rownames(est_tab) <- c("b_0", "b_1", "b_2", "b_3",
                       "b_4", "b_5", "sigma^2")
colnames(est_tab) <- c("GLMM Est", "MCMC Mean", "MCMC Median",
                       "MCMC Mode", "MCMC SD", "HPD Lower", "HPD Upper")
est_tab_show <- est_tab %>% select(-c(2:3))
xtable(est_tab_show, label = "tab:esttab",
       caption = "Final model parameter estimates on the log odds scale") %>%
  xtable2kable(booktabs = T, include.rownames = TRUE,
               table.placement = NULL, size="\\fontsize{9pt}{10pt}\\selectfont") %>%
  kable_styling(full_width = F, latex_options = "HOLD_position")
em1 <- emmeans(mod, ~ Treatment|Sex)
em1log <- regrid(em1, "log")
rrs1 <- contrast(em1log, interaction = "revpairwise", type = "response")
rrs1 <- as.data.frame(confint(rrs1)) %>%
  rename(Contrast = Treatment_revpairwise, Lower = asymp.LCL, Upper = asymp.UCL)
rrplot1 <- ggplot(rrs1, aes(x = ratio, y = Contrast, xmin = Lower, xmax = Upper)) +
  geom_errorbarh(aes(height = 0.2, color = Sex),
                 position = position_dodge(0.3), lwd = 1) +
  geom_point(aes(color = Sex), position = position_dodge(0.3)) +
  geom_text(aes(label = round(ratio, 2), color = Sex),
            position = position_dodge(0.7)) +
  theme_light() +
  geom_vline(aes(xintercept = 1), color = "#84BD00FF", lty = 2) +
  scale_y_discrete(labels = c("Gamma/Cont", "HZE/Cont","HZE/Gamma")) +
  theme(axis.text.y = element_text(angle=60)) +
  scale_color_startrek() +
  labs(x = "Relative Risk",
       title = "Within Sex")
em2 <- emmeans(mod, ~ Sex|Treatment)
em2log <- regrid(em2, "log")
rrs2 <- contrast(em2log, interaction = "revpairwise", type = "response")
rrs2 <- as.data.frame(confint(rrs2)) %>%
  rename(Contrast = Sex_revpairwise, Lower = asymp.LCL, Upper = asymp.UCL)
rrplot2 <- ggplot(rrs2, aes(x = ratio, y = Contrast, xmin = Lower, xmax = Upper)) +
  geom_errorbarh(aes(height = 0.2, color = Treatment),
                 position = position_dodge(0.3), lwd = 1) +
  geom_point(aes(color = Treatment), position = position_dodge(0.3)) +
  geom_text(aes(label = round(ratio, 2), color = Treatment),
            position = position_dodge(0.3), vjust = -1) +
  theme_light() +
  geom_vline(aes(xintercept = 1), color = "#FFCD00FF", lty = 2) +
  scale_color_startrek() +
  labs(x = "Relative Risk", y = "", title = "Within Treatment")
grid.arrange(rrplot1, rrplot2, ncol = 2)            
p_sigs <- est_tab[7,]
est <- as.numeric(p_sigs[1])
est <- exp(est) / (1 + exp(est))
psig <- as.numeric(p_sigs[4])
psig <- exp(psig) / (1 + exp(psig))
hpdl <- as.numeric(p_sigs[6])
hpdl <- exp(hpdl) / (1 + exp(hpdl))
hpdu <- as.numeric(p_sigs[7])
hpdu <- exp(hpdu) / (1 + exp(hpdu))
REs <- augment(ranef(mod,condVar = TRUE), ci.level = 0.95) %>%
  select(c(level, estimate, lb, ub)) %>%
  rename(Family = level) %>%
  mutate(Prob = exp(estimate)/(1+exp(estimate)),
         Lower = exp(lb)/(1+exp(lb)),
         Upper = exp(ub)/(1+exp(ub)))
colors <- c("Family Effect" = "#5C88DAFF", "GLMM Est" = "#CC0C00FF", "Bayes Mode" = "#84BD00FF",
            "HPD Interval" = "#FFCD00FF")
ggplot(REs, aes(x = Prob, y = Family, xmin = Lower, xmax = Upper)) +
  geom_errorbarh(aes(height = 0, color = "Family Effect")) +
  geom_point(aes(color = "Family Effect")) +
  geom_vline(aes(xintercept = est, color = "GLMM Est"), lwd = 1, lty = 2) +
  geom_vline(aes(xintercept = psig, color = "Bayes Mode"), lwd = 1, lty = 2) +
  geom_vline(aes(xintercept = hpdl, color = "HPD Interval"), lwd = 1, lty = 4) +
  geom_vline(aes(xintercept = hpdu, color = "HPD Interval"), lwd = 1, lty = 4) +
  theme_light() +
  scale_color_manual(values = colors) +
  labs(color = "")
wt <- ggplot(cats, aes(x = Weight, fill = Treatment)) +
  geom_histogram(alpha = 0.5) +
  theme_light() +
  theme(legend.position = "none") +
  labs(x = "Weight (grams)")
ag <- ggplot(cats, aes(x = Age, fill = Treatment)) +
  geom_histogram(alpha = 0.5) +
  theme_light() +
  labs(x = "Age (days)")
grid.arrange(wt, ag, ncol = 2)
cancers <- cats %>%
  group_by(MyeloidLeukemia, HarderianTumor, PreTLymphoma) %>%
  count()
xtable(cancers, caption = "Counts of cancer status") %>%
  xtable2kable(booktabs = T, include.rownames = FALSE, table.placement = NULL) %>%
  kable_styling(full_width = F, latex_options = "hold_position")
box1 <- ggplot(cats, aes(x = as.factor(Cataracts), y = Age, fill = Sex)) +
  geom_boxplot() + facet_wrap(vars(Sex)) +
  theme_light() + scale_fill_startrek() +
  labs(x = "Cataracts Status")
box2 <- ggplot(cats, aes(x = as.factor(Cataracts), y = Age, fill = as.factor(Cataracts))) +
  geom_boxplot() + facet_wrap(vars(Treatment)) +
  theme_light() + scale_fill_startrek() +
  labs(x = "Cataracts Status", fill = "Cataracts")
grid.arrange(box1, box2, ncol = 2)
library(equatiomatic)
# mixed model binomial logistic regression
full_mod <- glmer(Cataracts ~ Treatment + Sex + Weight + CoatColor + BCS + HarderianTumor +
                    MyeloidLeukemia + PreTLymphoma + (1|Family), data = cats, family = binomial,
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
no_interaction <- glmer(Cataracts ~ Treatment + Sex + (1|Family), data = cats, family = binomial)
final_mod <- mod
simple_mod <- glmer(Cataracts ~ Treatment + (1|Family), data = cats, family = binomial)
fixed_mod <- glm(Cataracts ~ Treatment*Sex, data = cats, family = binomial)
extract_eq(full_mod)
extract_eq(no_interaction)
extract_eq(simple_mod)
extract_eq(fixed_mod)
aics <- c(round(AIC(final_mod),2), round(AIC(full_mod),2), round(AIC(fixed_mod),2), round(AIC(no_interaction),2), round(AIC(simple_mod),2))
models <- c("Final Model", "Full Mixed Model",
            "Fixed Effects Model", "Mixed Model with no Interaction", "Base Model")
dt <- cbind(models,aics)
colnames(dt) <- c("Model", "AIC")
xtable(dt, caption = "Model selection via AIC comparision") %>%
  xtable2kable(booktabs = T, include.rownames = FALSE, table.placement = NULL) %>%
  kable_styling(full_width = F, latex_options = "hold_position")
## assumption of normally distributed random effect
reff <- as.data.frame(ranef(mod)$Family) %>% rename(re = `(Intercept)`)
ggplot(reff, aes(sample = re)) +
  stat_qq() + stat_qq_line() +
  theme_light() +
  labs(x = "qnorm", y = "random intercept")
## assumption of no over-dispersion
rp <- residuals(final_mod, type = "pearson")
rat <- sum(rp^2)/df.residual(final_mod)
baysum <- summary(mcmc.model)
baystats <- baysum$statistics
xtable(baystats, caption = "Bayesian Model output on log odds scale") %>%
  xtable2kable(booktabs = T, include.rownames = FALSE, table.placement = NULL) %>%
  kable_styling(full_width = F, latex_options = "hold_position")
params = c("b0", "b1", "b2", "b3", "b4", "b5", "sigma^2")
traplot(mcmc.model, parms = params)
acf(mcmc.model[[3]][,1], main = "b0")
acf(mcmc.model[[3]][,2], main = "b1")
acf(mcmc.model[[3]][,3], main = "b2")
acf(mcmc.model[[3]][,4], main = "b3")
acf(mcmc.model[[3]][,5], main = "b4")
acf(mcmc.model[[3]][,6], main = "b5")
acf(mcmc.model[[3]][,8], main = "sigma2")
gelman.diag(mcmc.model)
#gelman.plot(mcmc.model)
colors <- c("Density" = "#5C88DAFF", "Mode" = "#84BD00FF", "HPD" = "#CC0C00FF")
p0 <- ggplot(posts, aes(x = b0)) +
  geom_density(aes(color = "Density", fill = "Density"), alpha = 0.5) +
  geom_vline(aes(xintercept = est_tab[1,6], color = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = est_tab[1,7], color = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = est_tab[1,4], color = "Mode")) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_light() +
  theme(plot.title = element_text(size= 12)) +
  labs(x = "Log odds",
       title = "b_0")
p1 <- ggplot(posts, aes(x = b1)) +
  geom_density(aes(color = "Density", fill = "Density"), alpha = 0.5) +
  geom_vline(aes(xintercept = est_tab[2,6], color = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = est_tab[2,7], color = "HPD"), linetype = "dashed") +
  
  geom_vline(aes(xintercept = est_tab[2,4], color = "Mode")) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_light() +
  theme(plot.title = element_text(size= 12)) +
  labs(x = "Log odds",
       title = "b_1")
p2 <- ggplot(posts, aes(x = b2)) +
  geom_density(aes(color = "Density", fill = "Density"), alpha = 0.5) +
  geom_vline(aes(xintercept = est_tab[3,6], color = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = est_tab[3,7], color = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = est_tab[3,4], color = "Mode")) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_light() +
  theme(plot.title = element_text(size= 12)) +
  labs(x = "Log odds",
       title = "b_2")
p3 <- ggplot(posts, aes(x = b3)) +
  geom_density(aes(color = "Density", fill = "Density"), alpha = 0.5) +
  geom_vline(aes(xintercept = est_tab[4,6], color = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = est_tab[4,7], color = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = est_tab[4,4], color = "Mode")) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_light() +
  theme(plot.title = element_text(size= 12)) +
  labs(x = "Log odds",
       title = "b_3")
p4 <- ggplot(posts, aes(x = b4)) +
  geom_density(aes(color = "Density", fill = "Density"), alpha = 0.5) +
  geom_vline(aes(xintercept = est_tab[5,6], color = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = est_tab[5,7], color = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = est_tab[5,4], color = "Mode")) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_light() +
  theme(plot.title = element_text(size= 12)) +
  labs(x = "Log odds",
       title = "b_4")
p5 <- ggplot(posts, aes(x = b5)) +
  geom_density(aes(color = "Density", fill = "Density"), alpha = 0.5) +
  geom_vline(aes(xintercept = est_tab[6,6], color = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = est_tab[6,7], color = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = est_tab[6,4], color = "Mode")) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_light() +
  theme(plot.title = element_text(size= 12)) +
  labs(x = "Log odds",
       title = "b_5")
p6 <- ggplot(posts, aes(x = sigma2)) +
  geom_density(aes(color = "Density", fill = "Density"), alpha = 0.5) +
  geom_vline(aes(xintercept = est_tab[7,6], color = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = est_tab[7,7], color = "HPD"), linetype = "dashed") +
  geom_vline(aes(xintercept = est_tab[7,4], color = "Mode")) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_light() +
  theme(plot.title = element_text(size= 12)) +
  labs(x = "Log odds",
       title = "sigma^2")
grid.arrange(p0, p1, p2, p3, p4, p5, p6, ncol = 2)
cats_emms <- emmeans(mod, ~ Treatment | Sex, infer = TRUE, type = "response")
emmip(cats_emms, Treatment ~ Sex) +
  theme_light() + scale_color_startrek()
cats_emms <- emmeans(final_mod, ~ Treatment | Sex, infer = TRUE, type = "response")
podds <- pairs(cats_emms, reverse = TRUE)
podds <- as.data.frame(confint(podds)) %>%
  rename(Contrast = contrast, Lower = asymp.LCL, Upper = asymp.UCL)
or1_plot <- ggplot(podds, aes(x = odds.ratio, y = Contrast, xmin = Lower, xmax = Upper)) +
  geom_errorbarh(aes(height = 0.2, color = Sex),
                 position = position_dodge(0.3), lwd = 1) +
  geom_point(aes(color = Sex), position = position_dodge(0.3)) +
  geom_text(aes(label = round(odds.ratio, 2), color = Sex), position = position_dodge(0.7)) +
  geom_vline(aes(xintercept = 1), color = "#84BD00FF", lty = 2) +
  theme_light() +
  theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 0.5)) +
  scale_color_startrek() +
  labs(x = "Odds Ratio")
or1_plot
cats_emms1 <- emmeans(final_mod, ~Sex | Treatment, infer = TRUE, type = "response")
podds1 <- pairs(cats_emms1, reverse = TRUE)
podds1
podds1 <- as.data.frame(confint(podds1)) %>%
  rename(Contrast = contrast, Lower = asymp.LCL, Upper = asymp.UCL)
or2_plot <- ggplot(podds1, aes(x = odds.ratio, y = Contrast, xmin = Lower, xmax = Upper)) +
  geom_errorbarh(aes(height = 0.2, color = Treatment),
                 position = position_dodge(0.3), lwd = 1) +
  geom_point(aes(color = Treatment), position = position_dodge(0.3)) +
  geom_text(aes(label = round(odds.ratio, 2), color = Treatment),
            position = position_dodge(0.3), vjust = -1) +
  theme_light() +
  geom_vline(aes(xintercept = 1), color = "#FFCD00FF", lty = 2) +
  scale_color_startrek() +
  labs(x = "Odds Ratio")
or2_plot
options(width=100)







