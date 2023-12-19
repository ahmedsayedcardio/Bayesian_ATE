#Load necessary packages (please make sure to have these packages installed beforehand)
library(dplyr)
library(tidyr)
library(readxl)
library(bayesplot)
library(bayestestR)
library(MetaStan)
library(brms)
library(tidybayes)
library(stringr)
library(ggplot2)
library(ggridges)
library(ggpubr)

#Import data from excel file
x <- read_xlsx("data_cancer_anticoags.xlsx")
#Store as dataframe
x <- as.data.frame(x)

#Provide data in "long" format suitable for brms
y <- create_MetaStan_dat(
  dat = x,
  armVars = c(responders = "e", sampleSize = "n"),
  nArmsVar = "n_arms"
)

#Store "long" format data in a separate object
data <- y$data_long

#Assign treatment variables
data$ttt <- c(-0.5, 0.5)
#Assign study variables
data$study <- rep(x$study, each = 2)

#Set formula to be used for logistic regression
meta_formula <- bf(
  responders | trials(sampleSize) ~ 0 + factor(study) + ttt + (ttt - 1 | study)
)

#Set priors
priors <- c(prior(normal(0, 10), class = b), #Baseline log-odds of ATE
            prior(normal(0,1), class = b, coef = "ttt"), #Treatment effect
            prior(cauchy(0,0.5), class = sd)) #Heterogeneity

# Specify the control parameters (to prevent divergences)
control <- list(adapt_delta = 0.95)

#Run model
b_model <- brm(data = data,
               family = binomial,
    formula = meta_formula,
    seed = 100, #Set seed for reproducibility
    prior = priors,
    control = control)

#Get probability of any benefit (OR < 1)
p_direction(b_model)


#Get desired OR to achieve a 1% ARD assuming a 3.99% risk (TKI study)
arr = 0.01
baseline_risk = baseline_risk_tki = 0.0399
desired_or <- ((baseline_risk - arr) / (1 - (baseline_risk - arr))) / (baseline_risk / (1 - baseline_risk))

#Get the probability of the treatment effect being equal to or greater than the desired OR
sigp_tki <- p_significance(b_model, log(desired_or))


#Do the same thing but assuming a 1.3% risk (ICI study)
baseline_risk = baseline_risk_ic = 0.013
desired_or <- ((baseline_risk - arr) / (1 - (baseline_risk - arr))) / (baseline_risk / (1 - baseline_risk))
sigp_ic <- p_significance(b_model, log(desired_or))




#Sample the posterior distribution of study-level estimates and the overall estimate of treatment effect
study_es <- b_model %>%
  spread_draws(r_study[study, ], b_ttt) %>%
  mutate(b_ttt = r_study + b_ttt, #Create treatment effect estimates for each study
         type = "Study-level estimate", #Clarify that this is the treatment effect for each study
         study = study %>% str_replace_all("\\.", " ")) #Remove unnecessary dots added

pooled_es <- spread_draws(b_model, b_ttt) %>% 
  mutate(study = " Overall Effect Size", #Clarify that this is the pooled/overall treatment effect
         type = "Pooled estimate") #Same

#Exponentiate to get odds instead of log-odds ratio
fp_data <- bind_rows(study_es, pooled_es) %>%
  mutate(b_ttt = b_ttt %>% exp)



#Create title and subtitles for the plot
main_title <- "Reconstructed forest plot for effect of anticoagulants on the risk of arterial thrombotic events."
subtitle <- "We used a binomial-normal heirarchical model with the following priors: treatment effect (log odds ratio) ~ N(0, 1),\nbaseline risk (log odds) ~ N(0, 10), heterogeneity (\U03C4) ~ Half-Cauchy(0, 0.5)."


#Plot
ggplot(data = fp_data[fp_data$study != "Pooled Effect Size", ],
       aes(y = study,
           x = b_ttt
)) +
  #Add Density plots
  geom_density_ridges(col = NA,
                      scale = 0.5,
                      alpha = 0.5,
                      ) +
  stat_density_ridges(geom = "density_ridges_gradient",
                      quantile_lines = TRUE,
                      color = NA,
                      scale = 0.5,
                      quantiles = c(0.962),
                      calc_ecdf = TRUE,
                      data = fp_data[fp_data$study == " Overall Effect Size", ],
                      aes(fill = factor(after_stat(quantile)))
  ) +
  #Set colors
  scale_fill_manual(values = c("lightblue", "salmon")) +
  #Create title
  ggtitle(main_title,
          subtitle = subtitle) +
  geom_vline(xintercept = 1, color = "black", 
             lwd = 1.25, linetype = 2) +
  #X and Y axes aesthetics
  scale_y_discrete(name = "Study") +
  scale_x_continuous(name = "Treatment Effect (Odds Ratio)",
                     trans = "log",
                     breaks = c(0.25, 0.5, 1, 2, 4)) +
  #Set reasonable Y axis limits
  coord_cartesian(xlim = c(0.25, 2)) +
  #Set theme
  theme_pubclean() +
  theme(text = element_text(size = 23),
        plot.title=element_text(face = "bold",hjust = 0.0, size = 20),
        plot.subtitle = element_text(face = "bold", size = 15, hjust = 0.0, color = "grey45"),
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 25, face = "bold"),
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "black", size = 1.2),
        plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
        legend.background = element_rect(fill = "transparent"),
        legend.position = "none",
        legend.text = element_text(size = 12, face = "bold"),
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(0.75, "cm"))

