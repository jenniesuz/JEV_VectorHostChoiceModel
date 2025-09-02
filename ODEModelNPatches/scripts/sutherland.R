
library(tidyverse)
library(patchwork) 

#* Pipeline
#* 1. Prepare data
#* 2. Write model function and write different potential likelihoods
#* 3. Run and select model based on AIC values
#* 4. Predict data using the best model (lowest AIC)
#* 5. Make visualisations, including one with m_p = 1 and no aggregation (evenly distributed)

#************ prepare data **************
#*
tritaeHostsHH <- read.csv("tritaeHostsHH.csv") |> 
  filter(Females > 0) |>  # Remove households with zero mosquitoes (no need to add 1 to logs then?)
  mutate(
    TotalHostsCE = as.numeric(TotalHostsCE) # Ensure numeric type
  )

# Calculate totals
total_hosts <- sum(tritaeHostsHH$TotalHostsCE, na.rm = TRUE)
total_mosquitoes <- sum(tritaeHostsHH$Females, na.rm = TRUE)

tritaeHostsHH <- tritaeHostsHH |> 
  mutate(
    a_x = TotalHostsCE / total_hosts, # normalise host numbers
    delta_obs = Females / total_mosquitoes 
  )

#* Sutherland aggregation model
sutherland_model <- function(m_p, a_x) {
  numerator <- a_x^(1/m_p)
  numerator / sum(numerator)
}

#************ Define likelihoods **************
# Define different likelihoods

likelihoods <- list(
  
  # Log-normal likelihood
  lognormal = function(params) {
    m_p <- params[1]
    sigma <- params[2]
    pred_counts <- total_mosquitoes * sutherland_model(m_p, tritaeHostsHH$a_x)
    -sum(dnorm(log(tritaeHostsHH$Females), mean = log(pred_counts),
               sd = sigma, log = TRUE))
  },
  
  # Log10 transformed likelihood
  log10normal = function(params) {
    m_p <- params[1]
    sigma <- params[2]
    pred_counts <- total_mosquitoes * sutherland_model(m_p, tritaeHostsHH$a_x)
    -sum(dnorm(log10(tritaeHostsHH$Females), mean = log10(pred_counts),
               sd = sigma, log = TRUE))
  },
  
  # Negative binomial (overdispersed counts?)
  negbin = function(params) {
    m_p <- params[1]
    size <- params[2] # dispersion parameter
    pred_counts <- total_mosquitoes * sutherland_model(m_p, tritaeHostsHH$a_x)
    -sum(dnbinom(tritaeHostsHH$Females, mu = pred_counts, size = size, log = TRUE))
  }
)

#************ Fit and compare models **************

fits <- list(
  log10normal = optim(c(0.1, 1), fn = likelihoods$log10normal, method = "L-BFGS-B",
    lower = c(0.01, 0.1), upper = c(1, 10), hessian = TRUE),
  
  negbin = optim(c(0.1, 10), fn = likelihoods$negbin, method = "L-BFGS-B",
    lower = c(0.01, 0.1), upper = c(1, 100), hessian = TRUE)
)

# Calculate AIC
aic_values <- sapply(names(fits), function(model) {
  k <- length(fits[[model]]$par)
  2 * fits[[model]]$value + 2 * k
})

aic_values

# Model comparison table
data.frame(
  Model = names(fits),
  Parameters = sapply(fits, function(x) paste(round(x$par[1], 3), collapse = ", ")),
  AIC = round(aic_values, 1),
  row.names = NULL
) |> arrange(AIC) |> 
  mutate(DeltaAIC = AIC - min(AIC))

#************ Select best model and predict **************

best_model <- names(which.min(aic_values))
best_model

# Add predictions to data
tritaeHostsHH_long <- tritaeHostsHH |> 
  mutate(
    pred_log10normal = total_mosquitoes * sutherland_model(fits$log10normal$par[1], a_x),
    pred_negbin = total_mosquitoes * sutherland_model(fits$negbin$par[1], a_x),
    # Add extremes
    pred_m1 = total_mosquitoes * sutherland_model(1, a_x), # m_p=1
    pred_m0 = total_mosquitoes * rep(1/n(), n()) # m_p=0 (even distribution)
  ) |> 
  pivot_longer(cols = pred_log10normal:pred_m0, # c(30, 29, 28),
               names_to = "model", 
               values_to = "pred" )
  

model_labels <- c(
  pred_log10normal   = paste0("Log10-normal*', '*italic(m)[p]==", round(fits$log10normal$par[1], 2)),
  
  pred_negbin      = paste0("Negative*' '*Binomial*', '*italic(m)[p]==", round(fits$negbin$par[1], 2)),
  pred_m1          = "italic(m)[p]==1",
  pred_m0          = "Evenly~distributed"
)

# Convert model names to factor with expression labels
tritaeHostsHH_long <- tritaeHostsHH_long |> 
  mutate(
    model_label = factor(
      model,
      levels = names(model_labels),
      labels = model_labels
    )
  )


#************ Plots **************
#*
custom_colors <- c(RColorBrewer::brewer.pal(3, "Set1"), "grey50", "black")
custom_colors

# Residual plot
resid_plot <- tritaeHostsHH_long %>% 
  filter(!(model %in% c("pred_m1", "pred_m0") ) ) %>% 
  ggplot(aes(x = log1p(pred), y = log10(Females) - log10(pred), colour = model_label)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "log10(Predicted)", y = "log10(Residual)", title = "Model residuals") +
  scale_color_manual(values = custom_colors, name = "", labels = function(x) parse(text = x) ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        legend.text=element_text(size=10),
        strip.text=element_text(size=10),
        axis.title=element_text(size=12)
        )

ref_plot <- tritaeHostsHH_long %>% 
  filter(!(model %in% c("pred_m1", "pred_m0") ) ) %>% 
  ggplot(aes(x = log10(Females), y = log10(pred), colour = model_label)) +
  geom_point(size = 2, alpha = 0.8, show.legend = F) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "log10(Observed females)", y = "log10(Predicted females)",
       title = "Observed vs. predicted (fitted models)" ) +
  scale_color_manual(values = custom_colors, name = "", labels = function(x) parse(text = x) ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        legend.text=element_text(size=10),
        strip.text=element_text(size=10),
        axis.title=element_text(size=12)
        )

# Model comparison plot
main_plot <- ggplot(tritaeHostsHH_long, aes(x = TotalHostsCE, y = Females)) +
  geom_point(alpha = 1, size = 1.5) +
  geom_line(aes(y = pred + 1, colour = model_label), linewidth = 1) +
  scale_y_log10() +
  scale_color_manual(values = custom_colors, name = "", labels = function(x) parse(text = x) ) +
  labs(x = "Total hosts in cattle equivalents", y = "Female mosquitoes") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        legend.text=element_text(size=10),
        strip.text=element_text(size=10),
        axis.title=element_text(size=12)
        )

# Combine plots
main_plot/(ref_plot | resid_plot)  +
  plot_layout(heights = c(1.2, 1))

# ggsave("./outputs/mosAggregationData.pdf", width = 12, height = 12, units = "in")

main_plot + 
  theme(legend.position = c(0.8, 0.2),
        axis.text=element_text(size=14),
        legend.text=element_text(size=14),
        axis.title=element_text(size=12)
  )
# ggsave("./outputs/fittedLinesOnMosAggData.pdf", width = 12, height = 10, units = "in")

