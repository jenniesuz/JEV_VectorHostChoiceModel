
rm(list = ls())
# packages required
library(reshape)
library(patchwork)
library(parallel)
library(here)
library(lhs)
library(parallel)
library(tidyverse)
library(ggplot2)
#library(sensitivity)

#source(here("scripts/parameters.R"))
source(here("scripts/supportingFunctions.R"))
source(here("scripts/R0_function.R"))

set.seed(202510)

mc.cores <- detectCores() / 2

# RColorBrewer::brewer.pal(8, "Set1")
cols <- c("black",  "#377EB8", "#E41A1C", "#4DAF4A", "#377EB8", "#FF7F00")
shapes <- c(16, 15, 17, 18)

agg_labels <- c(
  "0" = "No aggregation",
  "0.9" = expression(m[p] == 0.9),
  "0.1" = expression(m[p] == 0.1)
)

agg_labels2 <- c(
  "0" = "No aggregation",
  "0.9" = expression(atop("Weak aggregation", m[p] == 0.9)),
  "0.1" = expression(atop("Strong aggregation", m[p] == 0.1))
 
)

agg_labels3 <- c(
  "0" = "No aggregation",
  "0.9" = expression(paste("Weak aggregation, ", m[p] == 0.9)),
  "0.1" = expression(paste("Strong aggregation, ", m[p] == 0.1))
)



#*################################################################
#******************** R0 Baseline scenario **************************
#*################################################################
#* Mosquitoes have equal access to all patches, and all patches 
#* are equally (equal no. of hosts in each patch) suitable 
#* (and no 'long-range' host-seeking behavior)
# a. Within-patch host choice is independent of vector density
# plot_heatmap <- expand_grid(prefComp = seq(0, 0.5, length.out = 300), 
#             de2Comp = seq(1, 50, length.out = 300),
#             n_patches = 20 ) |>  
#   rowwise() |> 
#   mutate(H_de = de2Comp*startingSh,
#          r0 = extendedR0(H_c_vect = rep(startingSh/n_patches, n_patches),
#                          H_de_vect = rep(H_de/n_patches, n_patches),
#                          prefComp = prefComp,
#                          infC = 0
#          ) ) |>  
#   ggplot(aes(x = prefComp, y = de2Comp, z = r0)) +
#   geom_raster(aes(fill = r0) ) +
#   geom_contour(color = "gray10", breaks = 1, linewidth = 1) +
#   scale_fill_distiller(palette = "Spectral", direction = -1) +
#   labs(x = expression("Preference for competent hosts, " * eta[x]),
#        y = "Dead-end to competent host ratio",
#        fill = expression(R[0])) +
#   annotate("text", x = 0.2, y = 20, label = "R[0] < 1", parse = TRUE, size = 8, color = "black") +
#   annotate("text", x = 0.4, y = 2.5, label = "R[0] > 1", parse = TRUE, size = 8, color = "black") +
#   theme_minimal() +
#   theme(axis.text=element_text(size=12),
#       legend.text=element_text(size=12),
#       axis.title=element_text(size=12))
# 
# plot_heatmap
# ggsave("./outputs/simpleR0heatmap.pdf", width = 12, height = 10, units = "in")


startingSh <- 500       # Number of susceptible hosts
startingSv <- 10000     # Number of susceptible vectors
mosquitoes_per_comp_host <- startingSv/startingSh


expand_grid(prefComp = 0.1, 
            de2Comp = 1,
            n_patches = 1:50 ) |>  
  rowwise() |> 
  mutate(H_de = de2Comp*startingSh,
         r0 = extendedR0(H_c_vect = rep(startingSh/n_patches, n_patches),
                         H_de_vect = rep(H_de/n_patches, n_patches),
                         prefComp = prefComp,
                         infC = 0
         ) ) |>  
  ggplot(aes(x = n_patches, y = round(r0, 3))) +
  geom_line() 

# Mosquitoes have equal access to all patches, but patches vary in suitability,
# with higher host densities attracting more mosquitoes via long-range cues.
# Scenario 2: Patch comparison

## TJ: loop over both preference values (0.1 = main text, 0.9 = supplementary)
pref_values <- c(0.1, 0.9)


scenario_2 <- expand_grid(
  prefComp = pref_values    ## Change hardcoded 0.9
  ,de2Comp = seq(1, 100, length.out = 100)
  ,hostDist = c("exp:hostx", "exp:hosty")
  ,nPatches = seq(1, 50, 1)
  ,infC = c(0, 0.1, 0.9)
  ,decay = c(0.2, 1)
  ) |> 
  mutate(
    scenario = case_when(
    # hostDist == "equal" ~ "Equal",
      hostDist == "exp:hostx"  ~ paste("Competent hosts aggregated, \n (decay = ", decay,")",sep = ""),
      hostDist == "exp:hosty"  ~ paste("Dead-end hosts aggregated, \n (decay = ", decay,")",sep = "") 
      )  ) |> 
  distinct() 
#  refactor_host_dist()
scenario_2


tictoc::tic()
scenario_2_list <- pbmcapply::pbmclapply(1:nrow(scenario_2), function(i) {
  
  # Generate scenario R0 replicates
  idx <- scenario_2[i, ]
  
  scenario_r0_list  <- replicate(200, {
    icHost <- initialCondsFunc(
      nPatches = idx$nPatches,
      totalSh  = startingSh,
      de2comp  = idx$de2Comp,
      hostDist = idx$hostDist,
      decay    = idx$decay
    )
    
    extendedR0(
      H_c_vect = unname(icHost[[1]][startsWith(names(icHost[[1]]), 'Sh')]), 
      H_de_vect = as.vector(icHost[[2]]),
      prefComp    = idx$prefComp,
      infC        = idx$infC
      )
  })
  
  # Compute baseline R0 (deterministic)
  baseline_icHost <- initialCondsFunc(
    nPatches = idx$nPatches,
    totalSh = startingSh,
    de2comp  = idx$de2Comp,
    hostDist = "equal",
    decay = NA
  )
  
  baseline_r0 <- extendedR0(
    H_c_vect = unname(baseline_icHost[[1]][startsWith(names(baseline_icHost[[1]]), 'Sh')]),
    H_de_vect = as.vector(baseline_icHost[[2]]),
    prefComp    = idx$prefComp,
    infC = 0
  )
  
  idx |> 
    mutate(
    baseline_r0 = baseline_r0,
    scenario_r0 = list(scenario_r0_list)  # Store replicates as list-column
  )

}, mc.cores = mc.cores, mc.set.seed = TRUE )
tictoc::toc()


scenario_res_df <- bind_rows(scenario_2_list) |> 
  unnest_longer(scenario_r0) |> 
  mutate(rel_change = scenario_r0 / baseline_r0,
         abs_change = scenario_r0 - baseline_r0 ) |> 
  group_by(prefComp, de2Comp, nPatches, hostDist, infC, scenario) |> 
  dplyr::summarise(mean_r0 = mean(scenario_r0),
                   
                   mean_rel = mean(rel_change),
                   lower_rel = Rmisc::CI(rel_change)["lower"],
                   upper_rel = Rmisc::CI(rel_change)["upper"],
                   
                   mean_abs = mean(abs_change),
                   lower_abs = Rmisc::CI(abs_change)["lower"],
                   upper_abs = Rmisc::CI(abs_change)["upper"],
                   .groups = "drop" ) 



##  function, generates plots for both preference values (main + supplementary)
plot_patches_list <- lapply(pref_values, function(pref) {
  scenario_res_df |>  
    filter(de2Comp == 1, prefComp == pref) |>  
    ggplot(aes(x = nPatches, y = mean_abs, ymin = lower_abs, ymax = upper_abs, 
               fill = factor(infC,levels=c("0","0.9","0.1")), 
               colour = factor(infC,levels=c("0","0.9","0.1")) ) ) +
    geom_ribbon(alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~scenario, scales = "fixed") +
    scale_colour_manual(values = cols, name = "", labels = agg_labels3) +
    scale_fill_manual(values = cols, name = "", labels = agg_labels3) +
    scale_linetype_manual(values = c("solid", "dashed") ) +
    labs(x = "Number of patches", y = expression(Delta~R[0])) +
    theme_minimal() +
    plotThemeFunc(leg.pos = "bottom")
})
names(plot_patches_list) <- pref_values

# Main figure (pref = 0.1)
plot_patches <- plot_patches_list[["0.1"]]
plot_patches
# ggsave("./outputs/R0PatchesbyDistn_main.pdf", width = 12, height = 10, units = "in", dpi = 500)

# Supplementary figure (pref = 0.9)
plot_patches_supp <- plot_patches_list[["0.9"]]
plot_patches_supp
# ggsave("./outputs/R0PatchesbyDistn_supp.pdf", width = 12, height = 10, units = "in", dpi = 500)

# plot_host_comp <- scenario_res_df |>  
#   filter(nPatches == 20) |>  
#   ggplot(aes(x = de2Comp, y = mean_abs, ymin = lower_abs, 
#              ymax = upper_abs, fill = factor(infC,levels=c("0","0.9","0.1")), 
#              colour = factor(infC,levels=c("0","0.9","0.1")) )) +
#   geom_ribbon(alpha = 0.2, colour = NA) +
#   geom_line(linewidth = 0.8) +
#   facet_wrap(~scenario, scales = "fixed") +
#   scale_colour_manual(values = cols, name = "", labels = agg_labels3) +
#   scale_fill_manual(values = cols, name = "", labels = agg_labels3) +
#   scale_linetype_manual(values = c("solid", "dashed") ) +
#   labs(x = "Dead-end to competent host ratio", y = expression(Delta~R[0])  ) +
#   theme_minimal() +
#   scale_x_log10() +
#   plotThemeFunc(leg.pos = "bottom")
# 
# plot_host_comp
# # ggsave("./outputs/R0HostCompByDistnPref0.9.pdf", width = 12, height = 10, units = "in", dpi = 500)


## TJ: function, now generates host composition plots for both preference values
plot_host_comp_list <- lapply(pref_values, function(pref) {
  scenario_res_df |>  
    filter(nPatches == 20, prefComp == pref) |>  
    ggplot(aes(x = de2Comp, y = mean_abs, ymin = lower_abs, ymax = upper_abs, 
               fill = factor(infC,levels=c("0","0.9","0.1")), 
               colour = factor(infC,levels=c("0","0.9","0.1")) )) +
    geom_ribbon(alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~scenario, scales = "fixed") +
    scale_colour_manual(values = cols, name = "", labels = agg_labels3) +
    scale_fill_manual(values = cols, name = "", labels = agg_labels3) +
    scale_linetype_manual(values = c("solid", "dashed") ) +
    labs(x = "Dead-end to competent host ratio", y = expression(Delta~R[0])) +
    theme_minimal() +
    scale_x_log10() +
    plotThemeFunc(leg.pos = "bottom")
})
names(plot_host_comp_list) <- pref_values

# Main figure (pref = 0.1)
plot_host_comp <- plot_host_comp_list[["0.1"]]
plot_host_comp
# ggsave("./outputs/R0HostCompByDistn_main.pdf", width = 12, height = 10, units = "in", dpi = 500)

# Supplementary figure (pref = 0.9)
plot_host_comp_supp <- plot_host_comp_list[["0.9"]]
plot_host_comp_supp
# ggsave("./outputs/R0HostCompByDistn_supp.pdf", width = 12, height = 10, units = "in", dpi = 500)


###############################################################################
#*########################## Sensitivity analysis with R0 #####################
###############################################################################
# Define parameter ranges
lower <- c(
  mu_m = 1/30,            # Mosquito mortality rate (per day)
  mu_h = 1/(5*365),     # Host mortality rate (per day), was 1/365 (1yr), now 1/(5*365) for 5-year lifespan
  alpha = 1/5,           # bite rate (per day)
  p_hv = 0.001,               # Probability of infection from host to mosquito
  p_vh = 0.001,               # Probability of infection from mosquito to host
  sigma = 1/14,           # Host recovery rate (per day)
  prefComp = 0.001,       # Mosquito preference for competent hosts
  mos_per_comp_host = 1,# Mosquito density per competent host
  infC = 0.001,               # Mosquito aggregation parameter
  de2Comp = 1,            # Dead-end to competent host ratio
  n_patches = 1,
  decay = 0.001
)

upper <- c(
  mu_m = 1/10,            # Mosquito mortality rate (per day)
  mu_h = 4/365,       # Host mortality rate (per day), was 10/365 (typo), now 4/365 for 3-months lifespan
  alpha = 1,              # bite rate (per day)
  p_hv = 1,               # Probability of infection from host to mosquito
  p_vh = 1,               # Probability of infection from mosquito to host
  sigma = 1/4,              # Host recovery rate (per day)
  prefComp = 1,         # Mosquito preference for competent hosts
  mos_per_comp_host = 100,# Mosquito density per competent host
  infC = 1,               # Mosquito aggregation parameter
  de2Comp = 100,          # Dead-end to competent host ratio
  n_patches = 50,
  decay = 2
)

# Setup
n_replicates <- 100      # number of equiprobable intervals
n_samples <- 1000         # Samples per replicate
host_dist_params <- c("equal", "exp:hostx", "exp:hosty")

host_dist_labels <- c(
  "equal" = "Equal",
  "exp:hostx" = "Exp: competent host x",
  "exp:hosty" = "Exp: dead-end host y"
)

parameter_labels <- c(
  "mu_m" = "bold('Mosquito')~bold('mortality')~(bold(mu[m]))", 
  "mu_h" = "bold('Host')~bold('mortality')~(bold(mu[h]))",
  "alpha" = "bold('Biting')~bold('rate')~(bold(alpha))",
  "p_hv" = "bold('Trans')~bold('prob:')~bold('host\U2013mos')~(bold(p[hv]))",
  "p_vh" = "bold('Trans')~bold('prob:')~bold('mos\U2013host')~(bold(p[vh]))",
  "sigma" = "bold('Host')~bold('recovery')~(bold(sigma))",
  "prefComp" = "bold('Host')~bold('preference.')~(bold(eta[x]))",
  "mos_per_comp_host" = "bold('Mosquitoes')~bold('per')~bold('host')~(bold(N[m]/H[x]))",
  "infC" = "bold('Interference')~(bold(m[p]))",
  "n_patches" = "bold('Number')~bold('of')~bold('patches')~(bold(n))",
  "de2Comp" = "bold('Host')~bold('ratio')~(bold(H[y]/H[x]))",
  "decay" = "bold('Exponential')~bold('decay')~(bold(lambda))"
)


# Main analysis function
run_replicate <- function(replicate_id, host_dist) {
  # Generate LHS samples for this replicate
  lhs_design <- randomLHS(n = n_samples, k = length(lower))
  
  # Transform to parameter space
  params.set <- t(t(lhs_design) * (upper - lower) + lower)
  colnames(params.set) <- names(lower)
                  
                  
  # Calculate R0 for all samples
  R0_values <- lapply(1:n_samples, function(i) {
    icHost <- initialCondsFunc(
      nPatches = params.set[i, "n_patches"],
      totalSh = startingSh,
      de2comp = params.set[i, "de2Comp"],
      hostDist = host_dist,
      decay = ifelse(host_dist == "equal", NA_real_, params.set[i, "decay"])
    )
    
    extendedR0(
      prefComp = params.set[i, "prefComp"],
      sigma = params.set[i, "sigma"],
      p_vh = params.set[i, "p_vh"],
      p_hv = params.set[i, "p_hv"],
      alpha = params.set[i, "alpha"],
      mu_h = params.set[i, "mu_h"],
      mu_m = params.set[i, "mu_m"],
      mos_per_comp_host = params.set[i, "mos_per_comp_host"],
      H_c_vect = unname(icHost[[1]][startsWith(names(icHost[[1]]), 'Sh')]),
      H_de_vect = as.vector(icHost[[2]]),
      infC = params.set[i, "infC"]
    ) |> replace_na(0)
  })
  
  # Return as data frame with replicate ID
  data.frame(params.set, R0 = unlist(R0_values), replicate = replicate_id)
}

# Run full analysis for all host distributions
full_results <- lapply(host_dist_params, function(host_dist) {
  replicate_results <- mclapply(1:n_replicates, function(r) {
    run_replicate(r, host_dist)
  }, mc.cores = mc.cores)
  
  # Combine and add host distribution label
  bind_rows(replicate_results) |> 
    mutate(hostDist = host_dist,
           scenario = factor(host_dist, levels = names(host_dist_labels), 
                             labels = host_dist_labels))
}) |> bind_rows()

# Uncertainty analysis (modified for replicates)
uncertainty_data <- full_results |>
  mutate(R0_category = ifelse(R0 < 1, "<1", ">1")) |>
  pivot_longer(cols = -c(hostDist, R0, R0_category, replicate, scenario), 
               names_to = "Parameter", values_to = "Value") |>
  mutate(Parameter_Math = factor(Parameter, levels = names(parameter_labels),
                                 labels = parameter_labels))

uncertainty_data |>  
  ggplot(aes(x = R0_category, y = Value, colour = scenario)) +
  geom_boxplot(width = .5, alpha = 1, outlier.alpha = 0.1, outlier.size = 0.1, notch = T) +
  facet_wrap(~Parameter_Math, scales = "free_y", labeller = label_parsed, nrow = 4) +
  scale_colour_manual(values = cols) +
  labs(x = R[0]~category,
       y = "Parameter value",
       colour = "Host distribution") +
  theme_minimal() +
  plotThemeFunc("top")
# ggsave("./outputs/uncertaintyPlot.pdf", width = 12, height = 12, units = "in", dpi = 500)

# Sensitivity analysis with PRCC variability
prcc_results <- bind_rows(
  lapply(host_dist_params, function(dist) {
    dist_data <- full_results |> filter(hostDist == dist)
    
    bind_rows(
      mclapply(unique(dist_data$replicate), function(r) {
        replicate_data <- dist_data |> filter(replicate == r)
        
        pcc_result <- pcc(replicate_data |> select(names(lower)), 
                          replicate_data$R0,
                          nboot = 1000, rank = TRUE)
        
        pcc_result$PRCC |> 
          as.data.frame() |> 
          rownames_to_column("Parameter") |> 
          mutate(replicate = r, hostDist = dist)
      }, mc.cores = mc.cores )
    )
  } )
)

# summarise results
prcc_summary <- prcc_results |> 
  group_by(hostDist, Parameter) |> 
  summarise(
    median_PCC = median(original),
    q05 = quantile(original, 0.05),
    q95 = quantile(original, 0.95),
    .groups = "drop"
  )

prcc_summary |> mutate(
  Parameter_Math = factor(Parameter, levels = names(parameter_labels), labels = parameter_labels),
  scenario = factor(hostDist, levels = names(host_dist_labels), labels = host_dist_labels) 
  ) |> 
  ggplot(aes(x = Parameter_Math, y = median_PCC, ymin = q05, ymax = q95, color = scenario)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Partial rank correlation coefficient") +
  scale_colour_manual(values = cols, name = "Host distribution") +
  theme_minimal() +
  plotThemeFunc(c(0.22, 0.3)) +
  theme(
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.key.size = unit(1, "cm")
    ) +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  coord_flip()

# ggsave("./outputs/prcc.pdf", width = 12, height = 12, units = "in", dpi = 500, device = cairo_pdf)
