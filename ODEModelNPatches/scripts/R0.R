
rm(list = ls())
# packages required
library(reshape)
library(patchwork)
library(parallel)
library(here)
library(lhs)
library(parallel)
library(sensitivity)
library(tidyverse)

# read in scripts
source(here("scripts/parameters.R"))
source(here("scripts/supportingFunctions.R"))

set.seed(202510)

mc.cores <- detectCores() / 2

# RColorBrewer::brewer.pal(8, "Set1")
cols <- c("black", "#E41A1C",  "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
shapes <- c(16, 15, 17, 18)

agg_labels <- c(
  "0" = "No aggregation",
  "0.1" = expression(m[p] == 0.1),
  "0.9" = expression(m[p] == 0.9)
)

agg_labels2 <- c(
  "0" = "No aggregation",
  "0.1" = expression(atop("Strong aggregation", m[p] == 0.1)),
  "0.9" = expression(atop("Weak aggregation", m[p] == 0.9))
)

agg_labels3 <- c(
  "0" = "No aggregation",
  "0.1" = expression(paste("Strong aggregation, ", m[p] == 0.1)),
  "0.9" = expression(paste("Weak aggregation, ", m[p] == 0.9))
)

###############################################################################
####################### Plots host distribution in patches ####################
###############################################################################

# Define parameter ranges
parameter_grid <- expand_grid(
  n_patches = 1:50                    # Fixed number of patches
  ,hostDist = c("equal", "exp:hostx", "exp:hosty")   # Distribution types
  ,infC = c(0, 0.1, 0.9)
  ,decay = c(0.2, 1)
  ,de2comp = seq(1, 50, 1) ) |> 
  mutate(decay = ifelse(hostDist == "equal", NA, decay) ) |> 
  distinct()

all_distributions <- bind_rows(
  pbmcapply::pbmclapply(1:nrow(parameter_grid), function(i){ # for each row of parameter_grid
    cbind(parameter_grid[i, ],                              # take the row
          generate_distribution_data(                       # apply the function
            n_patches = parameter_grid$n_patches[i],
            totalSh = startingSh,
            de2comp = parameter_grid$de2comp[i],
            hostDist = parameter_grid$hostDist[i],
            infC = parameter_grid$infC[i],
            Nm = startingSv,
            decay = parameter_grid$decay[i]
          )
    )
  }, mc.cores = mc.cores )
) |>  
  mutate(
    scenario = case_when(
    hostDist == "equal" ~ "Equal",
    hostDist == "exp:hostx"  ~ paste("Exp: comp. (x), decay =", decay),
    hostDist == "exp:hosty"  ~ paste("Exp: dead-end (y), decay =", decay) 
    ),
    mosAgg = paste("Interference = ", infC)
  ) 

# Apply to data
all_distributions <- all_distributions |>  
  refactor_host_dist()


# ----- Plot A: Host Distribution per Patch -----
# This plot shows the number of hosts in each patch, distinguishing competent (H_c) and dead-end (H_de) hosts.
plot_a <- all_distributions |>  
  filter(n_patches == 10, infC == 0, de2comp == 1) |>  
  pivot_longer(cols = c(H_c, H_de), names_to = "species", 
               values_to = "host_count") |>  
  mutate(species = factor(species, levels = c("H_c", "H_de"), 
                          labels = c("Competent hosts", "Dead-end hosts"))) |>  
  ggplot(aes(x = factor(patch), y = host_count, fill = species)) +
  geom_bar(stat = "identity") +
  facet_wrap(~scenario1, nrow = 3, as.table = F, drop = F) +
 # scale_fill_brewer(palette = "Dark2", name = "Host species") +
  scale_fill_manual(values = cols[c(3, 5)], name = "Host species") +
  labs(x = "Patch", y = "Number of hosts", 
       title = "Distribution of host composition across patches") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.9),
        legend.background = element_rect(fill = NA, colour = NA),
        axis.line = element_line(color = 'black'),
        axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        strip.text=element_text(size=12),
        axis.title=element_text(size=12),
        ,plot.tag.position = c(0, 1)
        )

# # ----- Plot B: Proportion of Blood Meals on Competent Hosts -----
# plot_rho <- all_distributions |>
#   filter(n_patches == 10, infC == 0, de2comp == 1) |>
#   ggplot(aes(x = factor(patch), y = rho)) +
#   geom_point(size = 2.5) +
#   facet_wrap(~scenario1, nrow = 3, as.table = F, drop = F) +
#   labs(x = "Patch", y = "Proportion",
#        title = "Bloodmeal on competent hosts in each patch") +
#   theme_minimal() +
#   theme(plot.title = element_text(face = "bold"),
#         axis.text=element_text(size=10),
#         legend.text=element_text(size=10),
#         strip.text=element_text(size=10),
#         axis.title=element_text(size=12),
#         axis.line = element_line(color = 'black'))

# ----- Plot C: Mosquito Distribution across Patches -----
# How biting vectors (delta) are distributed among patches across different interference levels.
plot_b <- all_distributions |>  
  filter(n_patches == 10, de2comp == 1) |>  
  ggplot(aes(x = factor(infC), y = delta_vect, fill = as.factor(patch))) +
 # geom_bar(stat = "identity", width = 0.85) +
  geom_bar(stat = "identity", width = 0.85) +
  facet_wrap(~scenario1, nrow = 3, as.table = F, drop = F) +
  scale_fill_brewer(palette = "Spectral", direction = 1) +
  scale_x_discrete(labels = agg_labels) +
  labs(x = "Level of interference between mosquitoes",
       y = "Distribution",
       title = "Distribution of mosquitoes across patches",
       fill = "Patch") +
  theme_classic(base_size = 14) +
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=9),
        legend.text = element_text(size=10),
        legend.position = "inside",
        legend.direction = "vertical",
        legend.byrow = T,
        legend.position.inside = c(0.75, 0.9),
        strip.text=element_text(size=12),
        axis.title=element_text(size=12),
        strip.background = element_blank()
        ,plot.tag.position = c(0, 1)
        ) + 
  guides(fill=guide_legend(ncol=5, theme = theme(legend.byrow = T)))

# # ----- Plot C: bites per Competent Hosts -----
plot_c <- all_distributions |>  
  filter(n_patches == 10, de2comp == 1, hostDist != "equal") |>  
  ggplot(aes(x = factor(patch), y = NmPatch*rho/(H_c + epsilon), 
             colour = factor(infC), shape = factor(infC) )) +
  geom_point(size = 3) +
  facet_wrap(~scenario1, nrow = 2, as.table = F) +
  scale_colour_manual(values = cols, name = "", labels = agg_labels3) +
  scale_shape_manual(values = shapes, name = "", labels = agg_labels3) +
  labs(x = "Patch", y = "Bites per competent host",
       title = "Patch-level bites per competent hos") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.25, 0.9),
        legend.background = element_rect(fill = NA, colour = NA),
        axis.line = element_line(color = 'black'),
        axis.text=element_text(size=12),
        legend.text=element_text(size=10),
        strip.text=element_text(size=12),
        axis.title=element_text(size=12),
        ,plot.tag.position = c(1, 1)
        )

plot_d <- all_distributions |>  
  filter(n_patches == 10, de2comp == 1, !(scenario1 == "Equal")) |> 
  group_by(scenario1, infC) %>%
  summarise(
    bites_on_competent = sum(NmPatch * rho),
    bites_on_de = sum(NmPatch * (1-rho)),
    overall_rho        = bites_on_competent / (bites_on_competent + bites_on_de),
    .groups            = "drop") |> 
  ggplot(aes(x = factor(infC), y = overall_rho, fill = factor(infC))) +
  geom_col(position = position_dodge2(), show.legend = F) +
  scale_fill_manual(values = cols) +
  facet_wrap(~scenario1, nrow = 2, as.table = F) +
  theme_minimal(base_size = 14) +
  scale_x_discrete(labels = agg_labels) +
  labs(
    x = "Level of interference between mosquitoes",
    y = "Proportion of bites",
    title = "Overall proportion of bites on competent hosts") +
  theme(axis.text=element_text(size=10),
        strip.text=element_text(size=12),
        axis.title=element_text(size=12)
        ,plot.tag.position = c(1, 1)
        )

# ggsave("./outputs/overallRho.pdf", width = 12, height = 12, units = "in", dpi = 500)

# ----- Arrange Plots using cowplot -----
(plot_a/plot_b | plot_c/plot_d) +
  plot_layout(widths = c(1, 1) ) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

# ggsave("./outputs/10patchComparisonTJ.pdf", width = 12, height = 12, units = "in", dpi = 500)


all_distributions |>  
  filter(infC == 0, de2comp == 1) |>  
  mutate(prop_competent = H_c / (H_c + H_de) ) |>  
  pivot_longer(cols = c(prop_competent, rho)) |>  
  mutate(name = factor(name, levels = c("prop_competent", "rho"), 
                       labels = c("Proportion of competent host x", "Blood meal on competent host x") )) |>  
  ggplot(aes(x = n_patches, y = value, group = patch) ) +
  geom_jitter(size = 1.5, width = 0.2, height = 0.01, alpha = 0.5 ) +
  geom_hline(yintercept = 0.1, colour = "red", linetype = "dashed", linewidth = 0.8) +
  facet_grid(scenario~name) +
  theme_bw() +
  labs(y = "Proportion",x = "Number of patches") +
  theme(strip.background = element_blank())
# ggsave("./outputs/propCompHostVsBloodmeal.pdf", width = 12, height = 12, units = "in")

all_distributions |>  
  filter(de2comp == 1) |> 
  ggplot(aes(x = n_patches, y = NmPatch*rho/(H_c + epsilon) *rho*delta_vect, 
             colour = H_c / (H_c + H_de), group = patch ) ) +
  geom_jitter(size = 1.5, width = 0.25, height = 0.015, alpha = 0.8 ) +
  facet_grid(scenario~mosAgg, scales = "free_y") +
  scale_color_viridis_c(option = "H", name = "Proportion\nof host x") +
  labs(
    y = "Bites per competent host x",
    x = "Number of patches" ) +
  theme_minimal() +
  theme(strip.background = element_blank(),
        text = element_text(size = 12),
        strip.text = element_text(size = 8) )
# ggsave("./outputs/propCompHostVsBloodmealByInfC.pdf", width = 12, height = 12, units = "in")  



###############################################################################
############################# R0 for full spatial model #######################
###############################################################################
extendedR0 <- function( 
    mu_m = 1 / 20,               # Mosquito mortality rate         
    mu_h = 1 / 365,              # Host mortality rate      
    alpha = 1 / 3,               # Mosquito bite rate
    p_hv = 0.5,                  # Probability of host-to-vector transmission
    p_vh = 0.5,                  # Probability of vector-to-host transmission
    sigma = 1 / 5,               # Host recovery rate
    prefComp = 0.1,              # Mosquito preference for competent hosts   
    H_c_vect = 500,              # Vector of competent host numbers
    H_de_vect = 500,             # Vector of dead-end host numbers
    mos_per_comp_host = 20,      # Mosquito density to competent host ratio
    infC = 0                     # Interference constant for host aggregation
    ,hostDefBeh = FALSE
    ,m_c = NA
    ,m_de = NA ){
  
  eps <- 1e-10
  
  H_c_vect  <- as.numeric(H_c_vect)
  H_de_vect <- as.numeric(H_de_vect)
  
  H_tot_vect <- H_c_vect + H_de_vect # Total hosts per patch
  
#  Nh <- sum(Nh_vect)
  
 # propNh_vect <- Nh_vect/Nh # proportion of hosts in each patch
  np <- length(H_tot_vect)          # Number of patches
  
  H_c_total <- sum(H_c_vect)
  H_de_total <- sum(H_de_vect)
  
  Nm <- mos_per_comp_host * H_c_total
  
  # Mosquito distribution
  if (infC == 0) {
    delta_vect <- rep(1/np, np)
  } else {
    propH <- H_tot_vect / sum(H_tot_vect)
    delta <- propH^(1/infC)
    delta_vect <- delta / sum(delta)
  }
  
  # Per-patch rho (bloodmeal on host x)
  if (!hostDefBeh) {
    rho <- prefComp * H_c_vect / (prefComp * H_c_vect + (1 - prefComp) * H_de_vect + eps)
  } else { # Ideal‑free behaviour
    numBiteV <- alpha * Nm / (H_c_vect + eps)
    term  <- ((1 - prefComp) / (prefComp)^(1/m_de) * numBiteV^(m_c/m_de) * H_de_vect + eps)
    rho <- (numBiteV * H_c_vect) / (numBiteV * H_c_vect + term )
  }
  

  term_vec_to_host <- alpha * p_vh * rho * delta_vect / mu_m
  term_host_to_vec <- alpha * p_hv * rho * delta_vect * Nm / ((sigma + mu_h) * H_c_vect + eps)
  
  R0 <- sqrt(sum(term_vec_to_host * term_host_to_vec))
  
  return(R0)
}

#*################################################################
#******************** Baseline scenario **************************
#*################################################################
#* Mosquitoes have equal access to all patches, and all patches 
#* are equally (equal no. of hosts in each patch) suitable 
#* (and no 'long-range' host-seeking behavior)
# a. Within-patch host choice is independent of vector density
plot_heatmap <- expand_grid(prefComp = seq(0, 0.5, length.out = 300), 
            de2Comp = seq(1, 50, length.out = 300),
            n_patches = 20 ) |>  
  rowwise() |> 
  mutate(H_de = de2Comp*startingSh,
         r0 = extendedR0(H_c_vect = rep(startingSh/n_patches, n_patches),
                         H_de_vect = rep(H_de/n_patches, n_patches),
                         prefComp = prefComp,
                         infC = 0
         ) ) |>  
  ggplot(aes(x = prefComp, y = de2Comp, z = r0)) +
  geom_raster(aes(fill = r0) ) +
  geom_contour(color = "gray10", breaks = 1, linewidth = 1) +
  scale_fill_distiller(palette = "Spectral", direction = -1) +
  labs(x = expression("Preference for competent hosts, " * eta[x]),
       y = "Dead-end to competent host ratio",
       fill = expression(R[0])) +
  annotate("text", x = 0.2, y = 20, label = "R[0] < 1", parse = TRUE, size = 8, color = "black") +
  annotate("text", x = 0.4, y = 2.5, label = "R[0] > 1", parse = TRUE, size = 8, color = "black") +
  theme_minimal() +
  theme(axis.text=element_text(size=12),
      legend.text=element_text(size=12),
      axis.title=element_text(size=12))

plot_heatmap
# ggsave("./outputs/simpleR0heatmap.pdf", width = 12, height = 10, units = "in")

# Jen: This is what we compare all other scenarios to? Does this yield the same
# result as if there was only one patch? Just thinking here to make clear that
# this is what is usually assumed in most models (or not)?
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
scenario_2 <- expand_grid(
  prefComp = 0.1
  ,de2Comp = seq(1, 100, length.out = 100)
  ,hostDist = c("exp:hostx", "exp:hosty")
  ,nPatches = seq(1, 50, 1)
  ,infC = c(0, 0.1, 0.9)
  ,decay = c(0.2, 1)
  ) |> 
  mutate(
    scenario = case_when(
    #  hostDist == "equal" ~ "Equal",
      hostDist == "exp:hostx"  ~ paste("Exp: comp. (x), decay =", decay),
      hostDist == "exp:hosty"  ~ paste("Exp: dead-end (y), decay =", decay) 
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
  group_by(de2Comp, nPatches, hostDist, infC, scenario) |> 
  dplyr::summarise(mean_r0 = mean(scenario_r0),
                   
                   mean_rel = mean(rel_change),
                   lower_rel = Rmisc::CI(rel_change)["lower"],
                   upper_rel = Rmisc::CI(rel_change)["upper"],
                   
                   mean_abs = mean(abs_change),
                   lower_abs = Rmisc::CI(abs_change)["lower"],
                   upper_abs = Rmisc::CI(abs_change)["upper"],
                   .groups = "drop" ) 


# no mutual mosquito interference on either hosts Vs no of patches, if deComp = 1
plot_patches <- scenario_res_df |>  
  filter(de2Comp == 1) |>  
  ggplot(aes(x = nPatches, y = mean_abs, ymin = lower_abs, 
             ymax = upper_abs, fill = factor(infC), colour = factor(infC) ) ) +
  geom_ribbon(alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~scenario, scales = "fixed") +
  scale_colour_manual(values = cols, name = "", labels = agg_labels3) +
  scale_fill_manual(values = cols, name = "", labels = agg_labels3) +
  scale_linetype_manual(values = c("solid", "dashed") ) +
  labs(x = "Number of patches", y = expression(Delta~R[0]) ) +
  theme_minimal() +
  plotThemeFunc(leg.pos = "top")

plot_patches
# ggsave("./outputs/R0PatchesbyDistnPref.pdf", width = 12, height = 10, units = "in", dpi = 500)

plot_host_comp <- scenario_res_df |>  
  filter(nPatches == 20) |>  
  ggplot(aes(x = de2Comp, y = mean_abs, ymin = lower_abs, 
             ymax = upper_abs, fill = as.factor(infC), 
             colour = as.factor(infC) )) +
  geom_ribbon(alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~scenario, scales = "fixed") +
  scale_colour_manual(values = cols, name = "", labels = agg_labels3) +
  scale_fill_manual(values = cols, name = "", labels = agg_labels3) +
  scale_linetype_manual(values = c("solid", "dashed") ) +
  labs(x = "Dead-end to competent host ratio", y = expression(Delta~R[0])  ) +
  theme_minimal() +
  scale_x_log10() +
  plotThemeFunc(leg.pos = "top")

plot_host_comp
# ggsave("./outputs/R0HostCompByDistnPref.pdf", width = 12, height = 10, units = "in", dpi = 500)



###############################################################################
#*########################## Sensitivity analysis with R0 #####################
###############################################################################
# Define parameter ranges
lower <- c(
  mu_m = 1/30,            # Mosquito mortality rate (per day)
  mu_h = 1/365,           # Host mortality rate (per day)
  alpha = 1/5,           # bite rate (per day)
  p_hv = 0.001,               # Probability of infection from host to mosquito
  p_vh = 0.001,               # Probability of infection from mosquito to host
  sigma = 1/10,           # Host recovery rate (per day)
  prefComp = 0.001,       # Mosquito preference for competent hosts
  mos_per_comp_host = 1,# Mosquito density per competent host
  infC = 0.001,               # Mosquito aggregation parameter
  de2Comp = 1,            # Dead-end to competent host ratio
  n_patches = 1,
  decay = 0.001
)

upper <- c(
  mu_m = 1/10,            # Mosquito mortality rate (per day)
  mu_h = 10/365,          # Host mortality rate (per day)
  alpha = 1,              # bite rate (per day)
  p_hv = 1,               # Probability of infection from host to mosquito
  p_vh = 1,               # Probability of infection from mosquito to host
  sigma = 1,              # Host recovery rate (per day)
  prefComp = 0.5,         # Mosquito preference for competent hosts
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
  geom_boxplot(width = .6, alpha = 1, outlier.alpha = 0.1, outlier.size = 0.1, notch = T) +
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
