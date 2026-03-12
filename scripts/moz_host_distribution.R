
rm(list = ls())
# packages required
library(reshape)
library(patchwork)
library(parallel)
library(here)
library(parallel)
library(tidyverse)
library(ggplot2)

source(here("scripts/supportingFunctions.R"))

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

###############################################################################
####################### Plots host distribution in patches ####################
###############################################################################

# fixed parameters
hostPreference <- 0.1 # mosquito preference for competent hosts
startingSh <- 500       # Number of susceptible hosts
startingSv <- 10000 

# Introduce a small constant to avoid division by zero.
epsilon <- 1e-10

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
            decay = parameter_grid$decay[i],
            prefComp = hostPreference
          )
    )
  }, mc.cores = mc.cores )
) 


# Apply to data
all_distributions2 <- all_distributions |>  
  refactor_host_dist()



# ----- Plot A: Host Distribution per Patch -----
# This plot shows the number of hosts in each patch, distinguishing competent (H_c) and dead-end (H_de) hosts.
plot_a <- all_distributions2 |>  
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
       title = "A) Distribution of hosts across patches") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.9),
        legend.background = element_rect(fill = NA, colour = NA),
        axis.line = element_line(color = 'black'),
        axis.text=element_text(size=11),
        legend.text=element_text(size=11),
        strip.text=element_text(size=11),
        axis.title=element_text(size=11),
        plot.tag.position = c(0, 1)
  )

# ----- Plot B: Mosquito Distribution across Patches -----
# How biting vectors (delta) are distributed among patches across different interference levels.
plot_b <- all_distributions2 |>  
  filter(n_patches == 10, de2comp == 1) |>  
  ggplot(aes(x = factor(infC,levels=c("0","0.9","0.1")), y = delta_vect, fill = as.factor(patch))) +
  # geom_bar(stat = "identity", width = 0.85) +
  geom_bar(stat = "identity", width = 0.85) +
  facet_wrap(~scenario1, nrow = 3, as.table = F, drop = F) +
  scale_fill_brewer(palette = "Spectral", direction = 1) +
  scale_x_discrete(labels = agg_labels) +
  labs(x = "Level of interference between mosquitoes",
       y = "Distribution",
       title = "B) Distribution of mosquitoes across patches",
       fill = "Patch") +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=9),
        legend.text = element_text(size=10),
        legend.position = "inside",
        legend.direction = "vertical",
        legend.byrow = T,
        legend.position.inside = c(0.75, 0.9),
        strip.text=element_text(size=11),
        axis.title=element_text(size=11),
        strip.background = element_blank(),
        plot.tag.position = c(0, 1)
  ) + 
  guides(fill=guide_legend(ncol=5, theme = theme(legend.byrow = T)))

# # ----- Plot C: bites per Competent Hosts -----
plot_c <- all_distributions2 |>  
  filter(n_patches == 10, de2comp == 1, hostDist != "equal") |>  
  ggplot(aes(x = factor(patch), y = NmPatch*rho/(H_c + epsilon), 
             colour = factor(infC,levels=c("0","0.9","0.1")), shape = factor(infC,levels=c("0","0.9","0.1")) )) +
  geom_point(size = 3) +
  facet_wrap(~scenario1, nrow = 2, as.table = F) +
  scale_colour_manual(values = cols, name = "", labels = agg_labels3) +
  scale_shape_manual(values = shapes, name = "", labels = agg_labels3) +
  labs(x = "Patch", y = "Bites per competent host",
       title = "C) Patch-level bites per competent host") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.25, 0.9),
        legend.background = element_rect(fill = NA, colour = NA),
        axis.line = element_line(color = 'black'),
        axis.text=element_text(size=11),
        legend.text=element_text(size=10),
        strip.text=element_text(size=11),
        axis.title=element_text(size=11),
        plot.tag.position = c(1, 1)
  )

#*------Plot D: overall proportion of bloodmeals on competent hosts-----
plot_d <- all_distributions2 |>  
  filter(n_patches == 10, de2comp == 1, !(scenario1 == "Baseline: hosts equally distributed")) |> 
  group_by(scenario1, infC) %>%
  summarise(
    bites_on_competent = sum(NmPatch * rho),
    bites_on_de = sum(NmPatch * (1-rho)),
    overall_rho        = bites_on_competent / (bites_on_competent + bites_on_de),
    .groups            = "drop") |> 
  ggplot(aes(x = factor(infC,levels=c("0","0.9","0.1")), y = overall_rho, fill = factor(infC,levels=c("0","0.9","0.1")))) +
  geom_col(position = position_dodge2(), show.legend = F) +
  scale_fill_manual(values = cols) +
  facet_wrap(~scenario1, nrow = 2, as.table = F) +
  theme_minimal(base_size = 12) +
  scale_x_discrete(labels = agg_labels) +
  labs(
    x = "Level of interference between mosquitoes",
    y = "Proportion of bites",
    title = "D) Overall proportion of bites on competent hosts") +
  theme(axis.text=element_text(size=10),
        strip.text=element_text(size=11),
        axis.title=element_text(size=11)
        ,plot.tag.position = c(1, 1)
  )


cowplot::plot_grid(plot_a,plot_b,plot_c,plot_d, nrow=2)

# ggsave("./outputs/10patchComparison.pdf", width = 12, height = 12, units = "in", dpi = 500)


# all_distributions2 |>  
#   filter(infC == 0, de2comp == 1) |>  
#   mutate(prop_competent = H_c / (H_c + H_de) ) |>  
#   pivot_longer(cols = c(prop_competent, rho)) |>  
#   mutate(name = factor(name, levels = c("prop_competent", "rho"), 
#                        labels = c("Proportion of competent host x", "Blood meal on competent host x") )) |>  
#   ggplot(aes(x = n_patches, y = value, group = patch) ) +
#   geom_jitter(size = 1.5, width = 0.2, height = 0.01, alpha = 0.5 ) +
#   geom_hline(yintercept = 0.1, colour = "red", linetype = "dashed", linewidth = 0.8) +
#   facet_grid(scenario1~name) +
#   theme_bw() +
#   labs(y = "Proportion",x = "Number of patches") +
#   theme(strip.background = element_blank())
# ggsave("./outputs/propCompHostVsBloodmeal.pdf", width = 12, height = 12, units = "in")

# all_distributions2 |>  
#   filter(de2comp == 1) |> 
#   ggplot(aes(x = n_patches, y = NmPatch*rho/(H_c + epsilon) *rho*delta_vect, 
#              colour = H_c / (H_c + H_de), group = patch ) ) +
#   geom_jitter(size = 1.5, width = 0.25, height = 0.015, alpha = 0.8 ) +
#   facet_grid(scenario1~infC, scales = "free_y") +
#   scale_color_viridis_c(option = "H", name = "Proportion\nof host x") +
#   labs(
#     y = "Bites per competent host x",
#     x = "Number of patches" ) +
#   theme_minimal() +
#   theme(strip.background = element_blank(),
#         text = element_text(size = 12),
#         strip.text = element_text(size = 8) )
# ggsave("./outputs/propCompHostVsBloodmealByInfC.pdf", width = 12, height = 12, units = "in")  

