
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
       title = "A) Distribution of hosts across 10 patches") +
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


plot_b <- all_distributions2 |>  
  filter(n_patches == 20, infC == 0, de2comp == 1) |>  
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
       title = "B) Distribution of hosts across 20 patches") +
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


cowplot::plot_grid(plot_a,plot_b, nrow=2)
ggsave("./outputs/SFile10_20patchComparison.pdf", width = 12, height = 12, units = "in", dpi = 500)

