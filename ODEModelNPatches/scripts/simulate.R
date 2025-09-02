# packages required
library(tidyverse)
library(reshape)
library(deSolve)
library(grid)
library(gridExtra)
library(here)

# read in scripts
source(here(".//scripts//model.R"))
source(here(".//scripts//parameters.R"))
source(here(".//scripts//supportingFunctions.R"))

##################################################################
#****************Run different scenarios**********************
#* 
set.seed(202510)

times <- seq(0, 2*365, 1) # Time steps for reporting
nPatches <- 20          # Number of patches

birdPrefs <- 0.1  # assume single 10% preference for competent hosts

#*************Hosts Uniformly distributed*******************
combs <- expand_grid(
  de2comp = 1
  ,nPatches = 20
  ,hostDist = c("equal", "exp:hostx", "exp:hosty")
  ,infC = c(0, 0.1, 0.9)
  ,decay = c(0.2, 1)
  ,rep = 1:100 ) |> 
  mutate(decay = ifelse(hostDist == "equal", NA_real_, decay),
        # rep = ifelse(hostDist == "equal", 1, rep),
         infC = ifelse(hostDist == "equal", 0, infC)
         ) |> 
  distinct() |> 
  mutate(scenario = case_when(
    hostDist == "equal" ~ "Equal",
    hostDist == "exp:hostx"  ~ paste("Exp: comp. (x), decay =", decay),
    hostDist == "exp:hosty"  ~ paste("Exp: dead-end (y), decay =", decay) 
  ),
  mosAgg = paste("Interference = ", infC)
  ) |> refactor_host_dist()

combs

sims <- mclapply(1:nrow(combs), function(x){
  
  ic <- initialCondsFunc(
    nPatches = combs$nPatches[x]
    ,totalSh = startingSh
    ,de2comp = combs$de2comp[x]   # ratio of dead-end to competent hosts
    ,hostDist = combs$hostDist[x]
    ,decay = combs$decay[x]
  )
  
  simDat <- as.data.frame(
    lsoda(y = ic[[1]],
          times = times,
          func = mod,
          parms = params(prefComp = birdPrefs,
                         np = combs$nPatches[x],
                         S_de = as.vector(ic[[2]]),
                         infC = combs$infC[x] )
          )
    ) 
  
  out <- simDat |> 
    pivot_longer(!c(time, Sv, Iv), 
                 names_to = c("compartment", "patch"),
                 names_sep = "_"
                 ) |> 
    group_by(time, Sv, Iv, compartment) |>
    summarise(value = sum(value), .groups = "drop") |>  
    pivot_wider(names_from = compartment,
                values_from = value)
  
  
  return( cbind(combs[x, ] ,out) )                 
         
}, mc.cores = detectCores() - 2 )

dyn <- bind_rows(sims)

#***************************Hosts Exp distributed*******************

dyn_summary <- dyn |> 
  mutate(propIh = Ih/(Sh + Ih + Rh)) |> 
  group_by(scenario1, infC, time) |> 
  summarise(meanIh = mean(propIh),
            lower = Rmisc::CI(propIh)[[1]],
            upper = Rmisc::CI(propIh)[[3]],
            .groups = "drop" )

plot_dynamics <- dyn_summary |> 
  filter(time <= 130) |> 
  ggplot(aes(x = time, y = meanIh, colour = as.factor(infC) )) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = as.factor(infC) ),
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8 ) +
  facet_wrap(~ scenario1, ncol = 2, scales = "fixed", as.table = F, drop = F ) +
  scale_colour_manual(values = cols, name = "", labels = agg_labels3) +
  scale_fill_manual(values = cols, name = "", labels = agg_labels3) +
  scale_y_continuous(n.breaks = 3) +
  labs(x = "Time (days)", y = expression(I[x]/H[x]) ) +
  theme_minimal() +
  plotThemeFunc(leg.pos = c(0.75, 0.85))

plot_dynamics
# ggsave("./outputs/dynByHostDistInfC.pdf", width = 12, height = 12, units = "in", dpi = 500)

dyn |>
  filter(time <= 250) |>
  ggplot(aes(x=time, y = Ih, colour = as.factor(infC), group = interaction(rep, infC) )) +
  geom_line(linewidth = 0.3, alpha = 1 ) +
  facet_wrap(~ paste(hostDist, decay), ncol = 2, scales = "fixed" ) +
  facet_wrap(~ scenario1, ncol = 2, scales = "fixed", as.table = F, drop = F ) +
  scale_colour_manual(values = cols, name = Interference~m[p]) +
  theme_minimal()


# ( (plot_heatmap | plot_host_comp) /
#     (plot_patches | plot_dynamics) ) +
#   plot_layout(
#     heights = c(1, 1), # Equal row heights
#     widths = c(1, 1),  # Equal column widths
#     guides = "collect" # Merge legends
#   ) +
#   plot_annotation(
#     tag_levels = "A",
#  #   tag_suffix = ")",
#     theme = theme(plot.margin = margin(10, 10, 10, 10))
#   ) &
#   theme(
#     plot.tag = element_text(face = "bold", size = 16),
#     plot.tag.position = "top", # Position tags at top-left
#     legend.position = "bottom",
#     legend.box = "horizontal")

( (plot_host_comp/plot_patches) | plot_dynamics ) +
  plot_layout(
    heights = c(1, 1), # Equal row heights
    widths = c(1, 1),  # Equal column widths
    guides = "collect" # Merge legends
  ) +
  plot_annotation(
    tag_levels = "A", tag_suffix = ")", tag_prefix = "(",
    theme = theme(plot.margin = margin(6, 6, 6))
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = "topleft", # Position tags at top-left
    plot.tag.location = "margin",
    legend.position = "bottom")


# ggsave("./outputs/combined_plt.pdf", width = 12, height = 12, units = "in", dpi = 500)

### copy files to LaTeX directory. Clumsy but works for now.
latex_dir   <- "../plots/"
plot_dir <- "./outputs"
plot_files <- list.files(plot_dir, pattern = "\\.pdf$", full.names = TRUE)

file.copy(plot_files, file.path(latex_dir, basename(plot_files)), overwrite = TRUE)

