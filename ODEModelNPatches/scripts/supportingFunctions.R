
#************************Distribute hosts into patches*********************
#* Patches can have zero hosts after distribution
distribute_hosts <- function(n = 1, size, prob) {
  prob <- prob / sum(prob)  # Normalize probabilities
  nPatches <- length(prob)
  
  # Distribute all hosts according to probabilities 
  finalDistribution <- rmultinom(n = n, size = size, prob = prob)
  return(as.vector(finalDistribution))
}

#************************SETUP INTIIAL CONDITIONS*********************
initialCondsFunc <- function(nPatches = 2    # number of patches
                             ,totalHosts = 1000 # Fixed total number of hosts
                             ,totalSh = 500  # number of susceptible hosts
                             ,de2comp = 1    # ratio of dead end to susceptible hosts
                             ,hostDist = "equal",
                             decay = NA){
  
  # check hostDist-decay supplied
  stopifnot(hostDist %in% c("equal", "exp:hostx", "exp:hosty"))
  if (hostDist == "equal" && !is.na(decay)) {
    stop("Decay parameter not applicable to equal host distribution. Remove decay argument.")
  }
  if (hostDist %in% c("exp:hostx", "exp:hosty") && (is.na(decay) || decay <= 0)) {
    stop("Decay must be a positive numeric value for exponential host distributions")
  }
  
  #**************State variables************
  Sh = matrix(paste("Sh", 1:nPatches, sep = "_"), nrow = nPatches, ncol = 1)
  Ih = matrix(paste("Ih", 1:nPatches, sep = "_"), nrow = nPatches, ncol = 1)
  Rh = matrix(paste("Rh", 1:nPatches, sep = "_"), nrow = nPatches, ncol = 1)

  #*********** Calculate dead-end hosts ***************
  if (is.na(de2comp)) {
    # Calculate deadEnds from totalHosts if de2comp is not provided
    totalDeadEnds <- totalHosts - totalSh
    if (totalDeadEnds < 0) {
      stop("Error: totalSh cannot be greater than totalHosts.")
    }
  } else {
    # Calculate deadEnds from de2comp 
    totalDeadEnds <- totalSh * de2comp
  }
  
  #***********distribution of hosts***************
  if (hostDist == "equal") {
    # Equal distribution explicitly ignores decay
    deadEnds <- matrix(rep(totalDeadEnds / nPatches, nPatches), nrow = nPatches)
    startingSh <- matrix(rep(totalSh / nPatches, nPatches), nrow = nPatches)
  } else if (hostDist == "exp:hostx") {
    # Exponential decay ONLY for competent hosts (host x)
    deadEnds <- matrix(rep(totalDeadEnds / nPatches, nPatches), nrow = nPatches)
    startingSh <- matrix(
      distribute_hosts(
        n = 1, 
        size = totalSh, 
        prob = rev(exp(-decay * 1:nPatches))
      ), nrow = nPatches)
  } else if (hostDist == "exp:hosty") {
    # Exponential decay ONLY for dead-end hosts (host y)
    deadEnds <- matrix(
      distribute_hosts(
        n = 1, 
        size = totalDeadEnds, 
        prob = rev(exp(-decay * 1:nPatches))
      ), nrow = nPatches)
    startingSh <- matrix(rep(totalSh / nPatches, nPatches), nrow = nPatches)
  }
  
  
  initial <- c(startingSv - 1               # susceptible mosquitoes
               ,1                           # infectious mosquitoes
               ,startingSh                  # susceptible hosts
               ,rep(0, nPatches)       # infected hosts
               , rep(0, nPatches)       # recovered hosts
  )
  names(initial) <- c("Sv","Iv",Sh, Ih, Rh)
  
  return( list(initial, deadEnds) )
}

#*******************************************************************
#*******************Function to run model**********************************
sim <- function(init=initial, tseq = times, modFunction=mod
                , parms = params()) {
  simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  return(simDat)
}


simV <- Vectorize(sim)
#**************************************************************************



outputsFunc <- function(modelOutput){
  
  if(is.data.frame(modelOutput[,grep("Sh",names(modelOutput))])){
    totSh <- rowSums(modelOutput[,grep("Sh",names(modelOutput))])
    totIh <- rowSums(modelOutput[,grep("Ih",names(modelOutput))])
    totRh <- rowSums(modelOutput[,grep("Rh",names(modelOutput))])
  }else{
    totSh <- modelOutput$Sh
    totIh <- modelOutput$Ih
    totRh <- modelOutput$Rh
  }
  N <-  rowSums(modelOutput[,c(grep("Sh",names(modelOutput)), grep("Ih",names(modelOutput)) , grep("Rh",names(modelOutput)) ) ])
  
  hostProps <- cbind.data.frame(time=rep(modelOutput$time,3)
                                ,class=c(rep("Susceptible",length(modelOutput[,1]))
                                         ,rep("Infected",length(modelOutput[,1]))
                                         ,rep("Resistant",length(modelOutput[,1])))
                                ,value=c(totSh/N*100
                                         ,totIh/N*100
                                         ,totRh/N*100)
  )
  
  mozProps <- cbind.data.frame(time=modelOutput$time
                               ,pIm1000=(modelOutput$Iv)/(modelOutput$Iv + modelOutput$Sv)*1000
                               ,mph=rowSums(modelOutput[,2:3])/N
  )
  
  
  return(list(hostProps,mozProps))
}


leg.pos<-c(0.5,0.5)

# plotThemeFunc <-  function(leg.pos){
#  return(theme_set(theme_bw()) +
#   theme(
#     axis.line = element_line(color = 'black')
#     ,text=element_text(size=6)
#     ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
#     ,axis.text=element_text(size=5)
#     ,legend.key.size = unit(0.8,"line")
#     ,legend.background = element_blank()
#     ,legend.text=element_text(size=4)
#    # ,legend.position =leg.pos
#     ,legend.title=element_text(size=4)
#     ,strip.background = element_rect(colour="white", fill="white")
#     ,strip.text.x = element_text(size = 4)
#     ,panel.border = element_blank()
#   )
#  )
# }


plotThemeFunc <-  function(leg.pos){
  return(theme_set(theme_bw()) +
           theme(
             axis.line = element_line(color = 'black')
             ,text=element_text(size=12)
             ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
             ,axis.text=element_text(size=12)
             ,legend.key.size = unit(0.8,"line")
             ,legend.background = element_blank()
             ,legend.text=element_text(size=12)
             ,strip.text=element_text(size=12)
             ,axis.title = element_text(size = 14)
             ,legend.position =leg.pos
             ,legend.title=element_text(size=12)
             ,strip.background = element_rect(colour="white", fill="white")
             ,strip.text.x = element_text(size = 12)
             ,panel.border = element_blank()
           )
  )
}


compute_stats <- function(x) {
  ci_x <- Rmisc::CI(x)
  c(
    mean = mean(x),
    lower = ci_x[["lower"]],
    upper = ci_x[["upper"]]
  )
}


# Function to refactor host distributions for plotting
refactor_host_dist <- function(df) {
  # Get unique sorted decays (numeric sort)
  decays <- df |> 
    filter(hostDist %in% c("exp:hostx", "exp:hosty")) |> 
    distinct(decay) |> 
    pull(decay) |> 
    sort()
  
  # Create ordered levels with concise labels
  level_order <- c(
    paste("Exp: dead-end (y), decay =", decays),
    paste("Exp: comp. (x), decay =", decays),
    "Equal"
  )
  
  df |> 
    mutate(
      scenario1 = case_when(
        hostDist == "equal"    ~ "Equal",
        hostDist == "exp:hostx"  ~ paste("Exp: comp. (x), decay =", decay),
        hostDist == "exp:hosty"  ~ paste("Exp: dead-end (y), decay =", decay)
      ),
      scenario1 = factor(scenario1, levels = level_order)
    ) |> 
    arrange(scenario1)
}





generate_distribution_data <- function(
    n_patches, totalSh, de2comp = 1, hostDist, infC, prefComp = 0.1, Nm, decay = 1,
    hostDefBeh = FALSE, m_c = NA, m_de = NA
){
  eps <- 1e-10
  # Generate host distributions
  ic <- initialCondsFunc(
    nPatches = n_patches,
    totalSh = totalSh,
    de2comp = de2comp,
    hostDist = hostDist,
    decay = decay
  )
  
  H_c <- unname(ic[[1]][startsWith(names(ic[[1]]), 'Sh')])
  H_de <- as.vector(ic[[2]])
  
  if (infC == 0) { # equally distributed independent of host density
    delta_vect <- rep(1/n_patches, n_patches) 
  } else { # Sutherland 1986
    propHosts <- (H_c + H_de)/sum(H_c, H_de) # proportion of total hosts in each patch
    delta_vect <- propHosts^(1/infC)
    delta_vect <- delta_vect/sum(delta_vect)
  }
  
  # proportion of blood meals on competent hosts:
  if(hostDefBeh) {
    numBiteV <- (1/3) * Nm / (H_c + eps) 
    rho_x <- (numBiteV * H_c) / (numBiteV * H_c + 
                                   ((1-prefComp)/prefComp)^(1/m_de) * numBiteV^(m_c/m_de) * H_de + eps)
  } else {
    rho_x <- prefComp * H_c / (prefComp * H_c + (1 - prefComp) * H_de + eps)
  }
  
  return(
    tibble(
      patch = 1:n_patches,
      H_c = H_c,
      H_de = H_de,
      delta_vect = delta_vect,
      NmPatch = Nm*delta_vect, # number of mosquitoes in each patch
      rho = rho_x
    )
  )
}
