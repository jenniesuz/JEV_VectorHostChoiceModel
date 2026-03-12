###############################################################################
############################# R0 for full model #######################
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
    # term  <- ((1 - prefComp) / (prefComp)^(1/m_de) * numBiteV^(m_c/m_de) * H_de_vect + eps) # parenthesis wrongly placed
    term  <- ((1 - prefComp) / prefComp)^(1/m_de) * numBiteV^(m_c/m_de) * H_de_vect
    rho <- (numBiteV * H_c_vect) / (numBiteV * H_c_vect + term )
  }
  
  
  term_vec_to_host <- alpha * p_vh * rho * delta_vect / mu_m
  term_host_to_vec <- alpha * p_hv * rho * delta_vect * Nm / ((sigma + mu_h) * H_c_vect + eps)
  
  R0 <- sqrt(sum(term_vec_to_host * term_host_to_vec))
  
  return(R0)
}