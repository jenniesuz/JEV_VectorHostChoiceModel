
library(deSolve)            

#********************Mosquito-borne virus transmission model***********************************

mod <- function(tt,yy,parms) with(c(parms, as.list(yy)), {
  
  eps <- 1e-10  # Small value to prevent division by zero
  
  #*****************state variables**********************
  Sh = matrix(yy[3:(2+np)], nrow = np, ncol = 1)             # np = number of patches
  Ih = matrix(yy[(3+np):(2+np+np)], nrow = np, ncol = 1)  
  Rh = matrix(yy[(3+np+np):(2+np+np+np)], nrow = np, ncol = 1)
 # Ihcum = matrix(yy[(4+np+np):(3+np+np+np)], nrow = np, ncol = 1)
  
  Nh <- Sh + Ih + Rh  # Total hosts per patch
  
  Nv <- Sv + Iv         # total vectors
  
  # proportion of blood meals on competent hosts:
  if (!hostDefBeh) {
    rho <- prefComp * Nh / (prefComp * Nh + (1 - prefComp) * S_de)  
  } else { # Ideal‑free defensive behaviour
    numBiteV <- alpha * Nv / (Nh + 1e-10)
    term  <- ((1 - prefComp) / prefComp)^(1/m_de) * numBiteV^(m_c/m_de) * S_de
    rho <- (numBiteV * Nh) / (numBiteV * Nh + term )
  }

  #****************Mosquito density in each patch************************
  
  if (infC == 0) { # equally distributed independent of host density
    delta_vect <- rep(1 / np, np) 
  } else { # Sutherland 1986
    propHosts <- (Nh + S_de) / (sum(Nh + S_de) + eps)  # Prevent zero total hosts
    delta_vect <- propHosts^(1 / infC)
    delta_vect <- delta_vect / (sum(delta_vect) )
  }
  
  # mosquito distributed into patches
  SvPatch <- Sv * delta_vect
  IvPatch <- Iv * delta_vect

  #*********************odes*****************************
  
  new_host_infections <- alpha * p_vh * rho * IvPatch * Sh / (Nh + eps)
  
  # mosquitoes
  dSv <- f_m * Nv - sum((alpha * p_hv * rho * Ih / (Nh + eps) * SvPatch)) - mu_v * Sv    
  dIv <- sum((alpha * p_hv * rho * Ih / (Nh + eps) * SvPatch)) - mu_v * Iv       
  
  # hosts
  dSh <- f_h * Nh - new_host_infections - mu_h * Sh           
  dIh <- new_host_infections - phi_h * Ih - mu_h * Ih                  
  dRh <- phi_h * Ih - mu_h * Rh
 # dIhcum <- new_host_infections
  return(list(c(dSv, dIv, dSh, dIh, dRh)))
})

#*******************Function to run model**********************************
sim <- function(init=initial, tseq = times, modFunction=mod
                , parms = params()) {
  simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  return(simDat)
}


simV <- Vectorize(sim)
