
#****************************PARAMETERS***************************
params <- function(          # Parameters used in the model
  f_m=fecundM                # mosquito reproduction rate
  ,f_h=fecundH               # host reproduction rate
  ,mu_v = muMoz              # mosquito mortality rate        
  ,mu_h = muHost             # host mortality rate      
  ,alpha = biteRate          # mosquito bite rate
  ,p_hv = phostVec           # probability of host to vector transmission
  ,p_vh = pVecHost           # probability of vector to host transmission
  ,phi_h = recRate           # host recovery rate
  ,prefComp = 0.1              # mosquito preference for competence hosts   
  ,S_de=deadEnds              # numbers of dead-end hosts 
  ,np=nPatches               # number of patches
  ,infC=0.6                  # interference constant for "hostDepAgg"
  ,hostDefBeh = F
  ,m_c = NA
  ,m_de = NA
)
  return(as.list(environment()))
#******************************************************************

biteRate <- 1/3
muMoz <- 1/20              
muHost <- 1/365  
fecundH <- muHost
fecundM <- muMoz
phostVec <- 0.5              
pVecHost <- 0.5                 
recRate <- 1/5               


startingSh <- 500       # Number of susceptible hosts
startingSv <- 10000     # Number of susceptible vectors
mosquitoes_per_comp_host <- startingSv/startingSh

# Introduce a small constant to avoid division by zero.
epsilon <- 1e-10

