# The impacts of host community composition, structuring and vector host choice on zoonotic mosquito-borne virus transmission

**Jennifer S. Lord** ([ORCID](https://orcid.org/0000-0001-6616-1526))<sup>1, *</sup>, **Michael B. Bonsall** ([ORCID](https://orcid.org/0000-0003-0250-0423))<sup>2</sup>, **Tijani A. Sulaimon** ([ORCID](https://orcid.org/0000-0002-9735-8548))<sup>1, *</sup>

<sup>1</sup> Department of Vector Biology, Liverpool School of Tropical Medicine, Liverpool, United Kingdom 
<sup>2</sup> Department of Biology, University of Oxford, Oxford, United Kingdom  

**Correspondence**: [jennifer.suzanne.lord@gmail.com](mailto:jennifer.suzanne.lord@gmail.com)

doi: tbc


## Abstract
Zoonotic mosquito-borne viruses, including West Nile, Rift Valley fever and Japanese encephalitis viruses, present a substantial public health burden in some countries, with transmission involving complex interactions between multiple vector and host species. While host community composition is known to affect transmission, the role of host community structure coupled with mosquito aggregation behaviour remains poorly understood. 

Motivated by furthering understanding of Japanese encephalitis virus transmission and persistence in contexts where dead-end hosts are dominant in the community, we developed a patch-based transmission model that integrates mosquito host-seeking across two scales: long-range aggregation toward host-dense patches and short-range host species choice determined by local host availability and innate vector preference. Using this model, we explored the potential effect of host clustering, patch fragmentation, and the strength of mosquito aggregation on the basic reproduction number and epidemic dynamics of zoonotic mosquito-borne viruses.

Our results show that identical host communities can yield very different transmission outcomes depending on how they are distributed and the level of mosquito aggregation. Under circumstances where dead-end hosts are preferred by vectors, if competent hosts are clustered and mosquitoes aggregate strongly, transmission may be amplified compared with uniformly distributed hosts and vectors. If dead-end hosts are clustered, transmission may be diluted as mosquitoes aggregate in patches where dead-end hosts are dominant in the host community. Fragmentation effects are context-dependent: under moderate aggregation, increasing patch number generally increases transmission, while under strong or absent aggregation, the effect reverses depending on which host type is clustered.

%Our findings demonstrate that ignoring poor mixing assumptions may underestimate or overestimate risk in landscapes where host communities are fragmented or unevenly distributed.
These findings highlight that variations in the local spatial context of host-vector interactions can lead to variations in disease risk. Locally-informed surveillance and targeted interventions may therefore be important for the effective control of zoonotic mosquito-borne viruses.

This paper has been submitted to PLoS Neglected Tropical Diseases.

## Contents
- `scripts/`: contains all code for analyzing the model. Brief descriptions of each file:
  - **moz_host_distribution.R**: demonstrates how host and vector distribution influences bites on competent hosts.

  - **R0.R**: produces the plots to demonstrate how host and vector distributions vary and the consequences for R0. 
  - **R0_function.R**: contains the R0 function.
  
  - **model.R**: the functions needed to carry out numerical simulations of the ordinary differential equation model.

  - **parameters.R**: the parameter values used in numerical simulations of the dynamic model.

  - **simulate.R**: runs the dynamic model.
  
  - **supportingFunctions.R**: additional functions needed to run the model.
  
  - **sutherland.R**: script to fit the aggregation function (Equation 1 in the main text) to data from Bangladesh (see Supplementary File 1). The data is contained within **tritaeHostsHH.csv**.

## Getting Started
To run the model:  
1. Open `simulate.R`.  
2. Source the model functions from `model.R`
3. Ensure all required libraries are installed

```R
PASTE R SESSION INFO HERE
