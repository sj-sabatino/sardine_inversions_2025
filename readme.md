# Sardine Genomics – Supporting Scripts and Analyses



Analysis scripts and supporting code used in the sardine genomics study : citation: [🔗](https://www.biorxiv.org/content/10.1101/2025.01.15.633256v1).  
The code here reflects the actual analyses conducted and is made available for transparency and reproducibility.  

  


---

## repository map
   
   
| location                 |  function                                             
|--------------------------------------------------------------------------------|  
| 🛠️`1_filtering_sync/`    |  Raw variant filtering and sync file prep           |  
| 🛠️`2_pcadapt/`           |  Population structure analysis using pcadapt        |  
| 🛠️`3_inversions/`        |  Analysis of putative chromosomal inversions        |  
| 🛠️`4_dadi_modeling/`     |  Demographic modeling using diffusion models (∂a∂i) |  
| 🛠️`5_ancova_boot/`       |  Bootstrapped ANCOVA analysis for trait association |  
| 🛠️`6_tajimas_pi/`        |  Sliding-window diversity metrics (Tajima's D, π)   |  
|--------------------------------------------------------------------------------|  

## notes

- Scripts are structured by analysis block rather than modular pipelines.  
- Dependencies are noted in header comments or `requirements.txt` where available.  
- Input data is not included due to size. Please contact the author if needed.  

---  

## usage  

Each folder contains independent R/Python/bash scripts. Example:  
For details on their usage please see the preprint on 🔗[bioarxiv](https://www.biorxiv.org/content/10.1101/2025.01.15.633256v1)  



This code is released for academic use. Please cite appropriately.
