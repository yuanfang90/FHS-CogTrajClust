# FHS-CogTrajClust
Code repository for the manuscript "Modeling heterogeneity in cognitive trajectories in the Framingham Heart Study"

- The folder R contains:
  + The code for model specification to fit piecewise linear latent class mixed effect models (LCMM) using the function hlme from the lcmm package.
  + The code for performing backward selection of prespecified change points based on the model results.
- The folde MockData contains:
  + Simulated toy data following the same structure of the datasets analyzed in the manuscript
  + The data used in the manuscript are available from the Framingham Heart Study through NHLBI’s Biologic Specimen and Data Repository Information Coordinating Center (BioLINCC). https://biolincc.nhlbi.nih.gov/studies/fhs/. 
