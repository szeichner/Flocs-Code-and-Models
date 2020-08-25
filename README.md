# Flocs-Code

This repository contains 3 distinct pieces of code: 
(1) Matlab code to process images from settling experiments, and convert them into settling velocities based on calibration data. The file "FlocImageAnalysis.m" processes photos from a designated folder where the photos are stored; "FlocDataAnalysis.m" uses exported matrices from the image analysis to compute least squares regressions and calculate settling velocities for the experiments.
(2) R code to perform t-tests on distinct sets of experiments to test for statistical significance
(3) Matlab and R code to model floodplain grain size specific deposition rates for flocculated and un-flocculated mud cases. The file "parametric_gsd.R" generates a log-normal grain size distribution to input into the model; "floodplain_profile_model.m" models grain size specific deposit thickness; "floodplain_average_mud_deposition_model.m" models floodplain-averaged mud deposition rate.
