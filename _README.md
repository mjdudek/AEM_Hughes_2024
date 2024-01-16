README / AtmosphericEffectsModel_Hughes_2024

------------------------------------------------------------------------------------------------------------------------------
AUTHOR(S): Marissa Hughes(1)*, R Shane McGary (2)

AFFILIATION(S): 
(1) University of North Carolina at Chapel Hill
(2) James Madison University
* Corresponding author, mjdudek@email.unc.edu

LANGAUGE: MATLAB (2020a), RStudio (R), Microsoft Excel (2021)

ABOUT: A repository of MATLAB scripts, workspaces, and data files associated with the Atmospheric Effect Model (Hughes, 2024). 
All files are .m format accesible through MATLAB 2020a and MATLAB 2022a. 

------------------------------------------------------------------------------------------------------------------------------
REPOSITORY DESCRIPTION

ATMOSPHERICEFFECTSMODEL_HUGHES_2024
/ _JGR__AEM_Manuscript_240110.pdf
/ _README
/ AEM_240108.m
/ FLAVRS_240105.m

/ FIGURE_CREATION
Contains a MATLAB script, workspace, and powerpoint file associated with figures created for the manuscript
.. / AEM_ManuscriptFigures / .PNGs of figures
.. / AEM_ManuscriptFigures.ppt
.. / atmosprofiles.m
.. / projevol.m
.. / terrestrialanalogs.m
.. / AEM_HeatMaps_v2.R
.. / barcharts.exw

/ MATLAB_WORKSPACES (available upon request)
Contains MATLAB workspaces that can be imported to MATLAB with model outputs and processing of the theoretical and 
terrestrial projectile population model runs. Import workspaces to see data without running the model itself. 
.. / Workspace_AEM_TerrestrialAnalogs_AblationFragmentation.m
.. / Workspace_AEM_TerrestrialAnalogs_Ablationonly.m
.. / Workspace_AEM_TheoreticalPopulation_AblationFragmentation.m
.. / Workspace_AEM_TheoreticalPopulation_AblationOnly.m
.. / Workspace_AEM_TheoreticalPopulation_Abridged_AblationFragmentation.m
.. / Workspace_AEM_TheoreticalPopulation_Abridged_AblationOnly.m

/ PROJECTILE_POPULATIONS
Contains FLAVRS matricies created using the FLAVRS script for the theoretical and terrestrial projectile populations
.. / FLAVRS_TerrestrialAnalogPopulation.m
.. / FLAVRS_TheoreticalPopulation.m
.. / FLAVRS_TheoreticalPopulation_Abridged.m

------------------------------------------------------------------------------------------------------------------------------
FIGURE CREATION

FIG1 - AEM Diagram
Created using PowerPoint

FIG2 - Planetary Atmospheric Conditions
(script) atmosprofiles.mat
(input) Workspace_AEM_TerrestrialAnalogs_AblationFragmentation.m

FIG3 - AEM Flowchart
Created used lucid.io website

FIG4 - Model Sensitivity to Atmospheric Effect Equations
...... Note: use options 1 for flag selection
(script) AEM_240108.m
(input) FLAVRS_MeteorCraterAnalog.m

FIG5 - Model Sensitivity to Atmospheric Effect Paramerts
...... Note: use options 1 for flag selection
(script) AEM_240108.m
(input) FLAVRS_MeteorCraterAnalog.m

FIG6 - Velocity & Mass Evolutions of Terrestrial Analogs
(script) terrestrialanalogs.m
(input) Workspace_AEM_TerrestrialAnalogs_AblationOnly.m
(input) Workspace_AEM_TerrestrialAnalogs_AblationFragmentation.m

FIG7 - Dispersion of Terrestrial Analogs
(script) terrestrialanalogs.m
(input) Workspace_AEM_TerrestrialAnalogs_AblationFragmentation.m

FIG8 - Initial Projectile Mass vs. Final Crater Diameter
(script) projevol.m
(input) Workspace_AEM_TheoreticalPopulation_AblationFragmentation.m
(input) Workspace_AEM_TheoreticalPopulation_Expanded_AblationFragmentation.m
(input) Workspace_AEM_TerrestrialAnalogs_AblationFragmentation.m

FIG9 - Theoretical Projectile Evolution with Ablation & Fragmentation
(script) projevol.m
(input) Workspace_AEM_TheoreticalPopulation_Abridged_AblationFragmentation.m

FIG10 - Initial Characteristics of Theoretical Impactors
(script) barcharts.exw
(input) ImpactStats_v2_empirical.csv

FIG11 - Projectile Characteristics vs. Final Crater Diameter
(script) AEM_HeatMaps_v2.R
(input) statsE_ablation.csv
(input) statsE_ablation-break.csv
(input) statsL_ablation.csv
(input) statsL_ablation-break.csv
(input) statsM_ablation.csv
(input) statsM_ablation-break.csv
(input) statsV_ablation.csv
(input) statsV_ablation-break.csv

