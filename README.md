# MascroscopicStateExperience
This is the code for analysising data for the manuscript of "Linking One-Dimensional Symbolic Dynamics of Resting-State Brain to Recent Experience"
1. Mfiles named started by "Step" is the scripts who build up for hippocampal replay analysis. Important scripts are listed as below.
   (1)Step02*.m for preprocessing task fMRI data.
   (2) processRest.m for preprocessing resting state fMRI data.
   (3)Step03*.m for GLM model for task of memory material processing, and
   (4)Step03*perm.m for GLM in permuted onsets.
   (5)Step06* is for hippocampla replay estimation.
3. Mfiles for whole brain dynamics model estimation as below
   (1)EnergyLandscape8ROIs.m building energy landscape for 8 ROIs.
   (2)Randomwalk.m for simulations on the randomly walk on energy landscape.
   (3)*statetransition* or *trans* stastics on major states frequencies, on real data or on simulation data.
   (4) statsOnSImulationFake.m stastical on perturbed energy lanscape states.
   (5)RegressModel.m regress model to control for age and gender for behaviral data.
