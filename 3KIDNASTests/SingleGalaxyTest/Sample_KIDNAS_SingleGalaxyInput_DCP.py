#           This file contains definitions needed for the WRKP_GalaxyFitDriver
#
#   First give the name of the cube that will be fit
CubeName="../TestData/WALLABY_Test_sources/WALLABY_J103538-484832/WALLABY_J103538-484832_VelCube.fits"
#   Next give the name of the base cube mask
MaskName="../TestData/WALLABY_Test_sources/WALLABY_J103538-484832/WALLABY_J103538-484832_mask.fits"
#   Also give a name for the object for naming purposes
ObjName="WALLABY_J103538-484832"
#   And give a name for the folder that will contain the fitting results
TargFolder="TestFits/"
#   And get the initial estimates for the inclination and position angle in degrees
PA_Estimate= 241.319
Inc_Estimate=89.00
#   Also set the number of bootstraps to make and run
nBootstraps= 500
#   And set the number of processors to use for the bootstraps
nProcessors_Bootstraps=8

#   ---- DCP bootstrap distribution ----
#   FixtureFile/DCPJobDir/DCPBootstrapLauncher are fixed locations under
#   js/ (see SetFileLocations.py).
#   Compute group and slice price are optional and come from the environment
#   at run time (DCP_COMPUTE_GROUPS, DCP_SLICE_PRICE) -- see RunBootstrapsDCP.py.
UseDCP = 1