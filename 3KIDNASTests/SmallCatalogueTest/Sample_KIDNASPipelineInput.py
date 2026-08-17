#           This file contains definitions needed for the WRKP_CatalogueFitDriver
#
#   First give the name of the Catalogue that will be fit
CatName="../TestData/WALLABY_Test_Catalogue.fits"
#   Next give the name of the source folder
SourceFolder="../TestData/WALLABY_Test_sources/"
#   And give a name for the folder that will contain the fitting results
TargFolder="WALLABY_Test_Models/"
#   Also set the number of bootstraps to make and run
nBootstraps= 50
#   Set the number of processors to use for different galaxies
nProcessors=5
#   And the number of processors to use for bootstraps for each individual galaxy
nProcessors_Bootstrap=5
#   Set the KinTR keyword -- the kinematic team release version
KinTR="WALLABY_TR1"
#   Set the SRCTR keyword -- the team release for the sources
SRCTR="WALLABY"
