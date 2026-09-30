# 3KIDNAS
 The full pipeline for the WALLABY Kinematics analysis


---
Requirements

There are a number of standard python packages used in 3KIDNAS -- numpy, scipy, astropy, matplotlib, multiprocessing.  
Additionally, the colourspace and CosmosCanvas packages are stored in third_party/

The code currently is fixed to use python3.9 (this will be adjusted in the future). (note: this doc's own JS/DCP comparison tooling -- see "JS/DCP pipeline" below -- has been run successfully against this same Fortran/Python pipeline under Python 3.12 in this environment; python3.9 was not independently available to verify against directly)



---
Quick Installation Guide

The first step is to compile the key 3rd party software -- the easiest approach is to go to terminal then

1) cd third_party/fftw-3.3.8/  
    make clean  
    ./configure  
    make  
2) cd third_party/cfitsio-4.6.3/  
    make clean  
    ./configure  
    make  
3) cd third_party/SoFiA-2-master_2_5_1/  
    make clean  
    ./compile.sh  
    
Now it should be possible to compile the main code

1) cd src/  
    make clean  
    make  
    
---
Running and testing the code

There is a set of test data and sample inputs currently available at:
https://www.dropbox.com/scl/fi/i3elzk37u2dd53mrn0tlc/3KIDNAS_SampleRuns.zip?rlkey=qub81w5wcfqbf968e8y8ubv7r&st=etunwr2i&dl=0

To run on a single galaxy go to the folder with the GalaxyInputFile.py and run:  
python '$(PATH)/WRKP_GalaxyFitDriver.py $(GalaxyInputFile)'

To run on a catalogue:  
python '$(PATH)/WRKP_CatalogueDriver.py $(CatalogueInputFile)'  
(note: WRKP_CatalogueDriver.py does not currently exist in this repo -- only WRKP_GalaxyFitDriver.py is present, plus a 3KIDNASTests/SmallCatalogueTest/ folder with sample catalogue input but no driver script to run it)

---
JS/DCP pipeline

Alongside the native Fortran/Python pipeline above, js/ contains a from-scratch
JS port of the same fitting pipeline, runnable locally in Node or dispatched
over the DCP network -- see js/ARCHITECTURE.md for the full design, and
js/index.html for a browser-based single-galaxy job launcher.

js/tools/run_both.js compares the Fortran pipeline against the JS pipeline on
the same galaxy/seed, either against a real DCP dispatch or a local (no
network, no credentials) JS run, and writes both a JSON report and an
auto-opening HTML report (js/tools/run_both_report.html) summarizing
numerical agreement and performance. Run it from js/tools/.

Example output, no need to clone and rerun the pipeline to see what it looks
like: [js/tools/examples/fortran-vs-js-dcp-fidelity-performance.pdf](js/tools/examples/fortran-vs-js-dcp-fidelity-performance.pdf)

Fortran vs. DCP (with optional --skip-fortran flag; remove it to also run
and compare against the native Fortran leg):  
```
node run_both.js \
  --objName WALLABY_J100336-262923 \
  --mask ../TestData/WALLABY_Test_sources/WALLABY_J100336-262923/SoFiA_J100336-262923_mask.fits \
  --pa 81.2713489724864 \
  --inc 31.490766615048877 \
  --seed 42 \
  --nBootstraps 1000 \
  --cloudDensity 20 \
  --apiKey $DCP_API_KEY \
  --computeGroups google,95rhwgha \
  --slicePrice 2.000 \
  --json 1000_bootstraps_google.json \
  --skip-fortran
```

Fortran vs. Local JS (no DCP network, no credentials, no compute spend):  
```
node run_both.js \
  --objName WALLABY_J100336-262923 \
  --mask ../TestData/WALLABY_Test_sources/WALLABY_J100336-262923/SoFiA_J100336-262923_mask.fits \
  --pa 81.2713489724864 \
  --inc 31.490766615048877 \
  --seed 42 \
  --nBootstraps 10 \
  --cloudDensity 20 \
  --local \
  --json 10_bootstraps_local.json
```

---
Known Bugs
1) There is a known issue with linux vs Mac architechtures.  To adjust the code for linux:  
    a) Change line 161 in src/Inputs/DataCubeInput.f from "integer\*8 nInts(3)" to "integer nInts(3)"
    b) Change line 118 in src/Outputs/DataCubeOutputs.f from "integer\*8 naxesT(3),naxisT" to "integer naxesT(3),naxisT" (note: this line already reads "integer naxesT(3),naxisT" with no \*8 in the current source -- this step may already be applied, or may no longer be needed; worth confirming before relying on it)
    
2) There is a potential incompatibility with the third_party/colourspace package and the current version of matplotlib.


