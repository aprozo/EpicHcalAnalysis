# EpicHcalAnalysis

## positionResolutionScan

### Run Batch Farm

use `./run.sh` on SDCC shell

It will run multiple jobs with angles specified - `tileMap.txt`, the execution script is defined in `runSimBatch.sh`.

The analysis macro for filling histograms + fitting is defined in `readHCalRecoReader.C`. This macro is executed on EICRecon files.

If Simulation files and EICRecon files are already present for these angles, the script `./run.sh` will executy ONLY `readHCalRecoReader.C`

### Run local analysis

After the jobs have finished one needs to run `hadd` command inside the `eic-shell`.

```bash
cd /gpfs02/eic/alpro/output/positionResolutionScan/
hadd -f -k -j output.root Phi*Theta*/outputFill.root
```
Then copy `output.root` to the working dir. One can further analyze it and plot 2D position resolution histogram using `plot2DResolution.C`



## mapAngleTiles

The macro `calculateMap.C` is used for calculating angles for a sector `tileMap.txt`  which is further passed as an input to `positionResolutionScan`.

## calculateTileSizes

This Macro `zToEtaEpicHcal.cxx` is used for matching of STAR EEMC 24 layers to EPIC nHCAL 10 layers and also calculating all tile dimensions as well as eta ranges.







