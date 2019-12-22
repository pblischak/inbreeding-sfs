# `data/`

### `cabbage/`

 - `cabbage.fs`: The main data file containing the SFS for 45 individuals of domesticated cabbage.
 - `boot/`: Folder containing 100 bootstrapped replicates of the SFS for calculating parameter uncertainties.
 - `run_cabbage_fits_folded_3epoch_F.py` and `run_cabbage_fits_folded_3epoch_noF.py`: Scripts for fitting
   demographic models for cabbage with or without inbreeding, respectively.
 - `run_cabbage_godambe_3epoch_F.py` and `run_cabbage_godambe_3epoch_noF.py`: Scripts for calculating parameter
   uncertainties for the models with and without inbreeding, respectively. These scripts use the bootstrapped
   frequency spectra in the `boot/` folder.
 - `plot_cabbage_residuals.py`: Script for generating residual plots comparing the observed and expected SFS for
   models with and without inbreeding.

### `puma/`

 - `puma.fs`: The main data file containing the 2D-SFS for five Texas puma individuals and two Florida panther individuals.
 - `boot/`: Folder containing 100 bootstrapped replicates of the SFS for calculating parameter uncertainties.
 - `run_puma_fits_folded_F.py` and `run_puma_fits_folded_noF.py`: Scripts for fitting
   demographic models for American pumas with or without inbreeding, respectively.
 - `run_puma_godambe_F.py` and `run_cabbage_godambe_noF.py`: Scripts for calculating parameter
   uncertainties for the models with and without inbreeding, respectively. These scripts use the bootstrapped
   frequency spectra in the `boot/` folder.
 - `plot_puma_residuals.py`: Script for generating residual plots comparing the observed and expected SFS for
   models with and without inbreeding.
