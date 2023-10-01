Size and power simulations for paper *Spectral Backtests Unbounded and Folded*. 

Note:

* The main scripts to run are *sizepower1.R*, *sizepower2.R*, *sizepower_tlsf.R*, and *sizepower1v.R*.
* Scripts *simSetup.R*, *DefineFmodels.R*, *DefineKernels.R*, and *DefineVtransforms.R* contain code that is shared across the main scripts, in part to ensure consistency. 
* Results are written to Folder *output*. Final publication-ready outputs are transferred manually to folder *output_publication* to protect from overwriting. 
* In place of the `choose_dist()` function, we now simply pass a named list of DGP functions. We pass v-transforms into the `doOne()` function in much the same way.

