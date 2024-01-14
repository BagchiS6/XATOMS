To run from anywhere, copy `XATOMS` to your `conda` environment's site packages (e.g., `<Your_ENV>/bin/python/3.11/site-packages/`). It assumes `<Your_env>` at least has `lammps` and `ovito` paths intalled (apart from regular packages 
like `numpy`, `json` etc. 

Within any local directory,  `driver.py`, `input_paramters.json` and `md_metadata` should be sufficient to submit the runs with `srun(mprirun) -n <num_procs> python driver.py input_parameters.json`
