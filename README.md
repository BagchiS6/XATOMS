To run globally, in new `conda` environment install dependecies from `requirements.txt`. 
(apart from standard python packages, you have to compile `lammps` as a shared library with MANYBODY and REPLICA packages included. Install ovito-python package as per the conda install guidelines from https://www.ovito.org/docs/current/python/introduction/installation.html)

To install `XATOMS` then do
`python setup install`
if you would like to develop and contribute to a different brach,
`python setup develop`
`git branch <your_branch>`
`git checkout`

Within any local directory,  `driver.py`, `input_paramters.json` and `md_metadata` (or your preferred dir-name) should be sufficient to submit the runs with `srun(mprirun) -n <num_procs> python driver.py input_parameters.json`
