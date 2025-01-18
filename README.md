# GRIT project directory

## Contents
* `./activate`:
  * Script to activate the environment for this project.
  * Usage: `source ./activate`

* `data/`
  * Any input data not created or modified by this project.

* `docs/`
  * Documentation and descriptions of the project results.

* `results/`
  * Any result/data created or modified by this project organised in
    subdirectories representing result components.

* `src/`
  * Code developed to create the results.

* `bin/`
  * Any installed software, e.g. with `--prefix=$PROJECT_ROOT/bin`

* see Wilson et al. 2017 (https://doi.org/10.1371/journal.pcbi.1005510) for details


## Version control
* All code and config files created by humans are committed to git.
* Documentation (also in binary formats) is also committed to git.
* To clone the repository with submodules use:
  `git clone --recurse-submodules <url>`
* Changes to result files should be mainly represented by the source
  code, config and build files that created them.


## Python and R environments
* This project uses Anaconda as a python and R package manager and python 3.6+.
* Create the environment using (last directory of prefix is name):
  ```
  conda env create -f src/conda_environment_base.yml -p bin/conda/ef-base
  ```

## Workflows
* All results are created by [snakemake](https://snakemake.readthedocs.io)
  workflows defined in the respective `Snakefile`.
* Executing a workflow locally on 10 CPUs:
  ```
  snakemake -j 10 <rule or file>
  ```
* Install snakemake slurm interface:
  ```
  pip install cookiecutter
  cookiecutter -f --no-input -o $PROJECT_ROOT/src/snakemake_slurm_profiles https://github.com/Snakemake-Profiles/slurm.git
  ```

## Config
* Paths: in Snakefiles
* Constants
* Research assumptions/decisions
* variable names/labels/units
* plot styles/colours