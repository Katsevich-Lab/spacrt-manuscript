
# Computationally efficient and statistically accurate conditional independence testing with spaCRT

This repository reproduces the results reported in arXiv version 2 the
following paper:

Z. Niu, J. Ray Choudhury, E. Katsevich. “Computationally efficient and
statistically accurate conditional independence testing with spaCRT.”
([arXiv](https://arxiv.org/pdf/2407.08911))

# Get started

First, clone the `spacrt-manuscript` repository onto your machine.

    git clone git@github.com:Katsevich-Lab/spacrt-manuscript.git

One can choose to either run the simulation or real data analysis and
obtain the figures, or directly download the results from Dropbox and
use our plotting code to reproduce the figures. We will present these
two routes separately.

# Download results data and create the figures

The data are stored in .rds format. Download the simulation results and
real data results from: [Dropbox simulation results
repository](https://www.dropbox.com/scl/fo/9f4762h4m5uyk5yqzmfg0/APfgWvyjGqaFPQF50TjRnEk?rlkey=u9239c7649y3iifzlacl7c0uy&dl=0)
and [Dropbox real data results
repository](https://www.dropbox.com/scl/fo/0c65ld2wi580828n6vfuz/AAdiZm4wBzRQ4fagb7DuaMY?rlkey=h1shu4eywp2zfa3j5i0pe9qu2&dl=0),
respectively. The following command could be used for reproducing the
plots for simulation and real data analysis respectively.

## Create the figures for real data analysis

One needs to change the `data_dir` in `realdata-code/plotting-code.R` to
the right directory where the downloaded results are. The value for
`max_cutoff` should be chosen to 100.

    Rscript realdata-code/plotting-code.R $max_cutoff

## Create the figures for simulation results

One could use the following code for reproducing the plots for
simulation results. Note the `path_rds` variable in these Rscripts
should be the path to the downloaded simulation results.

    Rscript -e 'source("simulation-code/plotting-code/assemble-plots-NB-disp-5e-2.R")'
    Rscript -e 'source("simulation-code/plotting-code/assemble-plots-NB-disp-1.R")'
    Rscript -e 'source("simulation-code/plotting-code/assemble-plots-NB-disp-10.R")'

If you would like to rerun the simulations from scratch, do not download
the results and instead follow the steps in the next section.

# Reproduce the results and figures for simulation and real data analysis

One needs to first download the `spacrt` package from
[Katsevich-lab](https://github.com/Katsevich-Lab/spacrt) using the
following R code.

    library(devtools)
    install_github("katsevich-lab/spacrt")

We used a config file to increase the portability of our code across
machines. Create a config file called `.research_config` in your home
directory.

    cd
    touch ~/.research_config

Define the following variable within this file:

- `LOCAL_SPACRT_DATA_DIR`: the location of the directory in which to
  store results.

The contents of the `.research_config` file should look like something
along the following lines.

    LOCAL_INTERNAL_DATA_DIR="/Users/ziangniu/Documents/Projects/HPCC/data/projects/"
    LOCAL_SPACRT_DATA_DIR=$LOCAL_INTERNAL_DATA_DIR"spacrt/"

Navigate to the spacrt-manuscript directory. All scripts below must be
executed from this directory. Figures will be automatically created if
one uses the following code to reproduce the results.

## Run simulation and create figures

Also, for the commands below, depending on the limits of your cluster,
you may need to set the max_gb and max_hours parameters differently. The
choice in `run_all_simulation.sh` is 16 and 4, respectively.

    qsub run_all_simulation.sh

## Run real data analysis and create figures

One can use the following command to reproduce the real data analysis
results.

    qsub run_all_realdata.sh

Table 3 in the paper can be created using
`realdata-code/sparsity_dataset.R`.
