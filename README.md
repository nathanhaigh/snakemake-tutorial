For the purpose of the webinar we will be working on the head node of an HPC cluster running Slurm. This is the most likely infrastructure that fellow bioinformaticians already find themselves using on a regular basis.

The execution of the Snakemake workflow will actually take place on the cluster head node with jobs being submitted to Slurm for queing and processing. From the head node, Snakemake will monitor the submitted jobs for their completion status and submit new jobs as dependent jobs complete sucessfully.

# Snakemake Installation

See also: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

The recommended installation route for Snakemake is through a conda environemnt. As such, you need Anaconda3, usually avaiable to you on your cluster via the module system.

```bash
# Load the Anaconda3 module on your cluster
# If it's unavailable contact the cluster sysadmin
module load \
  Anaconda3

# Create a new conda environment for snakemake v5.4.0 (currently the latest version)
# https://anaconda.org/search?q=snakemake
conda create \
  --name snakemake \
  --channel bioconda --channel conda-forge \
  --yes \
  snakemake=5.4.0
```

# Snakemake Usage

First, load the conda environment you set up:

```bash
# Load the Anaconda3 module on your cluster
# If it's unavailable contact the cluster sysadmin
module load \
  Anaconda3

# activate the conda environment you created above
source activate snakemake
```

Clone the webinar repository

```bash
# Clone the repo
git clone https://github.com/UofABioinformaticsHub/2019_EMBL-ABR_Snakemake_webinar

cd 2019_EMBL-ABR_Snakemake_webinar/snakemake-tutorial
```

Now you are ready to run Snakemake.

Without any arguments, Snakemake does not know how to submit jobs to Slurm or how to monitor them etc. It will just execute the jobs locally (i.e. on the head node itself). Snakemake requires a "profile" so that it knows how to interact with a compute backend such as Slurm. I have configured such a profile for Slurm and you can instruct Snakemake to use it with the `--profile` argument:

```bash
snakemake \
  --profile slurm \
  --cluster-config cluster-configs/default.yaml \
  --cluster-config cluster-configs/phoenix.yaml \
  --use-conda \
  --dryrun
```

Lets visalise the DAG of jobs:

```bash
snakemake \
  --dag \
  | dot -Tsvg \
  > dag.svg
```

Now lets run the workflow for real:

```bash
snakemake \
  --profile slurm \
  --cluster-config cluster-configs/default.yaml \
  --cluster-config cluster-configs/phoenix.yaml \
  --use-conda
```

# Cleanup the Enviroment

Once we have finished, we will deactive the conda environment:

```bash
source deactivate
``` 
